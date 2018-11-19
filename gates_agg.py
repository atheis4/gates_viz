import os
from datetime import datetime
from glob import glob
from multiprocessing import Pool
import pandas as pd
import shutil
import sys
import uuid

try:
    # Python3
    from itertools import zip_longest
except ImportError:
    # fallback to Python2
    from itertools import izip_longest as zip_longest

from adding_machine import agg_engine as ae
from adding_machine.checks import ModelChecker

from db_queries import get_location_metadata
from hierarchies import dbtrees
from transmogrifier import super_gopher


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    args = [iter(iterable)] * n
    return list(zip_longest(fillvalue=fillvalue, *args))


def pretty_now():
    return datetime.now().strftime('[%m/%d/%Y %H:%M:%S]')


def get_pop():
    indir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
             'splits/01_populations')
    filename = 'age_pop.h5'
    return pd.read_hdf(os.path.join(indir, filename), 'data')


def get_subpop(population_df, location_ids=None, year_ids=None,
               age_group_ids=None, sex_ids=None):
    query_str = []
    if location_ids is not None:
        if not hasattr(location_ids, "__iter__"):
            location_ids = [location_ids]
        location_ids = [str(l) for l in location_ids]
        query_str.append('(location_id in [{}])'.format(
            ",".join(location_ids)))
    if age_group_ids is not None:
        if not hasattr(age_group_ids, "__iter__"):
            age_group_ids = [age_group_ids]
        age_group_ids = [str(l) for l in age_group_ids]
        query_str.append('(age_group_id in [{}])'.format(
            ",".join(age_group_ids)))
    if year_ids is not None:
        if not hasattr(year_ids, "__iter__"):
            year_ids = [year_ids]
        year_ids = [str(l) for l in year_ids]
        query_str.append('(year_id in [{}])'.format(",".join(year_ids)))
    if sex_ids is not None:
        if not hasattr(sex_ids, "__iter__"):
            sex_ids = [sex_ids]
        sex_ids = [str(l) for l in sex_ids]
        query_str.append('(sex_id in [{}])'.format(",".join(sex_ids)))
    if len(query_str) > 0:
        query_str = " & ".join(query_str)
        return population_df.query(query_str)
    else:
        return population_df


def get_me_ids(drawdir):
    pattern = '*.csv'
    files = glob(os.path.join(drawdir, pattern))
    me_list = [_file[91:-14] for _file in files]
    return [int(me) for me in me_list]


def aggregate_all_locations_mp(
        drawdir,
        stagedir,
        location_set_id,
        index_cols,
        years=[1990, 1995, 2000, 2005, 2010, 2016],
        sexes=[1, 2],
        include_leaves=False,
        operator='wtd_sum',
        custom_file_pattern=None,
        h5_tablename=None,
        single_file=True,
        gbd_round=2016,
        square_data=False,
        expand=False):

    if custom_file_pattern:
        specs_to_try = [{
            'file_pattern': custom_file_pattern,
            'h5_tablename': h5_tablename}]
    else:
        specs_to_try = super_gopher.known_specs[0:2]
    sgs = []
    for s in specs_to_try:
        try:
            sgs.append(super_gopher.SuperGopher(s, drawdir))
        except Exception:
            pass
    assert len(sgs) > 0, (
        "{} Could not find files matching known file specs in {}"
        .format(pretty_now(), drawdir))

    # Call get_pop once to cache result
    print('{} Caching population'.format(pretty_now()))
    global pop
    pop = get_pop()
    print('{} Population cached'.format(pretty_now()))

    if location_set_id == 35:
        most_detailed_ages = None

        lsvid = dbtrees.get_location_set_version_id(location_set_id,
                                                    gbd_round=gbd_round)
        lt = dbtrees.loctree(lsvid)
        depth = lt.max_depth() - 1
        pool = Pool(10)

        outfile = "{}/all_draws.h5".format(stagedir)

        while depth >= 0:
            maplist = []

            for parent_loc in lt.level_n_descendants(depth):
                this_include_leaves = include_leaves
                if len(parent_loc.children) > 0:
                    for year in years:
                        for sex in sexes:
                            maplist.append(
                                ((
                                    drawdir, lt, parent_loc.id, year, sex,
                                    index_cols, sgs, pop, this_include_leaves,
                                    operator, gbd_round, square_data,
                                    most_detailed_ages, expand),
                                    {'single_file': True}))
            res = pool.map(aclw_mp, maplist)
            if len([r for r in res if isinstance(r, tuple)]) > 0:
                print(
                    '{} !!!! SAVE RESULTS FAILED !!!! Looks like there were '
                    "errors in aggregation. See the 'uh ohs' above"
                    .format(pretty_now()))
                sys.exit()

            res = pd.concat(res)
            for col in [
                    'measure_id', 'location_id', 'year_id', 'age_group_id',
                    'sex_id']:
                res[col] = res[col].astype(int)
            print('{} Writing depth: {}'.format(pretty_now(), depth))
            if single_file:
                if not os.path.isfile(outfile):
                    res.to_hdf(
                        outfile, 'draws', mode='w', format='table',
                        data_columns=[
                            'measure_id', 'location_id', 'year_id',
                            'age_group_id', 'sex_id'])
                    sgs.append(super_gopher.SuperGopher({
                        'file_pattern': 'all_draws.h5',
                        'h5_tablename': 'draws'}, stagedir))
                else:
                    hdfs = pd.HDFStore(outfile)
                    hdfs.append('draws', res)
                    hdfs.close()
            for thissg in sgs:
                thissg._refresh_all_files()
            depth = depth - 1
        pool.close()
        pool.join()
        return lsvid

    else:
        most_detailed_ages = None

        lsvid = dbtrees.get_location_set_version_id(location_set_id,
                                                    gbd_round=gbd_round)
        lt = dbtrees.loctree(lsvid, return_many=True)

        for tree in lt:
            depth = tree.max_depth() - 1
            pool = Pool(10)

            outfile = "{}/all_draws.h5".format(stagedir)

            while depth >= 0:
                maplist = []

                for parent_loc in tree.level_n_descendants(depth):
                    this_include_leaves = include_leaves
                    if len(parent_loc.children) > 0:
                        for year in years:
                            for sex in sexes:
                                maplist.append(
                                    ((
                                        drawdir, tree, parent_loc.id, year,
                                        sex, index_cols, sgs, pop,
                                        this_include_leaves,
                                        operator, gbd_round, square_data,
                                        most_detailed_ages, expand),
                                        {'single_file': True}))
                res = pool.map(aclw_mp, maplist)
                if len([r for r in res if isinstance(r, tuple)]) > 0:
                    print(
                        '{} !!!! SAVE RESULTS FAILED !!!! Looks like there '
                        "were errors in aggregation. See the 'uh ohs' above"
                        .format(pretty_now()))
                    sys.exit()

                res = pd.concat(res)
                res = res[res.location_id.isin(
                    [44634, 44635, 44636, 44637, 44639])]
                for col in [
                        'measure_id', 'location_id', 'year_id', 'age_group_id',
                        'sex_id']:
                    res[col] = res[col].astype(int)
                print('{} Writing depth: {}'.format(pretty_now(), depth))
                if single_file:
                    if not os.path.isfile(outfile):
                        res.to_hdf(
                            outfile, 'draws', mode='w', format='table',
                            data_columns=[
                                'measure_id', 'location_id', 'year_id',
                                'age_group_id', 'sex_id'])
                        sgs.append(super_gopher.SuperGopher({
                            'file_pattern': 'all_draws.h5',
                            'h5_tablename': 'draws'}, stagedir))
                    else:
                        hdfs = pd.HDFStore(outfile)
                        hdfs.append('draws', res)
                        hdfs.close()
                for thissg in sgs:
                    thissg._refresh_all_files()
                depth = depth - 1
            pool.close()
            pool.join()
        return lsvid


def aclw_mp(mapargs):
    try:
        args, kwargs = mapargs
        (drawdir, lt, parent_loc, year, sex, index_cols,
         sgs, pop, include_leaves, operator, gbd_round, square_data,
         most_detailed_ages, expand) = args
        aggdf = aggregate_child_locs(drawdir,
                                     lt,
                                     parent_loc,
                                     year,
                                     sex,
                                     index_cols,
                                     sgs,
                                     pop,
                                     include_leaves=include_leaves,
                                     operator=operator,
                                     gbd_round=gbd_round,
                                     square_data=square_data,
                                     most_detailed_ages=most_detailed_ages,
                                     expand=expand,
                                     **kwargs)
        return aggdf
    except Exception as e:
        print(
            '{ts} Uh oh, something went wrong trying to aggregate '
            'location_id: {lid}, year: {y}, sex: {s}. Are child locations '
            'present? Error: {e}'.format(
                ts=pretty_now(), lid=parent_loc, y=year, s=sex, e=str(e)))
        return (500, str(e))


def aggregate_child_locs(
        drawdir,
        lt,
        parent_location_id,
        year_ids,
        sex_ids,
        index_cols,
        sgs,
        population_df,
        include_leaves=False,
        operator='wtd_sum',
        gbd_round=2016,
        square_data=False,
        expand=False,
        weight_col='pop_scaled',
        normalize='auto',
        force_lowmem=False,
        chunksize=4,
        draw_filters={},
        most_detailed_ages=None,
        **kwargs):

    sex_map = {'male': 1, 'female': 2, 1: 1, 2: 2}

    if square_data:
        assert most_detailed_ages is not None

    location_list = [
        l.id for l in lt.get_node_by_id(parent_location_id).children]
    location_list = list(location_list)

    # ensure we're passing sex_ids, not sex_names
    if not hasattr(sex_ids, "__iter__"):
        sex_ids = [sex_ids]
    sex_ids = [sex_map[id] for id in sex_ids]

    pops = get_subpop(
        population_df,
        location_ids=location_list,
        year_ids=year_ids,
        sex_ids=sex_ids)

    if force_lowmem:
        chunksize = chunksize
    else:
        chunksize = len(location_list)
    df = None
    indf = []
    for this_ll in grouper(location_list, chunksize, None):
        chunkdf = []
        this_ll = [l for l in this_ll if l is not None]
        for sg in sgs:
            # if location aggregates have been retrieved from modeler's
            # supergopher, no need to look in other SGs
            # (ie staging directory that contains new aggregates)
            if not this_ll:
                continue

            try:
                if expand:
                    thisdf = sg.content(
                        location_id=this_ll, year_id=year_ids,
                        skip_refresh=True, **draw_filters)
                    # replace sex_id = 3 with arg sex_ids
                    thisdf = expand.gbdize_sex(thisdf, sex_id=sex_ids)
                    thisdf = thisdf.loc[thisdf["sex_id"].isin(sex_ids)]
                    # replace aggregate age groups w/ appropriate most detailed
                    thisdf = expand.gbdize_age(thisdf, groups=True)
                    # check that no demographics were duplicated
                    dups = thisdf.duplicated(
                        subset=['location_id', 'year_id',
                                'age_group_id', 'sex_id',
                                'measure_id'])
                    assert ~(dups.any()), (
                        "If specifying expand=True, your draws"
                        " must contain only mutually exclusive"
                        " sex/age groups. For example, draws should"
                        " not contain age_group_ids 7 and 22,"
                        " or sexes 2 and 3")
                else:
                    thisdf = sg.content(
                        location_id=this_ll, year_id=year_ids, sex_id=sex_ids,
                        skip_refresh=True, **draw_filters)
                if 'pop_scaled' in thisdf.columns:
                    thisdf.drop('pop_scaled', axis=1, inplace=True)
                this_ll = list(
                    set(this_ll) - set(thisdf.location_id.unique()))
                chunkdf.append(thisdf)
            except Exception as e:
                print(e)
                continue
        chunkdf = pd.concat(chunkdf)
        chunkdf.reset_index(drop=True, inplace=True)

        # Check for missing location_ids in this_ll
        missing_lids = set(this_ll) - set(chunkdf.location_id.unique())
        assert not missing_lids, ("One or more of the child locations could "
                                  "not be found {}".format(missing_lids))

        for idx in ['location_id', 'year_id', 'age_group_id', 'sex_id']:
            try:
                pops[idx] = pops[idx].astype('int')
                chunkdf[idx] = chunkdf[idx].astype('int')
            except Exception as e:
                print(e)
                pass
        # square data
        if square_data and operator == 'wtd_sum':
            # for data outside most detailed age ids, merge on pop as normal
            no_square = chunkdf[
                ~chunkdf.age_group_id.isin(most_detailed_ages)]
            no_square = no_square.merge(
                pops, on=['location_id', 'year_id', 'age_group_id', 'sex_id'])
            # for data within most detailed age ids, square during merge
            to_square = chunkdf[chunkdf.age_group_id.isin(most_detailed_ages)]
            to_square = to_square.groupby(
                ['measure_id'], as_index=False).apply(
                    lambda to_square: _square_data(
                        to_square, pops, most_detailed_ages))
            # concat most detailed age and other age dfs back together
            chunkdf = pd.concat([to_square, no_square])
        else:
            chunkdf = chunkdf.merge(
                pops, on=['location_id', 'year_id', 'age_group_id', 'sex_id'])
        if df is not None:
            chunkdf = pd.concat([chunkdf, df])
        agg_cols = list(set(index_cols) - set(['location_id']))
        for ac in agg_cols:
            try:
                chunkdf[ac] = chunkdf[ac].astype('int')
            except Exception as e:
                print(e)
                pass
        draw_cols = [c for c in chunkdf if 'draw' in c]
        if include_leaves:
            indf.append(chunkdf[chunkdf.location_id.isin(
                [l.id for l in lt.leaves()])].copy())
        if operator == 'wtd_sum':
            if weight_col == 'scaling_factor':
                chunkdf = read_scalars(chunkdf, chunkdf.location_id.unique(),
                                       chunkdf.year_id.unique())
            elif weight_col == 'pop_scaled':
                thispopdf = ae.aggregate(
                    chunkdf[index_cols + ['pop_scaled']],
                    ['pop_scaled'],
                    agg_cols,
                    'sum')
            chunkdf = ae.aggregate(
                chunkdf[index_cols + draw_cols + [weight_col]],
                draw_cols,
                agg_cols,
                'wtd_sum',
                weight_col=weight_col,
                normalize=normalize)
            if weight_col == 'pop_scaled':
                chunkdf = chunkdf.merge(thispopdf)
        elif operator == 'sum':
            chunkdf = ae.aggregate(
                chunkdf[index_cols + draw_cols + ['pop_scaled']],
                draw_cols,
                agg_cols,
                'sum')
        df = chunkdf.copy()
    df['location_id'] = parent_location_id
    if include_leaves:
        indf = pd.concat(indf)
        df = df.append(indf)
    if 'pop_scaled' in df.columns:
        df.drop('pop_scaled', axis=1, inplace=True)
    return df.reset_index(drop=True)


def main():
    # instantiate SuperGopher for custom age draws
    indir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
             'splits/02_exposure_data')
    drawdir = os.path.join(indir, 'draws')

    outdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
              'splits/03_exposure_loc_aggregates')

    file_pattern = '{}_all_draws.csv'
    h5_tablename = None

    operator = 'wtd_sum'
    years = [1990, 1995, 2000, 2005, 2010, 2016]
    sexes = [1, 2]
    gbd_round = 2016
    gbd_round_id = 4
    round_map = {2016: 4, 2015: 3}

    meid_list = get_me_ids(drawdir)

    index = os.environ['SGE_TASK_ID']
    index = int(index) - 1
    meid = meid_list[index]

    file_pattern = file_pattern.format(meid)

    print('{} Beginning custom aggregation for Gates Viz modelable entity {}'
          .format(pretty_now(), meid))

    stagedir = os.path.join(outdir, '{}/draws'.format(meid))

    try:
        os.makedirs(stagedir)
    except OSError:
        pass

    # Aggregate (excluding SDI)
    lsvid = aggregate_all_locations_mp(
        drawdir, stagedir, 35,
        ['year_id', 'age_group_id', 'sex_id', 'measure_id'], years,
        sexes, True, operator, custom_file_pattern=file_pattern,
        h5_tablename=h5_tablename, gbd_round=gbd_round,
        square_data=False, expand=False)

    # Aggregate SDI
    lsvid = aggregate_all_locations_mp(
        drawdir, stagedir, 40,
        ['year_id', 'age_group_id', 'sex_id', 'measure_id'], years,
        sexes, True, operator, custom_file_pattern=file_pattern,
        h5_tablename=h5_tablename, gbd_round=gbd_round,
        square_data=False, expand=False)

    print('{} Saving custom aggregates completed'.format(pretty_now()))


if __name__ == '__main__':
    main()
