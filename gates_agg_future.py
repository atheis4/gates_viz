import os
from datetime import datetime
from glob import glob
import pandas as pd

from adding_machine import agg_engine as ae
from hierarchies.tree import parent_child_to_tree
from db_queries import get_location_metadata


def pretty_now():
    return datetime.now().strftime('[%m/%d/%Y %H:%M:%S]')


def parse_me_name(file_pattern):
    return file_pattern[100:-4]


def agg_hierarchy(tree, df, index_cols, data_cols, dimension):
    # keep only the most detailed
    thisdf = df.copy()
    thisdf = thisdf[index_cols + data_cols]
    leaves = [node.id for node in tree.leaves()]
    thisdf = thisdf[thisdf[dimension].isin(leaves)]

    # loop through levels from the most detailed
    group_cols = [col for col in index_cols if col != dimension]
    md = tree.max_depth()
    lvl = md - 1
    while lvl >= 0:
        aggs = []
        for node in tree.level_n_descendants(lvl):
            child_ids = [c.id for c in node.children]
            if len(child_ids) > 0:
                agg = thisdf[thisdf[dimension].isin(child_ids)]
                group_cols = [col for col in index_cols if col != dimension]
                if not group_cols:
                    group_cols = [dimension]
                    agg[dimension] = node.id
                agg = agg.groupby(group_cols).sum().reset_index()
                agg[dimension] = node.id
                aggs.append(agg)
        aggs = pd.concat(aggs)
        thisdf = pd.concat([thisdf, aggs])
        lvl = lvl - 1

    # combine
    thisdf = thisdf.groupby(index_cols).sum().reset_index()
    return thisdf


def create_locs_to_keep(loc_metadata):
    locs_to_keep = loc_metadata[loc_metadata.level <= 3].location_id.tolist()
    nonsov = loc_metadata[
        loc_metadata.location_type == 'nonsovereign'].location_id.tolist()
    locs_to_keep = list(set(locs_to_keep) - set(nonsov))
    return locs_to_keep


def create_custom_tree(custom_loc_df):
    tree = parent_child_to_tree(
        custom_loc_df,
        'parent_id',
        'location_id',
        info_cols=['location_name', 'location_name_short', 'map_id',
                   'location_type', 'is_estimate'])
    return tree


def get_pop():
    indir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
             'splits/01_populations/')
    filename = 'future_pop.h5'
    return pd.read_hdf(os.path.join(indir, filename), 'data')


def agg_sexes(df, pop):
    sex_agg = df.copy()
    index_cols = ['location_id', 'age_group_id', 'year_id']
    draw_cols = ['lower', 'mean', 'upper']
    sex_agg = sex_agg.merge(
        pop, on=['location_id', 'year_id', 'sex_id', 'age_group_id'])
    sex_agg = ae.aggregate(
        sex_agg[index_cols + draw_cols + ['pop_scaled']],
        draw_cols,
        index_cols,
        'wtd_sum',
        weight_col='pop_scaled')
    sex_agg['sex_id'] = 3
    return sex_agg


def main():
    drawdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
               'splits/02_exposure_data/final_forecast')
    filepattern = '*.csv'
    files = glob(os.path.join(drawdir, filepattern))

    me_name_to_meid = {
        'stunting_mild': 10557, 'stunting_moderate': 10556,
        'stunting_severe': 8949, 'underweight_mild': 10561,
        'underweight_moderate': 10560, 'underweight_severe': 2540,
        'wasting_mild': 10559, 'wasting_moderate': 10558,
        'wasting_severe': 8945}

    # Location Metadata
    locs = get_location_metadata(location_set_id=35, gbd_round_id=4)
    locs = locs[['location_id', 'parent_id', 'location_name', 'level',
                 'location_name_short', 'map_id', 'location_type',
                 'is_estimate']]
    # Generate locations to keep: 188 + parents
    locs_to_keep = create_locs_to_keep(locs)
    custom_loc_df = locs[locs.location_id.isin(locs_to_keep)]
    # create main custom tree
    print('{} creating custom tree'.format(pretty_now()))
    custom_tree = create_custom_tree(custom_loc_df)
    index_cols = ['age_group_id', 'sex_id', 'year_id', 'location_id']
    data_cols = ['lower', 'mean', 'upper']
    # Create SDI trees
    sdi_locs = get_location_metadata(location_set_id=40, gbd_round_id=4)
    sdi_locs = sdi_locs[['location_id', 'parent_id', 'location_name', 'level',
                         'location_name_short', 'map_id', 'location_type',
                         'is_estimate']]
    sdi_ids = [44635, 44634, 44639, 44636, 44637]
    sdi_trees = []
    for _id in sdi_ids:
        print('{} creating sdi tree for {}'.format(pretty_now(), _id))
        thisdf = sdi_locs[sdi_locs.parent_id == _id]
        thisdf = thisdf[thisdf.location_id.isin(locs_to_keep + [_id])]
        sdi_trees.append(create_custom_tree(thisdf))

    # get population
    pops = get_pop()

    for _file in files:
        # Define me_name
        me_name = parse_me_name(_file)
        meid = me_name_to_meid[me_name]
        print('{} processing file: {}'.format(pretty_now(), meid))
        df = pd.read_csv(_file)
        df.rename(columns={'worse': 'lower',
                           'reference': 'mean',
                           'better': 'upper'},
                  inplace=True)
        df = df[['location_id', 'age_group_id', 'sex_id', 'year_id', 'lower',
                 'mean', 'upper']]
        df['modelable_entity_id'] = meid
        # Remove bad locations
        bad_locs = [298, 305, 349, 351, 376, 385, 422, 433, 434, 4636, 4749]
        df = df[~df.location_id.isin(bad_locs)]
        # convert to counts
        print('{} convert to counts before aggregation'.format(pretty_now()))
        df = df.merge(pops, on=index_cols, how='left')
        for i in data_cols:
            df[i] = df[i] * df['pop_scaled']
        # aggregate all trees
        print('{} agg custom loc tree'.format(pretty_now()))
        agg_results = agg_hierarchy(custom_tree, df, index_cols, data_cols,
                                    dimension='location_id')
        for sdi_tree in sdi_trees:
            print('{} agg sdi tree for: {}'
                  .format(pretty_now(), sdi_tree.root))
            this_agg = agg_hierarchy(sdi_tree, df, index_cols, data_cols,
                                     dimension='location_id')
            this_agg = this_agg[this_agg.location_id.isin(sdi_ids)]
            agg_results = agg_results.append(this_agg)

        # copy data and set metric id to 1 for counts
        print('{} copy counts to new df'.format(pretty_now()))
        agg_counts = agg_results.copy()
        agg_counts = agg_counts[agg_counts.sex_id.isin([1, 2])]
        sex_agg = agg_sexes(agg_counts, pops)
        agg_counts = agg_counts.append(sex_agg)
        agg_counts['metric_id'] = 1
        # convert back to rate space
        print('{} converting back to rate space'.format(pretty_now()))
        agg_results = agg_results.merge(pops, on=index_cols, how='left')
        for i in data_cols:
            agg_results[i] = agg_results[i] / agg_results['pop_scaled']
        agg_results['metric_id'] = 3
        print('{} append counts to rate df'.format(pretty_now()))
        agg_results = agg_results.append(agg_counts)

        outdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_'
                  'age_splits/03_exposure_loc_aggregates/02_forecast_'
                  'prevalence/{}'.format(meid))
        outfile = '{}_prevalence_estimates.csv'.format(meid)
        print('{} saving as csv'.format(pretty_now()))
        agg_results.to_csv(os.path.join(outdir, outfile), index=False)
        print('{} finished processing meid: {}'.format(pretty_now(), meid))


if __name__ == '__main__':
    main()
