import os
import pandas as pd

from hierarchies import dbtrees
from db_queries import get_location_metadata


def melt_age_cols(df):
    id_vars = ['ihme_loc_id', 'year', 'sex']
    value_vars = [c for c in df if c.startswith('pys')]
    var_name = 'age_group'
    value_name = 'population'
    df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name=var_name,
                 value_name=value_name)
    return df


def aggregate_ages(df):
    to_agg = df[df.age_group.isin([2, 3, 4])]
    dont_agg = df[df.age_group == 1]
    id_vars = ['ihme_loc_id', 'year', 'sex_id']
    temp = to_agg.groupby(id_vars,
                          as_index=False)['population'].aggregate('sum')

    dont_agg['age_group_id'] = 238
    temp['age_group_id'] = 34

    data = pd.concat([dont_agg, temp])
    return data[['age_group_id', 'ihme_loc_id', 'population', 'sex_id',
                 'year']]


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


def main():
    popdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
              'splits/01_populations')
    popfile = 'single_age_u5_pops.dta'

    popdf = pd.read_stata(os.path.join(popdir, popfile))

    age_ranges = ['pys_{}'.format(x + 1) for x in range(4)]
    age_map = {val: idx + 1 for idx, val in enumerate(age_ranges)}

    popdf = melt_age_cols(popdf)
    popdf = popdf[popdf.age_group.isin(age_ranges)]

    popdf['age_group'] = popdf.age_group.map(age_map)
    popdf['sex_id'] = popdf.sex.map({'male': 1, 'female': 2})

    popdf = aggregate_ages(popdf)
    popdf.rename(columns={'year': 'year_id'}, inplace=True)

    locs = get_location_metadata(location_set_id=35, gbd_round_id=4)
    locs = locs[['location_id', 'ihme_loc_id']]

    popdf = pd.merge(popdf, locs, on='ihme_loc_id', how='left')
    popdf.drop(labels=['ihme_loc_id'], inplace=True, axis=1)
    popdf.rename(columns={'population': 'pop_scaled'}, inplace=True)

    lsvid = dbtrees.get_location_set_version_id(35, gbd_round=2016)
    lt = dbtrees.loctree(lsvid)

    index_cols = ['age_group_id', 'sex_id', 'year_id', 'location_id']
    data_cols = ['pop_scaled']

    aggpop = agg_hierarchy(lt, popdf, index_cols, data_cols, 'location_id')

    sdi_lsvid = dbtrees.get_location_set_version_id(40, gbd_round=2016)
    sdi_lts = dbtrees.loctree(sdi_lsvid, return_many=True)

    for tree in sdi_lts:
        sdi_agg_df = agg_hierarchy(tree,
                                   popdf,
                                   index_cols,
                                   data_cols,
                                   'location_id')
        sdi_agg_df = sdi_agg_df[sdi_agg_df.location_id.isin(
            [44634, 44635, 44636, 44637, 44639])]
        aggpop = aggpop.append(sdi_agg_df)

    outdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
              'splits/01_populations')
    outfile = os.path.join(outdir, 'age_pop.h5')

    aggpop.to_hdf(outfile, 'data', mode='w', format='table',
                  data_columns=index_cols)
    print('fin')


if __name__ == '__main__':
    main()
