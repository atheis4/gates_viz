import os
from datetime import datetime
import pandas as pd

from adding_machine import agg_engine as ae
from db_queries import get_location_metadata
from hierarchies.tree import parent_child_to_tree


def pretty_now():
    return datetime.now().strftime("[%m/%d/%Y %H:%M:%S]")


def melt_age_cols(df):
    id_cols = ['location_id', 'year_id', 'sex_id']
    value_cols = [c for c in df if c.startswith('mean')]
    var_name = 'age_group'
    value_name = 'population'
    df = df.melt(id_vars=id_cols, value_vars=value_cols, var_name=var_name,
                 value_name=value_name)
    return df


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


def agg_sexes(df):
    sex_agg = df.copy()
    index_cols = ['location_id', 'age_group_id', 'year_id']
    data_cols = ['pop_scaled']
    sex_agg = ae.aggregate(
        sex_agg[index_cols + data_cols],
        data_cols,
        index_cols,
        'sum')
    sex_agg['sex_id'] = 3
    return sex_agg


def main():
    popdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
              'splits/01_populations')
    popfile = 'forecast_under_5_pops.csv'

    print('{} read in raw pop data from csv'.format(pretty_now()))
    popdf = pd.read_csv(os.path.join(popdir, popfile))

    age_ranges = ['mean_12_to_23', 'mean_2_to_4']
    age_map = {age_ranges[0]: 238, age_ranges[1]: 34}

    popdf = melt_age_cols(popdf)
    popdf = popdf[popdf.age_group.isin(age_ranges)]
    popdf['age_group_id'] = popdf.age_group.map(age_map)

    for col in ['location_id', 'sex_id', 'age_group_id']:
        popdf[col] = popdf[col].astype(int)

    popdf.rename(columns={'population': 'pop_scaled'}, inplace=True)

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
    data_cols = ['pop_scaled']

    # aggregate up standard custom tree
    print('{} aggregate pop from custom tree'.format(pretty_now()))
    aggpop = agg_hierarchy(custom_tree, popdf, index_cols, data_cols,
                           'location_id')

    # SDI locations
    sdi_locs = get_location_metadata(location_set_id=40, gbd_round_id=4)
    sdi_locs = sdi_locs[['location_id', 'parent_id', 'location_name', 'level',
                         'location_name_short', 'map_id', 'location_type',
                         'is_estimate']]
    sdi_ids = [44635, 44634, 44639, 44636, 44637]

    sdi_df_list = []
    for _id in sdi_ids:
        print('{} processing sdi: {}'.format(pretty_now(), _id))
        thisdf = sdi_locs[sdi_locs.parent_id == _id]
        thisdf = thisdf[thisdf.location_id.isin(locs_to_keep + [_id])]
        thistree = create_custom_tree(thisdf)
        print('{} aggregate pop from {} tree'.format(pretty_now(), _id))
        thisaggpop = agg_hierarchy(thistree, popdf, index_cols, data_cols,
                                   'location_id')
        thisaggpop = thisaggpop[thisaggpop.location_id == _id]
        aggpop = aggpop.append(thisaggpop)

    sexagg = agg_sexes(aggpop)
    aggpop = aggpop.append(sexagg)

    outdir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
              'splits/01_populations')
    outfile = os.path.join(outdir, 'future_pop.h5')

    print('{} output'.format(pretty_now()))
    aggpop.to_hdf(outfile, 'data', mode='w', format='table',
                  data_columns=index_cols)
    print('{} fin'.format(pretty_now()))


if __name__ == '__main__':
    main()
