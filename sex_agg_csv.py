import os
from glob import glob
import numpy as np
import pandas as pd

from adding_machine import agg_engine as ae


def get_pop():
    indir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
             'splits/01_populations')
    filename = 'age_pop.h5'
    return pd.read_hdf(os.path.join(indir, filename), 'data')


def agg_sexes(df, pop):
    sex_agg = df.copy()
    index_cols = ['location_id', 'age_group_id', 'year_id']
    draw_cols = [c for c in sex_agg if c.startswith('draw')]
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


def create_estimates(df):
    thisdf = df.copy()
    draw_cols = [c for c in thisdf if c.startswith('draw')]
    thisdf['mean'] = np.mean(thisdf[draw_cols], axis=1)
    thisdf['lower_ui'] = np.percentile(thisdf[draw_cols], q=2.5, axis=1)
    thisdf['upper_ui'] = np.percentile(thisdf[draw_cols], q=97.5, axis=1)
    return thisdf[['age_group_id', 'location_id', 'sex_id', 'year_id',
                   'modelable_entity_id', 'mean', 'upper_ui', 'lower_ui']]


def main():
    pop = get_pop()

    project_dir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_'
                   'age_splits')
    drawdirs = os.path.join(project_dir, '03_exposure_loc_aggregates/{}/draws')
    globpattern = drawdirs.format('*')

    me_dirs = glob(globpattern)
    for me in me_dirs:
        me_id = me.split('/')
        me_id = int(me_id[7])
        print('Converting ME: {} HDF to CSV'.format(me_id))

        out = os.path.join(project_dir,
                           '03_exposure_loc_aggregates/{}/estimates')
        outfile = out.format(me_id)

        df = pd.read_hdf(os.path.join(me, 'all_draws.h5'), 'draws')
        df['modelabel_entity_id'] = me_id

        outfile = os.path.join(outfile, 'estimates.csv')
        sex_agg = agg_sexes(df, pop)
        df = df.append(sex_agg)
        df.reset_index(inplace=True, drop=True)
        df = create_estimates(df)
        df.to_csv(outfile, index=None)


if __name__ == '__main__':
    main()
