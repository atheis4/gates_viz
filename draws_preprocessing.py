import os
from glob import glob
import pandas as pd


def set_age_id(row):
    age_end = row.age_upper
    age_group_id = 24 if age_end == 5 else 238
    return age_group_id


def main():
    """Reading in each individual draw csv, adding the appropriate age_group_id
    to the data, eliminating some unneeded columns and overwriting back into
    the original folder."""

    drawdir = ('/FILEPATH_TO/Gates_CGF_Viz/custom_age_'
               'splits/exposure_data/draws')
    filepattern = '*.csv'

    files = glob(os.path.join(drawdir, filepattern))

    age_map = {5: 34, 2: 238}

    for _file in files:
        print('processing file {}'.format(_file))
        df = pd.read_csv(_file)
        df['age_group_id'] = df.age_upper.map(age_map)
        draw_cols = [c for c in df if c.startswith('draw')]
        index_cols = ['location_id', 'measure_id', 'age_group_id', 'sex_id',
                      'year_id', 'modelable_entity_id', 'model_version_id']
        df = df[index_cols + draw_cols]
        df.to_csv(_file, index=False)
        print('finished processing {}'.format(_file))


if __name__ == '__main__':
    main()
