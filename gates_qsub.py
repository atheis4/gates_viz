import os
from glob import glob
import subprocess


def get_me_ids(drawdir):
    pattern = '*.csv'
    files = glob(os.path.join(drawdir, pattern))
    me_list = [_file[91:-14] for _file in files]
    return [int(me) for me in me_list]


if __name__ == '__main__':
    indir = ('/FILEPATH_TO/Child Growth Failure/Gates_CGF_Viz/custom_age_'
             'splits/02_exposure_data')
    drawdir = os.path.join(indir, 'draws')
    meid_list = get_me_ids(drawdir)
    limit = len(meid_list)

    dirname = os.path.dirname(os.path.abspath(__file__))
    py_file = 'gates_agg.py'
    py_file = os.path.join(dirname, py_file)

    jname = 'gates_viz_loc_aggregation'
    sh_script = '/FILEPATH_TO/atheis/utils/py_shell.sh'

    call = ('qsub -t 1:{limit} -cwd -pe multi_slot 50 -now no '
            '-e /FILEPATH_TO/atheis/errors '
            '-o /FILEPATH_TO/atheis/output '
            '-N {name} {shell} {python_file}'.format(limit=limit,
                                                     name=jname,
                                                     shell=sh_script,
                                                     python_file=py_file))

    subprocess.call(call, shell=True)
