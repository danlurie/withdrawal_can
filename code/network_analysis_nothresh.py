import os
import argparse
import matlab.engine
from glob import glob

OUT_DIR = '/home/despoB/dlurie/Projects/withdrawal_CAN/results/'
assert (os.path.exists(OUT_DIR)),"Output directory does not exist!"
MAT_FILE_FPT = '/home/despoB/dlurie/Projects/withdrawal_CAN/data/cFosCAN_alcohol_{0}_mat_20180526.csv'
THRESH_FILE_FPT = '/home/despoB/dlurie/Projects/withdrawal_CAN/data/cFosCAN_alcohol_withdrawal_mat_thresh{0}_CC_20180526.csv'
GAMMA = 1.0

# Initialize MATLAB interface
eng = matlab.engine.start_matlab()
eng.addpath('/home/despoB/dlurie/Projects/withdrawal_CAN/code', nargout=0)

conditions = ['withdrawal', 'air', 'naive']

for condition in conditions:
    print('Condition: {0}...'.format(condition))

    graph_out_dir = os.path.join(OUT_DIR, "{0}_nothresh".format(condition))
    if os.path.isdir(graph_out_dir) is not True:
        os.makedirs(graph_out_dir)

    mat_file_path = MAT_FILE_FPT.format(condition)
    
    # Run graph-theory pipeline
    print('...Running MATLAB graph pipeline...')
    out_prefix = os.path.join(graph_out_dir, 'gamma{0}'.format(str(GAMMA)))
    res = eng.network_analysis_nothresh(mat_file_path, out_prefix, GAMMA, nargout=0)

thresholds = ['075', '085', '095']

for threshold in thresholds:
    print('Withdrawal threshold: {0}...'.format(threshold))

    graph_out_dir = os.path.join(OUT_DIR, "withdrawal_thresh{0}".format(threshold))
    if os.path.isdir(graph_out_dir) is not True:
        os.makedirs(graph_out_dir)

    mat_file_path = THRESH_FILE_FPT.format(threshold)
    
    # Run graph-theory pipeline
    print('...Running MATLAB graph pipeline...')
    out_prefix = os.path.join(graph_out_dir, 'gamma{0}'.format(str(GAMMA)))
    res = eng.network_analysis_nothresh(mat_file_path, out_prefix, GAMMA, nargout=0)

print('...Done!')


