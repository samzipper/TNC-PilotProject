## MODFLOW_HTC_CheckFailure.py
# This script is intended to look at all the .out files after an HTC run
# and figure out which ones failed.

# imports
import os, glob
import numpy as np

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', 'SFR')
prefix_runs = 'mf'

# get filenames
files = glob.glob(os.path.join(dir_runs, prefix_runs+'*.out'))

# scroll through all files
num_all = []
succ_all = []
for f in np.arange(0, len(files)):
    # extract number
    num = files[f][len(os.path.join(dir_runs, prefix_runs)):-12]

    # open text file
    fo = open(files[f], "r")
    lines = fo.readlines()
    fo.close()

    # figure out if failed or success
    num_succ = False
    if "Normal termination of simulation" in lines[-1]:
        num_succ = True

    # add to overall output file
    num_all.append(num)
    succ_all.append(num_succ)

# combine and save output
output = np.column_stack((num_all,succ_all))
output_sort = output[output[:,0].argsort()]
np.savetxt(os.path.join(dir_runs, 'CheckFailure.csv'), output_sort,
        delimiter=",", fmt='%s', header="WellNum,Success", comments='')

print('successes: ', sum(succ_all))
print('fails: ', len(succ_all)-sum(succ_all))
