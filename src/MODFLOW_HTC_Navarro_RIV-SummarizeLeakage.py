## MODFLOW_HTC_Navarro_RIV-SummarizeLeakage.py
# This script will open RIV output files and summarize the leakage by SegNum.
# Needs the script MODFLOW_HTC_Navarro_CheckFailures.py to be run first.

import os
import numpy as np
import flopy.utils.binaryfile as bf
import pandas as pd
import glob

## set things up
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', 'RIV')
prefix_runs = 'mf'
modelname = 'Navarro-SteadyState'
time=1

## load RIV input
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')
iriv_ReachData = pd.read_table(os.path.join('modflow', 'input', 'iriv_ReachData.txt'), delimiter=' ')

## figure out which runs succeeded (output from CheckFailures script)
succ = pd.read_table(os.path.join(dir_runs, 'CheckFailure.csv'), delimiter=",")

## loop through and make output data frame
start_flag = True
for w in succ.WellNum:
    # check if model ran successfully
    converge = succ.Success.loc[(succ['WellNum']==w)].bool()
    if converge:
        ## make a copy of iriv to avoid altering original
        iriv_copy = iriv.copy()

        ## load RIV output
        rivout = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.riv.out'), verbose=False)
        rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
        iriv_copy['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
        rivout.close()

        ## join leakage to reach_data
        iriv_merge = pd.merge(iriv_ReachData[['SegNum', 'row', 'col','seg_proportion']], iriv_copy[['row', 'col', 'leakage']],
                              how='left', on=['row','col'])
                              
        ## summarize by segment number
        iriv_merge['leakage'] = iriv_merge['leakage']*iriv_merge['seg_proportion']
        iriv_out = iriv_merge.groupby('SegNum', as_index=False).agg({'leakage': 'sum'})
        iriv_out['WellNum'] = w

        ## save RIV output
        np.savetxt(os.path.join(dir_runs, prefix_runs+str(w), 'RIV_leakage.txt'), iriv_out,
                header="SegNum,leakage,WellNum", comments='')
        
        ## add to overall data frame
        if (start_flag):
            iriv_all = iriv_out
            start_flag = False
        else:
            iriv_all = iriv_all.append(iriv_out)

    ## status update
    print(w, ' complete')

## save all output
np.savetxt(os.path.join(dir_runs, 'RIV-SummarizeLeakage.csv'), iriv_all,
        delimiter=",", header="SegNum,leakage,WellNum", comments='')

## print total leakage across all segments
iriv_summary = iriv_all.groupby('WellNum', as_index=False).agg({'leakage': 'sum'})
print(iriv_summary.head())
