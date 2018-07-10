## MODFLOW_HTC_Navarro_RIV-SummarizeLeakage.py
# This script will open RIV output files and summarize the leakage by SegNum.
# Needs the script MODFLOW_HTC_Navarro_CheckFailures.py to be run first.

import os
import numpy as np
import flopy.utils.binaryfile as bf
import pandas as pd
import glob

# what stream BC and modflow version?
stream_BC = 'RIV'  # 'RIV' or 'SFR'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
timeType = 'Intermittent'  # 'SteadyState' or 'Transient' or 'Intermittent'

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', timeType, stream_BC, modflow_v)
prefix_runs = 'mf'
modelname = 'Navarro-'+timeType

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

        ## load output files
        rivout = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.riv.out'), verbose=False)
        mnwout = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.mnw2.out'), verbose=False)
        outputTimes = rivout.get_times()

        # scroll through times
        for time in outputTimes:
            ## make a copy of iriv to avoid altering original
            iriv_copy = iriv.copy()
            
            # get riv output data
            rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
            iriv_copy['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
            
            ## join leakage to reach_data
            iriv_merge = pd.merge(iriv_ReachData[['SegNum', 'row', 'col','seg_proportion']], iriv_copy[['row', 'col', 'leakage']],
                    how='left', on=['row','col'])
            
            ## summarize by segment number
            iriv_merge['leakage'] = iriv_merge['leakage']*iriv_merge['seg_proportion']        
            iriv_out = iriv_merge.groupby('SegNum', as_index=False).agg({'leakage': 'sum'})
            iriv_out['WellNum'] = w
            iriv_out['Time'] = time
            
            ## calculate net well pumping rate (close, but not identical, to qdes)
            mnwout_data = mnwout.get_data(totim=time, text='MNW2', full3D=False)
            iriv_out['MNW_net'] = sum(mnwout_data[0]['q'])
            
            ## add to overall data frame
            if (start_flag):
                iriv_all = iriv_out
                start_flag = False
            else:
                iriv_all = iriv_all.append(iriv_out)
                
            ## status update
            print('Well '+str(w)+' Time '+str(time)+ ' complete')

        # close binary files
        rivout.close()
        mnwout.close()

## save all output
iriv_all.to_csv(os.path.join(dir_runs, 'RIV-SummarizeLeakage.csv'), header="SegNum,leakage,WellNum,Time,MNW_net", index=False)
