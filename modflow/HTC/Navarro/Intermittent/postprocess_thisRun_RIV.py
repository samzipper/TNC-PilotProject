## postprocess_thisRun_RIV.py
# This script will open RIV output files and summarize the leakage by SegNum.
# It needs to be run within the model output folders.

import os
import numpy as np
import flopy.utils.binaryfile as bf
import pandas as pd
import glob

# set modelname
timeType = 'Intermittent'  # 'SteadyState' or 'Transient' or 'Intermittent'
modelname = 'Navarro-'+timeType

## load RIV input
back_to_base_dir = os.path.join('..', '..', '..', '..', '..', '..', '..')
iriv = pd.read_table(os.path.join(back_to_base_dir, 'modflow', 'input', 'iriv.txt'), delimiter=' ')
iriv_ReachData = pd.read_table(os.path.join(back_to_base_dir, 'modflow', 'input', 'iriv_ReachData.txt'), delimiter=' ')

## load output files
rivout = bf.CellBudgetFile(modelname+'.riv.out', verbose=False)
mnwout = bf.CellBudgetFile(modelname+'.mnw2.out', verbose=False)
outputTimes = rivout.get_times()

# scroll through times
start_flag = True
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
    iriv_out['Time'] = float(time)
    
    ## calculate net well pumping rate (close, but not identical, to qdes)
    mnwout_data = mnwout.get_data(totim=time, text='MNW2', full3D=False)
    iriv_out['MNW_net'] = sum(mnwout_data[0]['q'])
    
    ## add to overall data frame
    if (start_flag):
        iriv_all = iriv_out
        start_flag = False
    else:
        iriv_all = iriv_all.append(iriv_out)

# close binary files
rivout.close()
mnwout.close()

## save all output
iriv_all.round({'leakage':3, 'Time':3, 'MNW_net':3}).to_csv(modelname+'_postprocess.csv', header="SegNum,leakage,Time,MNW_net", index=False)
