## postprocess_thisRun_SFR.py
# This script will open SFR output files and summarize the leakage by SegNum.
# It needs to be run within the model output folders.

import os
import numpy as np
import flopy.utils.binaryfile as bf
import flopy.utils.sfroutputfile as sf
import pandas as pd
import glob

# set modelname
timeType = 'Intermittent'  # 'SteadyState' or 'Transient' or 'Intermittent'
modelname = 'Navarro-'+timeType

## load SFR input
back_to_base_dir = os.path.join('..', '..', '..', '..', '..', '..', '..')
isfr_ReachData = pd.read_table(os.path.join(back_to_base_dir, 'modflow', 'input', 'isfr_ReachData.txt'), delimiter=' ')
isfr_SegmentData = pd.read_table(os.path.join(back_to_base_dir, 'modflow', 'input', 'isfr_SegmentData.txt'), delimiter=' ')

## load SFR binary files
sfrout = sf.SfrFile(modelname+'.sfr.out')
sfr_df = sfrout.get_dataframe()

# summarize by segment
sfr_df_summary = sfr_df.groupby(['segment', 'kstpkper'], as_index=False).agg({'Qaquifer': 'sum'})

## make a copy of iriv to avoid altering original
isfr_SegmentData_copy = isfr_SegmentData.copy()

# join with isfr_SegmentData
sfr_merge = pd.merge(isfr_SegmentData_copy, sfr_df_summary,
                     left_on=['SFR_NSEG'], right_on=['segment'])

## calculate net well pumping rate (close, but not identical, to qdes)
mnwout = bf.CellBudgetFile(modelname+'.mnw2.out', verbose=False)
pump_steps = pd.DataFrame({'kstpkper': list(set(sfr_merge['kstpkper'])), 'MNW_net': [0.0]*len(set(sfr_merge['kstpkper']))})
for k in range(0,pump_steps.shape[0]):
    mnwout_data = mnwout.get_data(kstpkper=pump_steps['kstpkper'][k], text='MNW2', full3D=False)
    pump_steps['MNW_net'][k] = sum(mnwout_data[0]['q'])
mnwout.close()

# join MNW with sfr data
sfr_merge = pd.merge(sfr_merge, pump_steps, left_on=['kstpkper'], right_on=['kstpkper'])[['SFR_NSEG', 'Qaquifer',  'kstpkper', 'MNW_net']]

## save all output
sfr_merge.round({'Qaquifer':3, 'MNW_net':3}).to_csv(modelname+'_postprocess.csv', header="SFR_NSEG,leakage,kstpkper,MNW_net", index=False)
