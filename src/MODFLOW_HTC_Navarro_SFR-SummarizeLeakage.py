## MODFLOW_HTC_Navarro_SFR-SummarizeLeakage.py
# This script will open SFR output files and summarize the leakage by SegNum.
# Needs the script MODFLOW_HTC_Navarro_CheckFailures.py to be run first.

import os
import numpy as np
import flopy.utils.binaryfile as bf
import flopy.utils.sfroutputfile as sf
import pandas as pd
import glob

# what stream BC and modflow version?
stream_BC = 'SFR'  # 'RIV' or 'SFR'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
timeType = 'Transient'  # 'SteadyState' or 'Transient'

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', timeType, stream_BC, modflow_v)
prefix_runs = 'mf'
modelname = 'Navarro-'+timeType

## load RIV input
isfr_ReachData = pd.read_table(os.path.join('modflow', 'input', 'isfr_ReachData.txt'), delimiter=' ')
isfr_SegmentData = pd.read_table(os.path.join('modflow', 'input', 'isfr_SegmentData.txt'), delimiter=' ')

## figure out which runs succeeded (output from CheckFailures script)
succ = pd.read_table(os.path.join(dir_runs, 'CheckFailure.csv'), delimiter=",")

## loop through and make output data frame
start_flag = True
for w in succ.WellNum:
    # check if model ran successfully
    converge = succ.Success.loc[(succ['WellNum']==w)].bool()
    
    if converge:
        ## load SFR binary files
        sfrout = sf.SfrFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.sfr.out'))
        sfr_df = sfrout.get_dataframe()
        
        # summarize by segment
        sfr_df_summary = sfr_df.groupby(['segment', 'kstpkper'], as_index=False).agg({'Qaquifer': 'sum'})

        ## make a copy of iriv to avoid altering original
        isfr_SegmentData_copy = isfr_SegmentData.copy()

        # join with isfr_SegmentData
        sfr_merge = pd.merge(isfr_SegmentData_copy, sfr_df_summary,
                             left_on=['SFR_NSEG'], right_on=['segment'])
                             
        # add well number column
        sfr_merge['WellNum'] = w
        
        ## calculate net well pumping rate (close, but not identical, to qdes)
        mnwout = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.mnw2.out'), verbose=False)
        pump_steps = pd.DataFrame({'kstpkper': list(set(sfr_merge['kstpkper'])), 'MNW_net': [0.0]*len(set(sfr_merge['kstpkper']))})
        for k in range(0,pump_steps.shape[0]):
            mnwout_data = mnwout.get_data(kstpkper=pump_steps['kstpkper'][k], text='MNW2', full3D=False)
            pump_steps['MNW_net'][k] = sum(mnwout_data[0]['q'])
            print('Well ', str(w), ' kstpkper ', str(k), ' complete')
        mnwout.close()

        # join MNW with sfr data
        sfr_merge = pd.merge(sfr_merge, pump_steps, left_on=['kstpkper'], right_on=['kstpkper'])

        ## add to overall data frame
        if (start_flag):
            sfr_all = sfr_merge[['SFR_NSEG', 'Qaquifer', 'WellNum', 'kstpkper', 'MNW_net']]
            start_flag = False
        else:
            sfr_all = sfr_all.append(sfr_merge[['SFR_NSEG', 'Qaquifer', 'WellNum', 'kstpkper',  'MNW_net']])

    ## status update
    print(w, ' complete')

## save all output
sfr_all.to_csv(os.path.join(dir_runs, 'SFR-SummarizeLeakage.csv'), header="SFR_NSEG,leakage,WellNum,kstpkper,MNW_net", index=False)
