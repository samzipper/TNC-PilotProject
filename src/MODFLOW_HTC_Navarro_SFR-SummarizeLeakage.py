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

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v)
prefix_runs = 'mf'
modelname = 'Navarro-SteadyState'
time=1

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
        ## make a copy of iriv to avoid altering original
        isfr_SegmentData_copy = isfr_SegmentData.copy()
        
        ## SFR output
        sfrout = sf.SfrFile(os.path.join(dir_runs, prefix_runs+str(w), modelname+'.sfr.out'))
        sfr_df = sfrout.get_dataframe()
        
        # summarize by segment
        sfr_df_summary = sfr_df.groupby('segment', as_index=False).agg({'Qaquifer': 'sum'})
        
        # join with isfr_SegmentData
        sfr_merge = pd.merge(isfr_SegmentData_copy, sfr_df_summary,
                             left_on=['SFR_NSEG'], right_on=['segment'])
                             
        # add well number column
        sfr_merge['WellNum'] = w

        ## save RIV output
        np.savetxt(os.path.join(dir_runs, prefix_runs+str(w), 'SFR_leakage.txt'), 
                   sfr_merge[['SFR_NSEG', 'Qaquifer', 'WellNum']],
                   header="SFR_NSEG,leakage,WellNum", comments='')
        
        ## add to overall data frame
        if (start_flag):
            sfr_all = sfr_merge[['SFR_NSEG', 'Qaquifer', 'WellNum']]
            start_flag = False
        else:
            sfr_all = sfr_all.append(sfr_merge[['SFR_NSEG', 'Qaquifer', 'WellNum']])

    ## status update
    print(w, ' complete')

## save all output
np.savetxt(os.path.join(dir_runs, 'SFR-SummarizeLeakage.csv'), sfr_all,
           delimiter=",", header="SFR_NSEG,leakage,WellNum", comments='')

## print total leakage across all segments
sfr_all_summary = sfr_all.groupby('WellNum', as_index=False).agg({'Qaquifer': 'sum'})
print(sfr_all_summary.head())
