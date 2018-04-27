## MODFLOW_HTC_Navarro_RIV-SummarizeLeakage.py
# This script will open RIV output files and summarize the leakage by SegNum.
# Needs the script MODFLOW_HTC_Navarro_CheckFailures.py to be run first.

import os
import numpy as np
import flopy
import pandas as pd

# what stream BC and modflow version?
stream_BC = 'SFR'  # 'RIV' or 'SFR'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v)
prefix_runs = 'mf'
modelname = 'Navarro-SteadyState'
time=1

## figure out which runs succeeded (output from CheckFailures script)
succ = pd.read_table(os.path.join(dir_runs, 'CheckFailure.csv'), delimiter=",")

## loop through and make output data frame
start_flag = True
for w in succ.WellNum:
    # check if model ran successfully
    converge = succ.Success.loc[(succ['WellNum']==w)].bool()
    if converge:
        ## load budget
        mfl = flopy.utils.MfListBudget(os.path.join(dir_runs, prefix_runs+str(w), modelname+".list"))
        df_flux, df_vol = mfl.get_dataframes()

        # add column for well name
        df_flux['WellNum'] = 0
       
        ## add to overall data frame
        if (start_flag):
            df_flux_all = df_flux
            start_flag = False
        else:
            df_flux_all = df_flux_all.append(df_flux)

    ## status update
    print(w, ' complete')

## save and print head to screen
df_flux_all.to_csv(os.path.join(dir_runs, "SummarizeBudgetTerms.csv"), index=False)
print(df_flux_all.head())
