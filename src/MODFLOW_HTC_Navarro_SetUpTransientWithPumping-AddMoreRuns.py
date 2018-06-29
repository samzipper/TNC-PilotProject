## MODFLOW_HTC_Navarro_SetUpTransientWithPumping-AddMoreRuns.py
# This script creates a bunch of different MODFLOW input files with
# synthetic pumping wells.
#
# It is intended to add additional runs to those already created using 
# MODFLOW_HTC_Navarro_SetUpTransientWithPumping.py
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters):q

import os
import numpy as np
import flopy
import pandas as pd
import shutil
import platform

## which wells do you want to add?
# note that these are indices within the iwel data frame; the actual WellNum is w+1
every_5 = np.arange(0,787,5)
every_25 = np.arange(0,787,25)
w_to_add = [x for x in every_5 if x not in every_25]

# set up your model
modelname = 'Navarro-Transient'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
stream_BC = 'SFR'     # 'RIV' or 'SFR'

# where is your MODFLOW-2005 executable?
if (modflow_v=='mf2005'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'
    else: 
        path2mf = modflow_v
elif (modflow_v=='mfnwt'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MODFLOW-NWT_1.1.4/bin/MODFLOW-NWT.exe'
    else:
        path2mf = modflow_v

# check if model workspace exists; create if not
model_prefix = 'mf'
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v, model_prefix+'0')

print(modelname)
print(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=model_ws)

# grab some info that's needed later on
nper = mf.dis.nper

# set up itmp
mf.mnw2.itmp[1:(nper-1)] = [-1]*(nper-2)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')

### Now: scroll through wells and create WEL input files
# all wells are already set up in MODFLOW_Navarro-SteadyState script so 
# we just have to turn them on

# read in well data
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# when should it start pumping? (which stress period; 0-based indexing)
well_start_sp = 4

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

for w in w_to_add:
    WellNum = iwel['WellNum'][w]
    wellid = 'Well'+str(WellNum)

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

    # update pumping rate for this well stress period data; 
    for sp in range(well_start_sp, nper):
        mf.mnw2.mnw[wellid].stress_period_data['qdes'][sp] = Qw
        mf.mnw2.stress_period_data[sp]['qdes'][mf.mnw2.stress_period_data[sp]['wellid']==wellid] = Qw
    
    # update itmp
    mf.mnw2.itmp[well_start_sp] = mf.mnw2.itmp[0]
    
    # write input (NAM and MNW2 only)
    mf.write_input(SelPackList=['MNW2'])

    # turn off this well
    for sp in range(well_start_sp, nper):
        mf.mnw2.mnw[wellid].stress_period_data['qdes'][sp] = 0
        mf.mnw2.stress_period_data[sp]['qdes'][mf.mnw2.stress_period_data[sp]['wellid']==wellid] = 0
    mf.mnw2.itmp[well_start_sp] = -1
    
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Transient', 'launch_thisRun.sh'),
                 os.path.join(w_model_ws, 'launch_thisRun.sh'))

    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Transient', 'postprocess_thisRun_'+stream_BC+'.py'),
                 os.path.join(w_model_ws, 'postprocess_thisRun.py'))

    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v, modelname+'_Template.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))
