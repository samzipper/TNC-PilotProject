## MODFLOW_HTC_Navarro_SetUpIntermittentWithPumping.py
# This script creates a bunch of different MODFLOW input files with
# synthetic pumping wells.
#
# MAKE SURE YOU'VE RUN MODFLOW_HTC_Navarro_SetUpTransientWithPumping.py FIRST
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import shutil
import platform

# set up your model
modelname = 'Navarro-Intermittent'
modelname_transient = 'Navarro-Transient'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
stream_BC = 'SFR'     # 'RIV' or 'SFR'

# which wells to simulate?
every_n_wells = 7  # if you want all wells, just set this to 1

# which wells to simulate? there are 787 total
every_n_wells = 7  # if you want all wells, just set this to 1

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
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', stream_BC, modflow_v, model_prefix+'0')
model_ws_transient = os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v,  model_prefix+'0')
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname_transient+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=model_ws_transient)

# update model workspace
mf.name = modelname
mf.change_model_ws(model_ws)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')

### for the no pumping simulation (mf0), we can just write the input - no need to change anything
# for some reason, MNW2 package reverses node data when a model is loaded and rewritten
# need to flip node_data so that topmost node is listed first
mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

## write inputs for no-pumping scenario
mf.write_input()

# for some reason, FloPy doesn't update output files with new modelname
# do this manually by reading in .nam file
namfile = open(os.path.join(model_ws, modelname+'.nam'), 'r') 
namtext = namfile.read()
namfile.close()

namfile = open(os.path.join(model_ws, modelname+'.nam'), 'w') 
namtext = namtext.replace(os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v, model_prefix+'0', ''), '')  # fix ddn - not sure why this is happening
namtext = namtext.replace(modelname_transient, modelname)  # fix modelname
namfile.write(namtext)
namfile.close()

## copy launch and postprocessing scripts
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', 'launch_thisRun.sh'), 
        os.path.join(model_ws, 'launch_thisRun.sh'))

shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', 'postprocess_thisRun_'+stream_BC+'.py'), 
        os.path.join(model_ws, 'postprocess_thisRun.py'))

### now, go through the pumping well scenarios - 
### the only package that will change is MNW2 for pumping
### Now: scroll through wells and create WEL input files
# all wells are already set up in MODFLOW_Navarro-SteadyState script so 
# we just have to turn them on

# grab number of stress periods
nper = mf.dis.nper

# read in well data
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# which stress periods should pumping occur? (index; 0-based indexing)
well_pump_mo = [6, 7, 8, 9, 10]
well_pump_yr = [False]*12
for i in well_pump_mo:
    well_pump_yr[i-1] = True
well_pump_sp = well_pump_yr*int(nper/12)

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

for w in range(0,iwel.shape[0], every_n_wells):
    WellNum = iwel['WellNum'][w]
    wellid = 'Well'+str(WellNum)

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', stream_BC, modflow_v, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

    # update pumping rate for this well stress period data; 
    for sp in range(0, nper):
        if well_pump_sp[sp]:
            mf.mnw2.mnw[wellid].stress_period_data['qdes'][sp] = Qw
            mf.mnw2.stress_period_data[sp]['qdes'][mf.mnw2.stress_period_data[sp]['wellid']==wellid] = Qw
        
        # update itmp
        if (sp==0):
            mf.mnw2.itmp[sp] = mf.mnw2.itmp[0]
        elif (well_pump_sp[sp]==well_pump_sp[sp-1]):
            mf.mnw2.itmp[sp] = -1
        else:
            mf.mnw2.itmp[sp] = mf.mnw2.itmp[0]
    
    # write input (NAM and MNW2 only)
    mf.write_input(SelPackList=['MNW2'])

    # turn off this well
    for sp in range(0, nper):
        mf.mnw2.mnw[wellid].stress_period_data['qdes'][sp] = 0
        mf.mnw2.stress_period_data[sp]['qdes'][mf.mnw2.stress_period_data[sp]['wellid']==wellid] = 0
    
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', 'launch_thisRun.sh'),
                 os.path.join(w_model_ws, 'launch_thisRun.sh'))

    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', 'postprocess_thisRun_'+stream_BC+'.py'),
                 os.path.join(w_model_ws, 'postprocess_thisRun.py'))

    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'Intermittent', stream_BC, modflow_v, modelname+'_Template.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))
