## MODFLOW_HTC_Navarro-SetUp-SteadyState-WithPumping.py
# This script creates a bunch of different MODFLOW input files with
# synthetic pumping wells.
#
# MAKE SURE YOU'VE RUN MODFLOW_Navarro-SteadyState.py FIRST
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import shutil
import platform

# set up your model
modelname = 'Navarro-SteadyState'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
stream_BC = 'RIV'  # options: 'RIV' 'SFR'

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
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v, model_prefix+'0')
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=os.path.join('modflow', modelname, stream_BC, modflow_v))

# update model workspace
mf.change_model_ws(model_ws)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')
                
## write inputs for no-pumping scenario
mf.write_input()

## copy launch script
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'), 
        os.path.join(model_ws, 'launch_thisRun.sh'))

### Now: scroll through wells and create WEL input files
# all wells are already set up in MODFLOW_Navarro-SteadyState script so 
# we just have to turn them on

# read in well data
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

for w in range(0,3):#iwel.shape[0]):
    WellNum = iwel['WellNum'][w]
    wellid = 'Well'+str(WellNum)

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

    # update pumping rate for this well stress period data; 
    mf.mnw2.mnw[wellid].stress_period_data['qdes'][0] = Qw
    mf.mnw2.stress_period_data[0]['qdes'][mf.mnw2.stress_period_data[0]['wellid']==wellid] = Qw
    
    # for some reason, MNW2 package reverses node data when a model is loaded and rewritten
    # need to flip node_data so that topmost node is listed first
    mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

    # write input (NAM and WEL only)
    mf.write_input(SelPackList=['MNW2'])

    # turn off this well
    mf.mnw2.mnw[wellid].stress_period_data['qdes'][0] = 0
    mf.mnw2.stress_period_data[0]['qdes'][mf.mnw2.stress_period_data[0]['wellid']==wellid] = 0
        
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'),
            os.path.join(w_model_ws, 'launch_thisRun.sh'))

    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modelname+'_Template_'+modflow_v+'.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))