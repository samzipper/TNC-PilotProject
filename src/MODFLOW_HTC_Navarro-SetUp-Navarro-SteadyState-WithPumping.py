## MODFLOW_HTC_Navarro-SetUp-SteadyState-WithPumping.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
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

# where is your MODFLOW-2005 executable?
if (modflow_v=='mf2005'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'
    else:
        path2mf = modflow_v
elif (modflow_v=='mfnwt'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MODFLOW-NWT_1.1.3/bin/MODFLOW-NWT.exe'
    else:
        path2mf = modflow_v

# do you want to use RIV or SFR for stream features?
stream_BC = 'RIV'  # options: 'RIV' 'SFR'

# check if model workspace exists; create if not
model_prefix = 'mf'
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, model_prefix+'0')
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=os.path.join('modflow', modelname))

# update model workspace
mf.change_model_ws(model_ws)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')
   
## create WEL package, but don't pump anything
wel = flopy.modflow.mfwel.ModflowWel(mf, stress_period_data={0: [0,50,50,0]},
                                     ipakcb=71, filenames=[modelname+'.wel', modelname+'.wel.out'])

# set up solver depending on version of MODFLOW
tol_head = 1e-2
if (modflow_v=='MF2005'):
    pcg = flopy.modflow.ModflowPcg(mf, hclose=tol_head, rclose=tol_head)
elif (modflow_v=='MFNWT'):
    # linmeth has two matrix solver options (1 or 2)
    nwt = flopy.modflow.ModflowNwt(mf, headtol=tol_head, linmeth=2, options='COMPLEX')
                
## write inputs for no-pumping scenario
mf.write_input()

## copy launch script
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'), 
        os.path.join(model_ws, 'launch_thisRun.sh'))

### Now: scroll through wells and create WEL input files
## load wel boundary condition info
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

for w in range(0,iwel.shape[0]):
    WellNum = iwel['WellNum'][w]

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

    # set up stress period data
    wel.stress_period_data = {0: [iwel['lay'][w], iwel['row'][w], iwel['col'][w], Qw]}
    wel.filenames=[modelname+'.wel', modelname+'.wel.out']

    # write input (NAM and WEL only)
    mf.write_input(SelPackList=['WEL'])
        
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'),
            os.path.join(w_model_ws, 'launch_thisRun.sh'))


    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modelname+'_Template.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))

