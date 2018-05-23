## MODFLOW_HTC_Navarro_SetUpTransientWithPumping.py
# This script creates a bunch of different MODFLOW input files with
# synthetic pumping wells.
#
# MAKE SURE YOU'VE RUN MODFLOW_Navarro-Transient-SpinUp.py FIRST
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import shutil
import platform

# set up your model
modelname = 'Navarro-Transient'
modelname_NoPump = 'Navarro-Transient-SpinUp'
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
model_ws_NoPump = os.path.join('modflow', modelname_NoPump, stream_BC, modflow_v)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname_NoPump+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=model_ws_NoPump)

# update model workspace
mf.name = modelname
mf.change_model_ws(model_ws)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')

### the following packages have stress period data: DIS, RCH, OC, RIV/SFR, MNW2
## assumption is that the length of the spin-up simulation will always be longer
## than the pumping simulation, so we can just take a subset of the stress period 
## data from the transient spin-up.
##
## remember that the first SP of the transient spin-up is SS

## update DIS

# parameters controlling time discretization
numyears = 2                # number of years for transient simulation
sp_per_year = 12             # dry season and wet season for now
ts_length_days = 1          # number of days per timestep - will be approximate because nstp must be integer

# stress periods
nper = numyears*sp_per_year
perlen = mf.dis.perlen[1:(nper+1)]
steady = [False]*(nper)

# timesteps
nstp = [round(elem, 0) for elem in (np.array(perlen, dtype='f')/ts_length_days).tolist()]

# grab other dis info from old model
nlay = mf.dis.nlay
ncol = mf.dis.ncol
nrow = mf.dis.nrow
delr = mf.dis.delr
delc = mf.dis.delc
ztop = mf.dis.top
botm = mf.dis.botm

# make new dis object
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, 
                               delr=delr, delc=delc,
                               top=ztop, botm=botm,
                               nper=nper, perlen=perlen, 
                               nstp=nstp, steady=steady)

## update RCH
rech = {}
for sp in range(0, nper):
    rech[sp] = mf.rch.rech[sp+1][0,0]   # sp+1 because first sp of SpinUp is steady-state
mf.rch.rech = rech

## update OC
oc_spd = {}
ts_counter = 0
for sp in range(0,nper):
    for stp in range(0,int(nstp[sp])):
        ts_counter = ts_counter+1
        if (ts_counter % 5 == 0): oc_spd[sp,stp] = ['save budget']
oc = flopy.modflow.ModflowOc(mf, stress_period_data=oc_spd, compact=True)

## update starting conditions from spinup
strt = np.empty([nlay, nrow, ncol])
for l in range(0,nlay):
    strt[l,:,:] = np.array(pd.read_csv(os.path.join(model_ws_NoPump, 'head_layer'+str(l)+'.csv'),
        header=None, delim_whitespace=False))
mf.bas6.strt = strt

## update streams
if (stream_BC=='RIV'):
    ## don't need to do anything unless RIV properties (e.g. stage, conductance)
    ## change between stress periods
    print(stream_BC)

if (stream_BC=='SFR'):
    ## number of active streams per stress period
    sfr_d5 = {}
    for i in range(0,nper): sfr_d5[i] = mf.sfr.dataset_5[0]
    mf.sfr.dataset_5 = sfr_d5
    
    ## update segment data with quickflow
    sfr_segData = mf.sfr.segment_data.copy()
    sfr_segData_dict = {(k-1): sfr_segData[k] for k in sorted(sfr_segData.keys())[1:(nper+1)]}
    mf.sfr.segment_data = sfr_segData_dict

## MNW2: first we are doing the no-pumping simulation so just change stress period number
# grab data from transient spin-up run
itmp = mf.mnw2.itmp
mnw_node_data = mf.mnw2.node_data

mnw_spd = {}
for sp in range(0,nper):
    mnw_spd[sp] =  mf.mnw2.stress_period_data[0].copy()

mnw2 = flopy.modflow.ModflowMnw2(model=mf, mnwmax=itmp[0], 
                                 node_data=mnw_node_data,
                                 stress_period_data=mnw_spd,
                                 itmp=itmp[0:1]*nper, ipakcb=71,
                                 filenames=[modelname+'.mnw2', modelname+'.mnw2.out'])

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
namtext = namtext.replace(os.path.join('modflow', modelname_NoPump, stream_BC, modflow_v, modelname_NoPump+'.ddn'), modelname+'.ddn')  # fix ddn - not sure why this is happening
namtext = namtext.replace(modelname_NoPump, modelname)  # fix modelname
namfile.write(namtext)
namfile.close()

# run model for testing
#mf.run_model()

## copy launch script
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Transient', 'launch_thisRun.sh'), 
        os.path.join(model_ws, 'launch_thisRun.sh'))

### Now: scroll through wells and create WEL input files
# all wells are already set up in MODFLOW_Navarro-SteadyState script so 
# we just have to turn them on

# read in well data
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# when should it start pumping? (which stress period; 0-based indexing)
well_start_sp = 4

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

every_n_wells = 25  # if you want all wells, just set this to 1
for w in range(0,iwel.shape[0], every_n_wells):
#for w in range(0,3):
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
    
    # write input (NAM and MNW2 only)
    mf.write_input(SelPackList=['MNW2'])

    # turn off this well
    for sp in range(well_start_sp, nper):
        mf.mnw2.mnw[wellid].stress_period_data['qdes'][sp] = 0
        mf.mnw2.stress_period_data[sp]['qdes'][mf.mnw2.stress_period_data[sp]['wellid']==wellid] = 0
        
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'Transient', 'launch_thisRun.sh'),
            os.path.join(w_model_ws, 'launch_thisRun.sh'))

    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'Transient', stream_BC, modflow_v, modelname+'_Template.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))
