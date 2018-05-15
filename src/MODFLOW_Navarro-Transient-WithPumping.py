## Navarro-Transient-WithPumping.py
# This script loads a transient groundwater flow model for the 
# Navarro River Watershed in California and pumps a well. 
# Starting head will be the output from Navarro-Transient.
#
# You need to have already run the script MODFLOW_Navarro-SteadyState.py for 
# the modflow_v and stream_BC you are using.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import platform

# set up your model
modelname = 'Navarro-Transient-WithPumping'
modelname_NoPump = 'Navarro-Transient-SpinUp'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
stream_BC = 'RIV'     # 'RIV' or 'SFR'

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
model_ws_NoPump = os.path.join('modflow', modelname_NoPump, stream_BC, modflow_v)
model_ws = os.path.join('modflow', modelname, stream_BC, modflow_v)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname_NoPump+'.nam', 
        exe_name=path2mf, version=modflow_v, 
        model_ws=model_ws_NoPump)

# rename and update model workspace
mf.name = modelname
mf.change_model_ws(model_ws)

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
for sp in range(1, (nper+1)):
    rech[sp] = mf.rch.rech[sp][0,0]
mf.rch.rech = rech

## update OC
oc_spd = {}
for sp in range(0,nper):
    for stp in range(int(nstp[sp]-1),int(nstp[sp])):
        oc_spd[sp,stp] = ['save head', 'save budget']
oc = flopy.modflow.ModflowOc(mf, stress_period_data=oc_spd, filenames=modelname, compact=True)

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
    
    ## set up quickflow 
    # monthly quickflow [mm/mo] from the script Navarro_StreamflowData.R
    quickflow_mo = [119.61, 91.86, 70.45, 27.65, 5.87, 1.64, 0.52, 0.21, 0.21, 2.04, 16.92, 78.91]
    quickflow_mo_prc = [quickflow_mo[i]/sum(quickflow_mo) for i in range(0,len(quickflow_mo))]
    quickflow_mo_prc_sp = [1]+quickflow_mo_prc*numyears
    
    ## update SFR
    SS_segData = mf.sfr.segment_data[0].copy()
    sfr_segData = {}
    for i in range(0,nper):
        sfr_segData[i] = SS_segData.copy()
        sfr_segData[i]['runoff'] = sfr_segData[i]['runoff']*quickflow_mo_prc_sp[i]
    mf.sfr.segment_data = sfr_segData

## update MNW2 package
# which well do you want to pump?
WellNum = 493
wellid = 'Well'+str(WellNum)

# when should it start pumping? (which stress period; 0-based indexing)
well_start_sp = 4

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

# grab data from transient spin-up run
itmp = mf.mnw2.itmp
mnw_node_data = np.flipud(mf.mnw2.node_data)

mnw_spd = {}
for sp in range(0,nper):
    mnw_spd[sp] =  mf.mnw2.stress_period_data[0].copy()

# update pumping rate for this well stress period data; 
for sp in range(well_start_sp, nper):
    mnw_spd[sp]['qdes'][mnw_spd[sp]['wellid']==wellid] = Qw

mnw2 = flopy.modflow.ModflowMnw2(model=mf, mnwmax=itmp[0], 
                                 node_data=mnw_node_data,
                                 stress_period_data=mnw_spd,
                                 itmp=itmp[0:1]*nper, ipakcb=71,
                                 filenames=[modelname+'.mnw2', modelname+'.mnw2.out'])

# for some reason, MNW2 package reverses node data when a model is loaded and rewritten
# need to flip node_data so that topmost node is listed first
mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

## write inputs and run model
# write input datasets
mf.write_input()

# for some reason, FloPy doesn't update output files with new modelname
# do this manually by reading in .nam file
namfile = open(os.path.join(model_ws, modelname+'.nam'), 'r') 
namtext = namfile.read()
namfile.close()

namfile = open(os.path.join(model_ws, modelname+'.nam'), 'w') 
namtext = namtext.replace('modflow\\'+modelname_NoPump+'\\'+stream_BC+'\\'+modflow_v+'\\', '')  # fix ddn - not sure why this is happening
namtext = namtext.replace(modelname_NoPump, modelname)  # fix modelname
namfile.write(namtext)
namfile.close()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')