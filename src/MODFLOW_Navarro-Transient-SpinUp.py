## Navarro-Transient-SpinUp.py
# This script loads the Navarro-SteadyState groundwater model and converts it
# to a transient model (monthly stress periods) and runs a spin-up simulation.
#
# You need to have already run the script MODFLOW_Navarro-SteadyState.py for 
# the modflow_v and stream_BC you are using.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import flopy.utils.binaryfile as bf
import platform
import copy

# set up your model
modelname = 'Navarro-Transient-SpinUp'
modelname_SS = 'Navarro-SteadyState'
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
model_ws_SS = os.path.join('modflow', modelname_SS, stream_BC, modflow_v)
model_ws = os.path.join('modflow', modelname, stream_BC, modflow_v)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Load steady-state model
mf = flopy.modflow.Modflow.load(modelname_SS+'.nam', 
        exe_name=path2mf, version=modflow_v, 
        model_ws=model_ws_SS)

# rename and update model workspace
mf.name = modelname
mf.change_model_ws(model_ws)

### the following packages have stress period data: DIS, UPW, RCH, OC, RIV/SFR, MNW2
## update DIS

# parameters controlling time discretization
numyears = 20                # number of years for transient simulation (1st year will always be SS)
sp_per_year = 12             # dry season and wet season for now
sp_length_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # 5 months (wet), 7 months (dry)
sp_season = ['wet']*4 + ['dry']*7 + ['wet']  # must be same length as sp_length_days
ts_length_days = 5          # number of days per timestep - will be approximate because nstp must be integer

# define stress period data
nper = 1+numyears*sp_per_year
perlen = [365]+sp_length_days*numyears
nstp_yr = [round(elem, 0) for elem in (np.array(sp_length_days, dtype='f')/ts_length_days).tolist()]
nstp = [1]+nstp_yr*numyears
steady = [True]+[False]*(nper-1)

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

## update UPW
sy = 0.10   # specific yield (using 50% of domain mean porosity for now)
ss = 1e-5   # specific storage
hani = 1.0  # ratio of hk along columns to hk along rows
mf.upw.sy = sy
mf.upw.ss = ss
mf.upw.hani = hani

## update RCH
rchrate_yr = 150/(1000*365)  # [mm/yr] --> [m/d]
rch_days = sum([sp_length_days[i] for i,x in enumerate(sp_season) if x=='wet'])    # number of days over which recharge occurs
rchrate_wet = rchrate_yr*365/rch_days  # condense annual recharge rate into seasonal recharge [m/d]
rchrate_dry = 0.

# set recharge rate for each stress period
sp_season_all = ['SS']+sp_season*numyears
rech = {}
rech[0] = rchrate_yr
for sp in np.arange(1, len(sp_season_all)):
    if sp_season_all[sp]=='wet': rech[sp] = rchrate_wet
    if sp_season_all[sp]=='dry': rech[sp] = rchrate_dry
mf.rch.rech = rech

## update OC - save at end of every stress period
#oc = flopy.modflow.ModflowOc(mf, save_every=True, compact=True)
oc_spd = {}
oc_spd[(0,0)] = ['save head', 'save budget']
for sp in range(1,nper):
    oc_spd[(sp,(nstp[sp]-1))] = ['save head', 'save budget']
oc = flopy.modflow.ModflowOc(mf, stress_period_data=oc_spd, compact=True)

## update MNW2
# for some reason, MNW2 package reverses node data when a model is loaded and rewritten
# need to flip node_data so that topmost node is listed first
mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

# update stress period data
mf.mnw2.nper = nper
#mnw_spd = mf.mnw2.stress_period_data[0]
itmp = mf.mnw2.itmp*nper
mf.mnw2.itmp = itmp

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

## write inputs and run model
mf.write_input()

# for some reason, FloPy doesn't update output files with new modelname
# do this manually by reading in .nam file
namfile = open(os.path.join(model_ws, modelname+'.nam'), 'r') 
namtext = namfile.read()
namfile.close()

namfile = open(os.path.join(model_ws, modelname+'.nam'), 'w') 
namtext = namtext.replace(modelname_SS, modelname)
namfile.write(namtext)
namfile.close()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## save budget output to check equilibrium
# load budget
mfl = flopy.utils.MfListBudget(os.path.join(model_ws, modelname+".list"))
df_flux, df_vol = mfl.get_dataframes()

# extract net leakage
if (stream_BC == 'RIV'):
    stream_prefix = 'RIVER'
if (stream_BC == 'SFR'):
    stream_prefix = 'STREAM'

# summary columns
df_flux['datetime'] = df_flux.index[:]
df_flux['leakage'] = df_flux[stream_prefix+'_LEAKAGE_IN'] - df_flux[stream_prefix+'_LEAKAGE_OUT']

# save budget output
df_flux.to_csv(os.path.join(model_ws, modelname+'_SummarizeBudget.csv'), index=False)

## load heads
h = bf.HeadFile(os.path.join(model_ws, modelname+'.hds'), text='head')
times = h.get_times()
heads = h.get_data(totim=times[-1])
h.close()

# save head output
for l in range(0,nlay):
    np.savetxt(os.path.join(model_ws, 'head_layer'+str(l)+'.csv'), heads[l,:,:], fmt='%6.2f', delimiter=',')