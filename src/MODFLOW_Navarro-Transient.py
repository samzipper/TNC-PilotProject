## Navarro-Transient.py
# This script loads the Navarro-SteadyState groundwater model and converts it
# to a transient model (monthly stress periods).
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
import flopy.utils.sfroutputfile as sf
import flopy.utils.postprocessing as pp
import platform
import matplotlib.pyplot as plt

# set up your model
modelname = 'Navarro-Transient'
modelname_SS = 'Navarro-SteadyState'
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

### the following packages have stress period data: DIS, RCH, OC, RIV/SF, MNW2
## update DIS

# parameters controlling time discretization
numyears = 2                # number of years for transient simulation (1st year will always be SS)
sp_per_year = 2             # dry season and wet season for now
sp_length_days = [150, 215] # 5 months (wet), 7 months (dry)
sp_season = ['wet', 'dry']  # must be same length as sp_length_days
ts_length_days = 5          # number of days per timestep - should divide evenly
                            #   into all numbers in sp_length_days

# define stress period data
nper = 1+numyears*sp_per_year
perlen = [365]+sp_length_days*numyears
nstp_yr = (np.array(sp_length_days, dtype='f')/ts_length_days).tolist()
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

## update OC - just save at every timestep for now
oc = flopy.modflow.ModflowOc(mf, save_every=True, compact=True)
#oc_spd = {(0, 0): ['save head'],
#          (1, 2): ['save head'],
#          (2, 2): ['save head']}
#mf.oc.stress_period_data = oc_spd

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

if (stream_BC=='SFR'):
    ## update SFR

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

## look at budget outputs
mfl = flopy.utils.MfListBudget(os.path.join(model_ws, modelname+".list"))
df_flux, df_vol = mfl.get_dataframes()

plt.plot(range(1,62), df_flux['RECHARGE_IN'])
plt.plot(range(1,62), df_flux['RIVER_LEAKAGE_IN'])
plt.plot(range(1,62), df_flux['RIVER_LEAKAGE_OUT'])
plt.plot(range(1,62), df_flux['RIVER_LEAKAGE_IN']-df_flux['RIVER_LEAKAGE_OUT'])


## note: for some reason, FloPy doesn't change the prefix of output files, so they
# will all be in the model_ws folder bot have the modelname_NoPump prefix

# figure out timestep
time = 1

## head output
# Create the headfile object
h = bf.HeadFile(os.path.join(model_ws, modelname+'.hds'), text='head')
times = h.list_records()

# extract head and calculate wte
head_SS = h.get_data(kstpkper=(0,0))
head_SS[head_SS <= mf.bas6.hnoflo] = np.nan
wte_SS = pp.get_water_table(head_SS, nodata=mf.bas6.hnoflo)

head_sp1 = h.get_data(kstpkper=(29,1))
head_sp1[head_sp1 <= mf.bas6.hnoflo] = np.nan
wte_sp1 = pp.get_water_table(head_sp1, nodata=mf.bas6.hnoflo)

head_sp2 = h.get_data(kstpkper=(29,2))
head_sp2[head_sp2 <= mf.bas6.hnoflo] = np.nan
wte_sp2 = pp.get_water_table(head_sp2, nodata=mf.bas6.hnoflo)

plt.imshow(wte_SS); plt.colorbar()
plt.imshow(wte_SS-wte_sp1); plt.colorbar()
plt.imshow(wte_SS-wte_sp2); plt.colorbar()
plt.imshow(wte_sp1-wte_sp2); plt.colorbar()