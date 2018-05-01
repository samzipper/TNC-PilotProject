## Navarro-SteadyState-WithPumping.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California. Starting head will be the 
# output from Navarro-SteadyState
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
modelname = 'Navarro-SteadyState-WithPumping'
modelname_NoPump = 'Navarro-SteadyState'
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

## update MNW2 package
# which well do you want to pump?
WellNum = 493
wellid = 'Well'+str(WellNum)

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

# update pumping rate for this well stress period data; 
mf.mnw2.mnw[wellid].stress_period_data['qdes'][0] = Qw
mf.mnw2.stress_period_data[0]['qdes'][mf.mnw2.stress_period_data[0]['wellid']==wellid] = Qw

# for some reason, MNW2 package reverses node data when a model is loaded and rewritten
# need to flip node_data so that topmost node is listed first
mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

## write inputs and run model
# write input datasets
mf.write_input()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## look at budget outputs
mfl = flopy.utils.MfListBudget(os.path.join(model_ws, modelname+".list"))
df_flux, df_vol = mfl.get_dataframes()

mfl_NoPump = flopy.utils.MfListBudget(os.path.join(model_ws_NoPump, modelname_NoPump+'.list'))
df_flux_NoPump, df_vol_NoPump = mfl_NoPump.get_dataframes()

# change in net fluxes
(df_flux['MNW2_IN'] - df_flux['MNW2_OUT']) - (df_flux_NoPump['MNW2_IN'] - df_flux_NoPump['MNW2_OUT'])
(df_flux['RIVER_LEAKAGE_IN'] - df_flux['RIVER_LEAKAGE_OUT']) - (df_flux_NoPump['RIVER_LEAKAGE_IN'] - df_flux_NoPump['RIVER_LEAKAGE_OUT'])

#### plot results ####
## look at output

## note: for some reason, FloPy doesn't change the prefix of output files, so they
# will all be in the model_ws folder bot have the modelname_NoPump prefix

# figure out timestep
time = 1

## head output
# Create the headfile object
h = bf.HeadFile(os.path.join(model_ws, modelname_NoPump+'.hds'), text='head')
h_NoPump = bf.HeadFile(os.path.join(model_ws_NoPump, modelname_NoPump+'.hds'), text='head')

# extract data matrix
head = h.get_data(totim=time)
head[head <= mf.bas6.hnoflo] = np.nan

head_NoPump = h_NoPump.get_data(totim=time)
head_NoPump[head_NoPump <= mf.bas6.hnoflo] = np.nan

# calculate WTE and DDN
wte = pp.get_water_table(head, nodata=mf.bas6.hnoflo)
wte_NoPump = pp.get_water_table(head_NoPump, nodata=mf.bas6.hnoflo)

ddn = wte_NoPump - wte

## load MNW output
mnwout = bf.CellBudgetFile(os.path.join(model_ws, modelname_NoPump+'.mnw2.out'), verbose=False)
mnwout_data = mnwout.get_data(totim=time, text='MNW2', full3D=True)[0]
mnwout.close()

# well of interest
row_wel = 410  # Well 493
col_wel = 360  # Well 493
mnwout_data[:,row_wel,col_wel]

# well
head[:,row_wel,col_wel]
head_NoPump[:,row_wel,col_wel]
wte[row_wel, col_wel]
wte_NoPump[row_wel, col_wel]

## make plots
plt.imshow(ddn)
plt.colorbar()

plt.imshow(ddn<0)

## process RIV data
# load RIV input
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')
iriv_ReachData = pd.read_table(os.path.join('modflow', 'input', 'iriv_ReachData.txt'), delimiter=' ')

# load RIV output
rivout = bf.CellBudgetFile(os.path.join(model_ws, modelname_NoPump+'.riv.out'))
rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
rivout.close()

rivout_NoPump = bf.CellBudgetFile(os.path.join(model_ws_NoPump, modelname_NoPump+'.riv.out'))
rivout_3D_NoPump = rivout_NoPump.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage_NoPump'] = rivout_3D_NoPump[0][iriv['lay'],iriv['row'],iriv['col']]
rivout_NoPump.close()

# join leakage to reach_data
iriv_ReachData = pd.merge(iriv_ReachData[['SegNum', 'row', 'col']], iriv[['row', 'col', 'leakage', 'leakage_NoPump']],
                          how='left', on=['row','col'])
                      
# summarize by segment number
iriv_out = iriv_ReachData.groupby('SegNum', as_index=False).agg({'leakage': 'sum', 'leakage_NoPump': 'sum'})

# calculate depletion and depletion.prc for each SegNum
MNW_net = (df_flux['MNW2_IN'] - df_flux['MNW2_OUT']) - (df_flux_NoPump['MNW2_IN'] - df_flux_NoPump['MNW2_OUT'])
iriv_out['depletion'] = iriv_out['leakage_NoPump'] - iriv_out['leakage']
iriv_out['depletion_prc'] = iriv_out['depletion']/MNW_net[0]

sum(iriv_out['depletion_prc'])
max(iriv_out['depletion_prc'])
min(iriv_out['depletion_prc'])