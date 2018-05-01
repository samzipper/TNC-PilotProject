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
model_ws = os.path.join('modflow', modelname, stream_BC, modflow_v)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load('Navarro-SteadyState.nam', 
        exe_name=path2mf, version=modflow_v, 
        model_ws=os.path.join('modflow', 'Navarro-SteadyState', stream_BC, modflow_v))

## set up initial conditions
h = bf.HeadFile(os.path.join(mf.model_ws, 'Navarro-SteadyState.hds'), text='head')
head = h.get_data(totim=1)
mf.bas6.start = head

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
# reverse node_data so that topmost node is listed first
for w in range(1,(mf.mnw2.mnwmax+1)):
    mf.mnw2.mnw['Well'+str(w)].node_data = np.flipud(mf.mnw2.mnw['Well'+str(w)].node_data)
mf.mnw2.node_data = np.flipud(mf.mnw2.node_data)

## write inputs and run model
# write input datasets
mf.write_input()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## look at budget outputs
mfl_pumped = flopy.utils.MfListBudget(os.path.join(model_ws, modelname+".list"))
df_flux_pumped, df_vol_pumped = mfl_pumped.get_dataframes()
print(df_flux_pumped)

mfl = flopy.utils.MfListBudget(os.path.join('modflow', 'Navarro-SteadyState', stream_BC, modflow_v, 'Navarro-SteadyState.list'))
df_flux, df_vol = mfl.get_dataframes()
print(df_flux)

df_flux_pumped - df_flux

#### plot results ####
## look at output
# figure out timestep
time = 1

## head output
# Create the headfile object
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')

# extract data matrix
head_pumped = h.get_data(totim=time)
head_pumped[head_pumped <= mf.bas6.hnoflo] = np.nan

# calculate WTE and DDN
wte_pumped = pp.get_water_table(head_pumped, nodata=mf.bas6.hnoflo)
ddn = wte - wte_pumped

## load MNW output
mnwout = bf.CellBudgetFile(os.path.join(model_ws, modelname+'.mnw2.out'), verbose=False)
mnwout_data = mnwout.get_data(totim=time, text='MNW2', full3D=True)[0]
mnwout_data[:,row_wel,col_wel]

# well
head[:,row_wel,col_wel]
head_pumped[:,row_wel,col_wel]
wte[row_wel, col_wel]
wte_pumped[row_wel, col_wel]

# closest point on RIV SegNum 220; RIV elevation is 243.28 m
head[:,410,359]
head_pumped[:,410,359]
wte[410,359]
wte_pumped[410,359]

mnwout.close()

## make plots
plt.imshow(ddn)
plt.colorbar()

plt.imshow(ddn<0)
