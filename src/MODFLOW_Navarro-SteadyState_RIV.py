## Navarro-SteadyState_RIV.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd

# set up your model
modelname = 'Navarro-SteadyState'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'

# where is your MODFLOW-2005 executable?
if (modflow_v=='mf2005'):
    path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'
elif (modflow_v=='mfnwt'):
    path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MODFLOW-NWT_1.1.3/bin/MODFLOW-NWT.exe'

# check if model workspace exists; create if not
model_ws = os.path.join('modflow', modelname)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

if not os.path.isdir(os.path.join(model_ws, 'output')):
    os.makedirs(os.path.join(model_ws, 'output'))

# Assign name and create modflow model object
mf = flopy.modflow.Modflow(modelname, exe_name=path2mf, 
                           model_ws=model_ws, version=modflow_v)

## Set up DIS and BAS
# read in text output from R for ibound
ibound = np.array(pd.read_csv(os.path.join('modflow', 'input', 'ibound.txt'),
                              header=None, delim_whitespace=True), dtype=np.int32)
                              
# discretization (space) - these should be the same as in your R script                         
nlay = 1
nrow = ibound.shape[0]
ncol = ibound.shape[1]
delr = 100
delc = 100

# discretization (time)
nper = 1
perlen = [1]
nstp = [1]
steady = [True]

# read in text output from R for ztop
ztop = np.array(pd.read_csv(os.path.join('modflow', 'input', 'ztop.txt'),
                            header=None, delim_whitespace=True))

# define starting head (set to land surface elevation for now)
strt = ztop

# define bottom elevation (-100 m everywhere)
zbot = -100

# make MODFLOW objects
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, 
                               delr=delr, delc=delc,
                               top=ztop, botm=zbot,
                               nper=nper, perlen=perlen, 
                               nstp=nstp, steady=steady)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=0)

## flow properties
# properties
hk = 1e-12*1e7*86400     # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
vka = 1.    # anisotropy
sy = 0.10    # specific yield (using 50% of domain mean porosity for now)
ss = 1e-5  # specific storage
laytyp = 1  # layer type

# make MODFLOW objects
lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp)

# set up solver depending on version of MODFLOW
tol_head = 1e-2
if (modflow_v=='MF2005'):
    pcg = flopy.modflow.ModflowPcg(mf, hclose=tol_head, rclose=tol_head)
elif (modflow_v=='MFNWT'):
    # linmeth has two matrix solver options (1 or 2)
    nwt = flopy.modflow.ModflowNwt(mf, headtol=tol_head, linmeth=2, options='COMPLEX')

## recharge
# long-term average baseflow is 150 mm/yr
rchrate = 150/(1000*365)  # [mm/yr] --> [m/d]
rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=3)

## output control
spd = {(0, 0): ['save head', 'save budget', 'save drawdown', 'print head', 'print budget', 'print drawdown']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

## river boundary condition
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')

# constant domain parameters
depth = 4  # river depth?
riverbed_K = hk/10
river_width = 10
riverbed_thickness = 1

# estimate conductance based on: river width, river length, riverbed thickness, riverbed K
iriv['cond'] = round(riverbed_K*river_width*iriv['totalLength_m']*riverbed_thickness)   # river bottom conductance? 

# empty list to hold stress period data
riv_list = []

# populate list
for r in range(0,iriv.shape[0]):
    riv_list.append([iriv['lay'][r], iriv['row'][r], iriv['col'][r], 
                     iriv['elev_m_min'][r], iriv['cond'][r], iriv['elev_m_min'][r]-depth])    
riv_spd = {0: riv_list}

# make MODFLOW object
riv = flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd, ipakcb=61, 
                               filenames=[modelname+'.riv', modelname+'.riv.out'])

## write inputs and run model
# write input datasets
mf.write_input()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')