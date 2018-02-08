## Navarro-SteadyState-WithPumping.py
# This script loads the steady-state groundwater flow model for the 
# Navarro River Watershed in California created by the script
# Navarro-SteadyState.py and pumps a series of synthetic wells.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

# where is your MODFLOW-2005 executable?
path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'

# what HUC is Navarro River?
HUC10 = 1801010804

# check if model workspace exists; create if not
modelname = 'Navarro-SteadyState-WithPumping'
model_ws = os.path.join('modflow', modelname)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

if not os.path.isdir(os.path.join(model_ws, 'output')):
    os.makedirs(os.path.join(model_ws, 'output'))

# load existing model
modelname_ss = 'Navarro-SteadyState'
mf = flopy.modflow.Modflow.load(modelname_ss+'.nam', 
                                model_ws=os.path.join('modflow', modelname_ss), 
                                verbose=False,
                                check=False, 
                                exe_name=path2mf)

# change model name and workspace
mf.change_model_ws(model_ws, reset_external=True)

# some discretization info
nlay = mf.dis.nlay
delr = mf.dis.delr[0]
delc = mf.dis.delc[0]
cells_in_navarro = 13037  # this is from R script: MODFLOW_Navarro_InputPrepData.R
area_navarro = cells_in_navarro*delr*delc # [m2]

# extract boundary conditions
ibound = mf.bas6.ibound[0,:,:]

## load river boundary condition info
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')

# figure out which reaches are in the Navarro River HUC10
iriv['Navarro'] = False
iriv.loc[(iriv['HUC'] >= HUC10*100) & (iriv['HUC'] < (HUC10+1)*100), 'Navarro'] = True

# create WEL package, but don't pump anything
wel = flopy.modflow.mfwel.ModflowWel(mf, stress_period_data={0: [0,50,50,0]},
                                     ipakcb=71, filenames=[modelname_ss+'.wel', modelname_ss+'.wel.out'])

# write input datasets
mf.write_input()

# run model without pumping
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## look at output
# figure out timestep
time = mf.dis.perlen[0]

# load RIV output
rivout_nopump = bf.CellBudgetFile(os.path.join(model_ws, modelname_ss+'.riv.out'), verbose=True)
rivout_3D_nopump = rivout_nopump.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D_nopump[0][iriv['lay'],iriv['row'],iriv['col']]
rivout_nopump.close()

# net leakage, calculated as (Infiltration - Discharge)
# units are [m3/d]
# leakage convention: negative value = gaining stream
leakage_mm_yr_nopump = 1000*365*sum(iriv.loc[iriv['Navarro'], 'leakage'])/area_navarro

#### cycle through pumping wells ####
## load wel boundary condition info
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# define pumping rate
Qw = -5000  # [m3/d]
Qw_mm_yr = 1000*365*Qw/area_navarro

# set up empty column
iwel['leakage_mm_yr'] = np.nan
for w in range(0,iwel.shape[0]):
    # set up stress period data
    wel.stress_period_data = {0: [iwel['lay'][w], iwel['row'][w], iwel['col'][w], Qw]}

    # write input and run model
    mf.write_input()
    success, mfoutput = mf.run_model(silent=True)
    if not success:
        raise Exception('MODFLOW did not terminate normally:', w)

    # load RIV output
    rivout_pump = bf.CellBudgetFile(os.path.join(model_ws, modelname_ss+'.riv.out'), verbose=False)
    rivout_3D_pump = rivout_pump.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
    leakage_pump = rivout_3D_pump[0][iriv['lay'],iriv['row'],iriv['col']].data
    rivout_pump.close()

    iwel['leakage_mm_yr'][w] = 1000*365*sum(leakage_pump*iriv['Navarro'])/area_navarro
    
    # status update
    print(w, ' complete')

# calculate depletion
iwel['depletion_prc'] = (leakage_mm_yr_nopump - iwel['leakage_mm_yr'])/Qw_mm_yr

# save output
np.savetxt(os.path.join(model_ws, 'output', 'iwel_out.txt'), iwel)

### code snippets for troubleshooting


# open list file
mfl = flopy.utils.MfListBudget(os.path.join(model_ws,modelname_ss+'.list'))
df_flux_pump, df_vol_pump = mfl.get_dataframes()

# load WEL output
wel_out = bf.CellBudgetFile(os.path.join(model_ws, modelname_ss+'.wel.out'), verbose=True)
wel_out_data = wel_out.get_data(totim=time, text='WELLS')
wel_out.close()