## Navarro-SteadyState.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd

# where is your MODFLOW-2005 executable?
path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'

# Assign name and create modflow model object
modelname = 'Navarro-SteadyState'
mf = flopy.modflow.Modflow(modelname, exe_name=path2mf, 
                           model_ws=os.path.join('modflow', '.'))

# what HUC is Navarro River?
HUC10 = 1801010804

## Set up DIS and BAS
# read in text output from R for ibound
ibound = np.array(pd.read_csv(os.path.join('modflow', 'input', 'ibound.txt'),
                              header=None, delim_whitespace=True), dtype=np.int32)
                              
# discretization (space) - these should be the same as in your R script                         
nlay = 1
nrow = ibound.shape[0]
ncol = ibound.shape[1]
delr = 250
delc = 250
cells_in_navarro = 13037  # this is from R script
area_navarro = cells_in_navarro*delr*delc # [m2]

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
hk = (1e-5)*86400     # horizontal K, convert [m/s] to [m/d]
vka = 1.    # anisotropy
sy = 0.2    # specific yield
ss = 1e-5  # specific storage
laytyp = 1  # layer type

# make MODFLOW objects
lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp)
pcg = flopy.modflow.ModflowPcg(mf, hclose=1e-2, rclose=1e-2)

## recharge
# long-term average baseflow is 150 mm/yr
rchrate = 150/(1000*365)  # [mm/yr] --> [m/d]
rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=3)

## output control
spd = {(0, 0): ['save head', 'save budget', 'save drawdown', 'print budget']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

## river boundary condition
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')

# constant domain parameters
depth = 4  # river depth?

# estimate conductance based on: river width, river length, riverbed thickness, riverbed K
river_width = 10
riverbed_thickness = 1
riverbed_K = hk/10
iriv['cond'] = round(riverbed_K*river_width*iriv['length_m']*riverbed_thickness)   # river bottom conductance? 

# figure out which reaches are in the Navarro River HUC10
iriv['Navarro'] = False
iriv['Navarro'][(iriv['HUC'] >= HUC10*100) & (iriv['HUC'] < (HUC10+1)*100)] = True

# empty list to hold stress period data
riv_list = []

# populate list
for r in range(0,iriv.shape[0]):
#    riv_list.append([iriv['lay'][r], iriv['row'][r], iriv['col'][r], 
#                     0, cond, 0-depth]) 
    riv_list.append([iriv['lay'][r], iriv['row'][r], iriv['col'][r], 
                     iriv['stage_m'][r], iriv['cond'][r], iriv['stage_m'][r]-depth])    
riv_spd = {0: riv_list}

# make MODFLOW object
riv = flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd, ipakcb=1, 
                               filenames=[modelname+'.riv', modelname+'.riv.out'])

## write inputs and run model
# write input datasets
mf.write_input()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#### plot results ####
# Imports
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import flopy.utils.sfroutputfile as sf
from flopy.utils.postprocessing import get_transmissivities, get_water_table, get_gradients, get_saturated_thickness

## look at input
# plot the top of the model
mf.dis.plot()

mf.lpf.hk.plot(colorbar=True)
mf.rch.rech.plot(colorbar=True)
mf.riv.plot()

## look at head output
# figure out timestep
time = perlen[0]

# Create the headfile object
h = bf.HeadFile(os.path.join('modflow', modelname+'.hds'), text='head')

h.plot(totim=time, colorbar=True, contour=True,
       masked_values=[bas.hnoflo])

# extract data matrix
head = h.get_data(totim=time)
head[head <= bas.hnoflo] = np.nan
plt.imshow(head[0,:,:]); plt.colorbar()

# calculate WTD
wtd = ztop - head[0,:,:]
plt.imshow(wtd); plt.colorbar()
plt.imshow(wtd<0)

# load RIV output
rivout = bf.CellBudgetFile(os.path.join('modflow', modelname+'.riv.out'), verbose=True)
rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]

# total leakage
# units are [m3/d]? convert to [mm/yr]
# leakage convention: negative value = groundwater discharge to stream
leakage_mm_yr = 1000*365*sum(iriv['leakage'][iriv['Navarro']])/area_navarro

rivout_dat[0]['node']

nodes = 3282, 3010, 3283
q = -1254.076, 352.000, 640


## save data to plot in R
# head and water table depth
np.savetxt(os.path.join('modflow', 'output', 'head.txt'), head[0,:,:])
np.savetxt(os.path.join('modflow', 'output', 'wtd.txt'), wtd)

## make some plots in Python if you want
plt.subplot(1,2,1)
plt.imshow(head[0,:,:])
plt.title('Head [m]')
plt.viridis()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar()

plt.subplot(1,2,2)
plt.imshow(wtd)
plt.title('Water Table Depth [m]')
plt.viridis()
plt.colorbar()

plt.savefig(os.path.join('modflow', 'output', 'head+WTD'), 
            bbox_inches='tight', dpi=300)

mf.dis.top.plot(colorbar=True, filename_base=)
mf.dis.botm.plot(colorbar=True, filename_base=os.path.join('modflow', 'output', 'dis_botm'))

### experimental
head[head <= -9999] = -100
plt.imshow(head[0,:,:])
plt.colorbar()

# water table
wt = get_water_table(head, mf, nodata=bas.hnoflo)
plt.imshow(wt)
plt.colorbar(label='Water Table')

# saturated thickness
st = get_saturated_thickness(head, mf, nodata=head[0,0,0])
plt.imshow(st)
plt.colorbar(label='Saturated Thickness')