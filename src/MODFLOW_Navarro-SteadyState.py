## Navarro-SteadyState.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import flopy.utils.sfroutputfile as sf

# where is your MODFLOW-2005 executable?
path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'

# what HUC is Navarro River?
HUC10 = 1801010804

# do you want to use RIV or SFR2 for surface water BC features?
surfWatBC = "SFR2"

# check if model workspace exists; create if not
modelname = 'Navarro-SteadyState'
model_ws = os.path.join('modflow', modelname)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

if not os.path.isdir(os.path.join(model_ws, 'output')):
    os.makedirs(os.path.join(model_ws, 'output'))

# Assign name and create modflow model object
mf = flopy.modflow.Modflow(modelname, exe_name=path2mf, 
                           model_ws=model_ws)

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
spd = {(0, 0): ['save head', 'save budget', 'save drawdown', 'print head', 'print budget', 'print drawdown']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

if (surfWatBC=="RIV"):
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
    iriv.loc[(iriv['HUC'] >= HUC10*100) & (iriv['HUC'] < (HUC10+1)*100), 'Navarro'] = True
    
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
    riv = flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd, ipakcb=61, 
                                   filenames=[modelname+'.riv', modelname+'.riv.out'])
            
elif (surfWatBC=="SFR2"):
    ## SFR2 boundary condition   
    # constant domain parameters
    depth = 4  # river depth?
    riverbed_K = hk/10

    # load R output
    isfr_ReachData = pd.read_table(os.path.join('modflow', 'input', 'isfr_ReachData.txt'), delimiter=' ')
    isfr_SegmentData = pd.read_table(os.path.join('modflow', 'input', 'isfr_SegmentData.txt'), delimiter=' ')

    # set up stream reach data (Dataset 2)
    reach_data = np.array(
              [(0, isfr_ReachData['row'][0], isfr_ReachData['col'][0], isfr_ReachData['SFR_NSEG'][0], 
                isfr_ReachData['SFR_IREACH'][0], isfr_ReachData['length_m'][0], 
                isfr_ReachData['elev_m'][0]-depth,isfr_ReachData['SLOPE'][0], 
                1.0, riverbed_K)], 
    dtype=[('k', '<f8'), ('i', '<f8'), ('j', '<f8'), ('iseg', '<f8'), 
           ('ireach', '<f8'), ('rchlen', '<f8'),
           ('strtop', '<f8'), ('slope', '<f8'), ('strthick', '<f8'), ('strhc1', '<f8')])
    for r in range(1,isfr_ReachData.shape[0]):
        reach_data=np.vstack([reach_data, 
                              np.array(
              [(0, isfr_ReachData['row'][r], isfr_ReachData['col'][r], isfr_ReachData['SFR_NSEG'][r], 
                isfr_ReachData['SFR_IREACH'][r], isfr_ReachData['length_m'][r], 
                isfr_ReachData['elev_m'][r]-depth,isfr_ReachData['SLOPE'][r], 
                1.0, riverbed_K)], 
                dtype=[('k', '<f8'), ('i', '<f8'), ('j', '<f8'), ('iseg', '<f8'), ('ireach', '<f8'), ('rchlen', '<f8'),
                       ('strtop', '<f8'), ('slope', '<f8'), ('strthick', '<f8'), ('strhc1', '<f8')])
                       ])
    reach_data=reach_data[:,0]
    
    # segment data (Dataset 6a-c)
    seg_data_array = np.array(
              [(isfr_SegmentData['SFR_NSEG'][0], 1, isfr_SegmentData['SFR_OUTSEG'][0], 
                0, 0, 0, 0, 0, 0.03, 3, 3)], 
    dtype=[('nseg', '<f8'), ('icalc', '<f8'), ('outseg', '<f8'), ('iupseg', '<f8'), 
           ('flow', '<f8'), ('runoff', '<f8'), ('etsw', '<f8'), ('pptsw', '<f8'), ('roughch', '<f8'),
           ('width1', '<f8'), ('width2', '<f8')])
    for s in range(1,isfr_SegmentData.shape[0]):
        seg_data_array=np.vstack([seg_data_array, 
                              np.array(
              [(isfr_SegmentData['SFR_NSEG'][s], 1, isfr_SegmentData['SFR_OUTSEG'][s], 
                0, 0, 0, 0, 0, 0.03, 3, 3)], 
                dtype=[('nseg', '<f8'), ('icalc', '<f8'), ('outseg', '<f8'), ('iupseg', '<f8'), 
                       ('flow', '<f8'), ('runoff', '<f8'), ('etsw', '<f8'), ('pptsw', '<f8'), 
                       ('roughch', '<f8'), ('width1', '<f8'), ('width2', '<f8')])
             ])
    segment_data = {0: seg_data_array[:,0]}
    
    # constants (dataset 1c)
    nstrm = -len(reach_data) # number of reaches  # negative value so no stream parameters are needed
    nss = len(seg_data_array) # number of segments
    nsfrpar = 0 # number of parameters (not supported)
    nparseg = 0
    const = 86400    # constant for manning's equation, units of m3/d
    dleak = 0.01 # closure tolerance for stream stage computation
    ipakcb = 53 # ISTCB1= flag for writing SFR output to cell-by-cell budget (on unit 53)
    istcb2 = 81 # flag for writing SFR output to text file
    isfropt = 1  # no vertical unsat flow beneath streams
    irtflg = 0
    
    # dataset 5
    dataset_5 = {0: [nss, 0, 0]} # [itmp, irdflag, iptflag]
    
    sfr = flopy.modflow.ModflowSfr2(mf, nstrm=nstrm, nss=nss, const=const, 
                                    dleak=dleak, ipakcb=ipakcb, istcb2=istcb2, 
                                    reach_data=reach_data,
                                    segment_data=segment_data,
                                    isfropt=isfropt,
                                    irtflg=irtflg,
                                    dataset_5=dataset_5,
                                    unit_number=16)
             

print('Using ', surfWatBC, ' for stream features')

## write inputs and run model
# write input datasets
mf.write_input()

# run model
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#### plot results ####
## look at output
# figure out timestep
time = perlen[0]

# open list file
mfl = flopy.utils.MfListBudget(os.path.join(mf.model_ws,modelname+'.list'))
df_flux, df_vol = mfl.get_dataframes()

# load RIV output
rivout = bf.CellBudgetFile(os.path.join(model_ws, modelname+'.riv.out'), verbose=True)
rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
rivout.close()

# net leakage, calculated as (Infiltration - Discharge)
# units are [m3/d]
# leakage convention: negative value = gaining stream
leakage_mm_yr = 1000*365*sum(iriv.loc[iriv['Navarro'], 'leakage'])/area_navarro

# Create the headfile object
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')

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

## save data to plot in R
# head and water table depth
np.savetxt(os.path.join(model_ws, 'output', 'head.txt'), head[0,:,:])
np.savetxt(os.path.join(model_ws, 'output', 'wtd.txt'), wtd)

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

plt.savefig(os.path.join(model_ws, 'output', 'head+WTD'), 
            bbox_inches='tight', dpi=300)
            
## look at input
# plot the top of the model
mf.dis.plot()

mf.lpf.hk.plot(colorbar=True)
mf.rch.rech.plot(colorbar=True)
mf.riv.plot()
