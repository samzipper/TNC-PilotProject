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
delr = 100
delc = 100
cells_in_navarro = 81533  # this is from R script
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
            isfr_ReachData['elev_m_min'][0]-depth,isfr_ReachData['SLOPE'][0], 
            1.0, riverbed_K)], 
dtype=[('k', '<f8'), ('i', '<f8'), ('j', '<f8'), ('iseg', '<f8'), 
       ('ireach', '<f8'), ('rchlen', '<f8'),
       ('strtop', '<f8'), ('slope', '<f8'), ('strthick', '<f8'), ('strhc1', '<f8')])
for r in range(1,isfr_ReachData.shape[0]):
    reach_data=np.vstack([reach_data, 
                          np.array(
          [(0, isfr_ReachData['row'][r], isfr_ReachData['col'][r], isfr_ReachData['SFR_NSEG'][r], 
            isfr_ReachData['SFR_IREACH'][r], isfr_ReachData['length_m'][r], 
            isfr_ReachData['elev_m_min'][r]-depth,isfr_ReachData['SLOPE'][r], 
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

## look at sfr output
# example: https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_sfrpackage_example.ipynb
sfr.plot(key='iseg', vmin=1, vmax=max(isfr_SegmentData.SFR_NSEG))
sfr.plot(key='strtop', vmin=0, vmax=max(isfr_ReachData.elev_m_min))

# data frame of water balace
sfrout = sf.SfrFile(os.path.join(model_ws, modelname+'.sfr.out'))
sfr_df = sfrout.get_dataframe()
sfr_df.head()

# plot streamflow and stream/aquifer interactions for a segment
inds = sfr_df.segment == 479
ax = sfr_df.ix[inds, ['Qin', 'Qaquifer', 'Qout']].plot(x=sfr_df.reach[inds])
ax.set_ylabel('Flow, in cubic meters per day')
ax.set_xlabel('SFR reach')

# look at stage, model top, and streambed top
streambed_top = mf.sfr.segment_data[0][mf.sfr.segment_data[0].nseg == 479][['elevup', 'elevdn']][0]
sfr_df['model_top'] = mf.dis.top.array[sfr_df.row.values - 1, sfr_df.column.values -1]
fig, ax = plt.subplots()
plt.plot([1, 6], list(streambed_top), label='streambed top')
ax = sfr_df.ix[inds, ['stage', 'model_top']].plot(ax=ax, x=sfr_df.reach[inds])
ax.set_ylabel('Elevation, in feet')
plt.legend()

# get SFR leakage results from cell budget file
bpth = os.path.join(model_ws, modelname+'.cbc')
cbbobj = bf.CellBudgetFile(bpth)
cbbobj.list_records()
sfrleak = cbbobj.get_data(text='  STREAM LEAKAGE')[0]
sfrleak.q[sfrleak.q == 0] = np.nan # remove zero values

# plot of leakage, plan view
im = plt.imshow(sfrleak.q, interpolation='none', cmap='coolwarm')
cb = plt.colorbar(im, label='SFR Leakage, in cubic meters per day');

# plot total streamflow
sfrQ = sfrleak[0].copy()
sfrQ[sfrQ == 0] = np.nan
sfrQ[sfr_df.row.values-1, sfr_df.column.values-1] = sfr_df[['Qin', 'Qout']].mean(axis=1).values
im = plt.imshow(sfrQ, interpolation='none')
plt.colorbar(im, label='Streamflow, in cubic meters per day');

## look at head output
# Create the headfile object
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')

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