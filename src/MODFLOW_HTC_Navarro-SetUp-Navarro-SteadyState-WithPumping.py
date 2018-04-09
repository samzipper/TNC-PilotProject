## MODFLOW_HTC_Navarro-SetUp-SteadyState-WithPumping.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import pandas as pd
import shutil

# where is your MODFLOW-2005 executable?
path2mf = 'mf2005'   # this was compiled with pymake and added to path so should be accessible everywhere

# do you want to use RIV or SFR for stream features?
stream_BC = 'SFR'  # options: 'RIV' 'SFR'

# check if model workspace exists; create if not
model_prefix = 'mf'
modelname = 'Navarro-SteadyState'
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, model_prefix+'0')
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

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
# total recharge = rchrate*np.sum(ibound == 1)*delr*delc = 914428.767 m3/d over domain
rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=3)

## output control
spd = {(0, 0): ['save head', 'save budget', 'save drawdown', 'print head', 'print budget', 'print drawdown']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

## set up stream boundary condition - RIV or SFR
# parameters needed for both packages
depth = 4  # river depth?
riverbed_K = hk/10
river_width = 10
riverbed_thickness = 1

## SFR2 boundary condition
if (stream_BC == 'SFR'):
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

    ## segment data (Dataset 6a-c)
    # width1, width2 only used if icalc=3
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

    sfr = flopy.modflow.ModflowSfr2(mf, nstrm=nstrm, nss=nss, const=const, 
                                    dleak=dleak, ipakcb=ipakcb, istcb2=istcb2, 
                                    reach_data=reach_data,
                                    segment_data=segment_data,
                                    isfropt=isfropt,
                                    irtflg=irtflg,
                                    unit_number=16)

    ## stream gaging station
    # read in output from R
    gage_data_in = pd.read_table(os.path.join('modflow', 'input', 'gage_data.txt'), delimiter=' ')

    # total number of gages
    numgage = gage_data_in.shape[0]

    # set up gage data
    gage_data=[[gage_data_in.SFR_NSEG[0], gage_data_in.SFR_IREACH[0], 
                90+gage_data_in.GageNum[0], 2]]

    gage = flopy.modflow.ModflowGage(mf, numgage=numgage,
                                     gage_data=gage_data, unitnumber=90)

if (stream_BC=='RIV'):
    ## RIV boundary condition
    iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')

    # estimate conductance based on: river width, river length, riverbed thickness, riverbed K
    iriv['cond'] = round(riverbed_K*river_width*iriv['totalLength_m']*riverbed_thickness)   # river bottom conductance? 
    
    # empty list to hold RIV BC data
    riv_list = []

    # populate list
    for r in range(0,iriv.shape[0]):
            riv_list.append([iriv['lay'][r], iriv['row'][r], iriv['col'][r], 
                                     iriv['elev_m_min'][r], iriv['cond'][r], iriv['elev_m_min'][r]-depth])    
            riv_spd = {0: riv_list}

    # make MODFLOW object
    riv = flopy.modflow.ModflowRiv(mf, stress_period_data=riv_spd, ipakcb=61, 
                                       filenames=[modelname+'.riv', modelname+'.riv.out'])

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')
   
## create WEL package, but don't pump anything
wel = flopy.modflow.mfwel.ModflowWel(mf, stress_period_data={0: [0,50,50,0]},
                                     ipakcb=71, filenames=[modelname+'.wel', modelname+'.wel.out'])

## write inputs for no-pumping scenario
mf.write_input()

## copy launch script
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun.sh'), model_ws)

### Now: scroll through wells and create WEL input files
## load wel boundary condition info
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

for w in range(0,iwel.shape[0]):
    WellNum = iwel['WellNum'][w]

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

    # set up stress period data
    wel.stress_period_data = {0: [iwel['lay'][w], iwel['row'][w], iwel['col'][w], Qw]}
    wel.filenames=[modelname+'.wel', modelname+'.wel.out']

    # write input (NAM and WEL only)
    mf.write_input(SelPackList=['WEL'])
        
    # copy launch script
    shutil.copy2(os.path.join(model_ws, 'launch_thisRun.sh'), w_model_ws)

    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modelname+'_Template.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))

