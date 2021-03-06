## Navarro-SteadyState.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
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

# set up your model
modelname = 'Navarro-SteadyState'
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
model_ws = os.path.join('modflow', modelname, stream_BC, modflow_v)
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow(modelname, exe_name=path2mf, 
                           model_ws=model_ws, version=modflow_v)

## Set up DIS and BAS
# read in text output from R for ibound
ibound = np.array(pd.read_csv(os.path.join('modflow', 'input', 'ibound.txt'),
                              header=None, delim_whitespace=True), dtype=np.int32)
                              
# discretization (space) - these should be the same as in your R script                         
nlay = 5
nrow = ibound.shape[0]
ncol = ibound.shape[1]
delr = 100
delc = 100
delv = 20

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

# set up bottom of each layer
botm = np.empty([nlay, nrow, ncol])
botm[0,:,:] = ztop-delv
for l in np.arange(1,(nlay-1)):
    botm[l,:,:] = botm[(l-1),:,:]-delv

# for lowermost layer, set botm equal to zbot everywhere
botm[(nlay-1),:,:] = zbot

# make MODFLOW objects
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, 
                               delr=delr, delc=delc,
                               top=ztop, botm=botm,
                               nper=nper, perlen=perlen, 
                               nstp=nstp, steady=steady)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=0)

## set up flow properties and solver depending on version of MODFLOW
hk = 1e-12*1e7*86400   # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
layvka = 1             # if layvka != 0, Kv = Kh/vka
vka = 10.              # anisotropy
laytyp = 1             # layer type
tol_head = 1e-1
if (modflow_v=='mf2005'):
    lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, layvka=layvka, laytyp=laytyp)
    pcg = flopy.modflow.ModflowPcg(mf, hclose=tol_head, rclose=tol_head)
elif (modflow_v=='mfnwt'):
    upw = flopy.modflow.ModflowUpw(mf, hk=hk, vka=vka, layvka=layvka, laytyp=laytyp)
    nwt = flopy.modflow.ModflowNwt(mf, headtol=tol_head, linmeth=2, options='COMPLEX')

## recharge
# long-term average baseflow is 150 mm/yr
rchrate = 150/(1000*365)  # [mm/yr] --> [m/d]
# total recharge = rchrate*np.sum(ibound == 1)*delr*delc = 914428.767 m3/d over domain
rch = flopy.modflow.ModflowRch(mf, rech=rchrate, nrchop=3)

## output control
spd = {(0, 0): ['save head', 'save budget', 'save drawdown', 'print head', 'print budget', 'print drawdown']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

## stream boundary condition: RIV or SFR  
# constant domain parameters
depth = 5  # river depth?
riverbed_K = hk/10
riverbed_thickness = 1

if (stream_BC=='RIV'):
    ## river boundary condition
    iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')
    
    # estimate conductance based on: river width, river length, riverbed thickness, riverbed K
    iriv['cond'] = round(riverbed_K*iriv['width_m']*iriv['totalLength_m']*riverbed_thickness)   # river bottom conductance? 
    
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

if (stream_BC=='SFR'):
    # load R output
    isfr_ReachData = pd.read_table(os.path.join('modflow', 'input', 'isfr_ReachData.txt'), delimiter=' ')
    isfr_SegmentData = pd.read_table(os.path.join('modflow', 'input', 'isfr_SegmentData.txt'), delimiter=' ')
    
    # set up stream reach data (Dataset 2)
    reach_data = np.array(
              [(0, isfr_ReachData['row'][0], isfr_ReachData['col'][0], isfr_ReachData['SFR_NSEG'][0], 
                isfr_ReachData['SFR_IREACH'][0], isfr_ReachData['length_m'][0], 
                isfr_ReachData['elev_m_min'][0]-depth,isfr_ReachData['SLOPE'][0], 
                riverbed_thickness, riverbed_K)], 
    dtype=[('k', '<f8'), ('i', '<f8'), ('j', '<f8'), ('iseg', '<f8'), 
           ('ireach', '<f8'), ('rchlen', '<f8'),
           ('strtop', '<f8'), ('slope', '<f8'),
           ('strthick', '<f8'), ('strhc1', '<f8')])
    for r in range(1,isfr_ReachData.shape[0]):
        reach_data=np.vstack([reach_data, 
                              np.array(
              [(0, isfr_ReachData['row'][r], isfr_ReachData['col'][r], isfr_ReachData['SFR_NSEG'][r], 
                isfr_ReachData['SFR_IREACH'][r], isfr_ReachData['length_m'][r], 
                isfr_ReachData['elev_m_min'][r]-depth,isfr_ReachData['SLOPE'][r],
                riverbed_thickness, riverbed_K)], 
                dtype=[('k', '<f8'), ('i', '<f8'), ('j', '<f8'), ('iseg', '<f8'), ('ireach', '<f8'), ('rchlen', '<f8'),
                       ('strtop', '<f8'), ('slope', '<f8'), ('strthick', '<f8'), ('strhc1', '<f8')])
                       ])
    reach_data=reach_data[:,0]
    
    ## segment data (Dataset 6a-c)
    # icalc=1 --> stream depth calculates using Manning's equation assuming wide rectangular channel
    # width1, width2 only used if icalc=3 so not in this case
    # roughch=Manning's n=0.4 which seemed like a good ballpark for natural channels: 
    #   http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm
    
    # monthly quickflow [mm/mo] from the script Navarro_StreamflowData.R
    quickflow_mo = [119.61, 91.86, 70.45, 27.65, 5.87, 1.64, 0.52, 0.21, 0.21, 2.04, 16.92, 78.91]
    
    seg_data_array = np.array(
              [(isfr_SegmentData['SFR_NSEG'][0], 1, isfr_SegmentData['SFR_OUTSEG'][0], 0, 
                0, isfr_SegmentData['DirectDASqKM'][0]*1000*1000*sum(quickflow_mo)/(1000*365), 0, 0, 
                0.04, isfr_SegmentData['width_m'][0], isfr_SegmentData['width_m'][0])], 
    dtype=[('nseg', '<f8'), ('icalc', '<f8'), ('outseg', '<f8'), ('iupseg', '<f8'), 
           ('flow', '<f8'), ('runoff', '<f8'), ('etsw', '<f8'), ('pptsw', '<f8'), 
           ('roughch', '<f8'), ('width1', '<f8'), ('width2', '<f8')])
    for s in range(1,isfr_SegmentData.shape[0]):
        seg_data_array=np.vstack([seg_data_array, 
                              np.array(
              [(isfr_SegmentData['SFR_NSEG'][s], 1, isfr_SegmentData['SFR_OUTSEG'][s], 0, 
                0, isfr_SegmentData['DirectDASqKM'][s]*1000*1000*sum(quickflow_mo)/(1000*365), 0, 0, 
                0.04, isfr_SegmentData['width_m'][s], isfr_SegmentData['width_m'][s])], 
                dtype=[('nseg', '<f8'), ('icalc', '<f8'), ('outseg', '<f8'), ('iupseg', '<f8'), 
                       ('flow', '<f8'), ('runoff', '<f8'), ('etsw', '<f8'), ('pptsw', '<f8'), 
                       ('roughch', '<f8'), ('width1', '<f8'), ('width2', '<f8')])
             ])
    segment_data = {0: seg_data_array[:,0]}
    
    # constants (dataset 1c)
    nstrm = -len(reach_data)  # number of reaches  # negative value so no stream parameters are needed
    nss = len(seg_data_array) # number of segments
    nsfrpar = 0      # number of parameters (not supported)
    nparseg = 0
    const = 86400    # constant for manning's equation, units of m3/d
    dleak = 0.01   # closure tolerance for stream stage computation
    ipakcb = 53      # ISTCB1= flag for writing SFR output to cell-by-cell budget (on unit 53)
    istcb2 = 81      # flag for writing SFR output to text file
    isfropt = 1      # no vertical unsat flow beneath streams
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

## write inputs and run model (the first time - will be run again with MNW2 package)
mf.write_input()
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## need to figure out steady-state WTE to define wells
# use head output from the SS situation without boreholes
time = perlen[0]
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')
head = h.get_data(totim=time)
wte = pp.get_water_table(head, nodata=bas.hnoflo)

## create MNW2 package
# Based on: https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_mnw2package_example.ipynb
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# define well parameters
losstype = 'THIEM'
pumploc = 0  # 0=intake location assumed to occur above first active node
qlimit = 0   # 0=do not constrain pumping rate based on head
ppflag = 1   # 1=adjust head for partial penetration
pumpcap = 0  # 0=do not adjust pumping rate for changes in lift
rw = 0.25
screen_length = 50  # [m] - will start at SS WTE

start_flag = True
for w in range(0,iwel.shape[0]):
    WellNum = iwel['WellNum'][w]
    wellid = 'Well'+str(WellNum)

    # extract WTE from SS results for top of screen
    row_wel = iwel['row'][w]
    col_wel = iwel['col'][w]
    screen_top = wte[row_wel, col_wel]
    screen_bot = screen_top - screen_length
    
    # figure out which layer well screen starts and ends in
    if (screen_top > mf.dis.botm[0, row_wel, col_wel]):
        screen_top_k = 0
    elif (screen_top > mf.dis.botm[1, row_wel, col_wel]):
        screen_top_k = 1
    elif (screen_top > mf.dis.botm[2, row_wel, col_wel]):
        screen_top_k = 2
    elif (screen_top > mf.dis.botm[3, row_wel, col_wel]):
        screen_top_k = 3
    elif (screen_top > mf.dis.botm[4, row_wel, col_wel]):
        screen_top_k = 4

    if (screen_bot > mf.dis.botm[0, row_wel, col_wel]):
        screen_bot_k = 0
    elif (screen_bot > mf.dis.botm[1, row_wel, col_wel]):
        screen_bot_k = 1
    elif (screen_bot > mf.dis.botm[2, row_wel, col_wel]):
        screen_bot_k = 2
    elif (screen_bot > mf.dis.botm[3, row_wel, col_wel]):
        screen_bot_k = 3
    elif (screen_bot > mf.dis.botm[4, row_wel, col_wel]):
        screen_bot_k = 4 
    
    # build node data
    if (screen_top_k==screen_bot_k):
        # when top and bottom of screen are in same layer, only 1 node necessary
        node_data = pd.DataFrame([[wellid, screen_top_k, row_wel, col_wel, 
                               screen_top, screen_bot, 
                               losstype, pumploc, qlimit, ppflag, pumpcap, rw]], 
                 columns=['wellid', 'k', 'i', 'j', 'ztop', 'zbotm', 'losstype', 
                 'pumploc', 'qlimit', 'ppflag', 'pumpcap', 'rw'])
    else:
        # build data frame with nodes in all intersected layers
        node_data = pd.DataFrame([[wellid, screen_top_k, row_wel, col_wel, 
                                   screen_top,
                                   mf.dis.botm[screen_top_k, row_wel, col_wel], 
                                   losstype, pumploc, qlimit, ppflag, pumpcap, rw]], 
                     columns=['wellid', 'k', 'i', 'j', 'ztop', 'zbotm', 'losstype', 
                     'pumploc', 'qlimit', 'ppflag', 'pumpcap', 'rw'])
        for lay_wel in range(screen_top_k+1, screen_bot_k+1):
            if (lay_wel==screen_bot_k):
                lay_bot = screen_bot
            else:
                lay_bot = mf.dis.botm[lay_wel, row_wel, col_wel]
                
            node_data_lay = pd.DataFrame([[wellid, lay_wel, row_wel, col_wel, 
                                       mf.dis.botm[lay_wel-1, row_wel, col_wel], lay_bot, 
                                       losstype, pumploc, qlimit, ppflag, pumpcap, rw]], 
                         columns=['wellid', 'k', 'i', 'j', 'ztop', 'zbotm', 'losstype', 
                         'pumploc', 'qlimit', 'ppflag', 'pumpcap', 'rw'])
            node_data = node_data.append(node_data_lay)
            
    # set up stress period data
    mnw_spd = pd.DataFrame([[0, wellid, 0]], columns=['per', 'wellid', 'qdes'])
    
    # make master output node dat
    if (start_flag):
        node_data_all = node_data
        mnw_spd_all = mnw_spd
        start_flag = False
    else:
        node_data_all = node_data_all.append(node_data)
        mnw_spd_all = mnw_spd_all.append(mnw_spd)

# convert node data to recarray to work with python
node_data_all = node_data_all.to_records()

# set up stress period data
stress_period_data = {0: mnw_spd_all.to_records()}

# mnw2
n_wel = iwel.shape[0]
mnw2 = flopy.modflow.ModflowMnw2(model=mf, mnwmax=n_wel, 
                                 node_data=node_data_all,
                                 stress_period_data=stress_period_data,
                                 itmp=[n_wel], ipakcb=71,
                                 filenames=[modelname+'.mnw2', modelname+'.mnw2.out'])

## update starting head
head[head <= bas.hnoflo] = 0
bas.strt = head

## rerun model with MNW boreholes (but no pumping)
mf.write_input()
success, mfoutput = mf.run_model()
if not success:
    raise Exception('MODFLOW did not terminate normally.')

## look at budget outputs
mfl = flopy.utils.MfListBudget(os.path.join(model_ws, modelname+".list"))
df_flux, df_vol = mfl.get_dataframes()
print(df_flux)

#### plot results ####
## look at output
# figure out timestep
time = perlen[0]

## head output
# Create the headfile object
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')

# extract data matrix
head = h.get_data(totim=time)
head[head <= bas.hnoflo] = np.nan

# calculate WTD
wte = pp.get_water_table(head, nodata=bas.hnoflo)
wtd = ztop - wte

## save data to plot in R
# head and water table depth
np.savetxt(os.path.join(model_ws, 'wte.csv'), wte, fmt='%6.2f', delimiter=',')
np.savetxt(os.path.join(model_ws, 'wtd.csv'), wtd, fmt='%6.2f', delimiter=',')

## stream output
if (stream_BC=='SFR'):
    ## SFR output
    sfrout = sf.SfrFile(os.path.join(model_ws, modelname+'.sfr.out'))
    sfr_df = sfrout.get_dataframe()
    sfr_df.to_csv(os.path.join(model_ws, 'sfr.csv'), index=False)
