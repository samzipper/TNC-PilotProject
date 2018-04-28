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
# which well do you want to use?
WellNum = 493
wellid = 'Well'+str(WellNum)

# load wel boundary condition info
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

# get index
w = iwel.index[iwel['WellNum']==WellNum].tolist()[0]

# define pumping rate
Qw = -6*100*0.00378541  # [m3/d]  6 gal/plant/day*100 plants*0.00378541 m3/gal

# define well screen length (will start at SS WTE)
screen_length = 50  # [m]

# load existing MNW2 package
mnw2 = mf.Mnw2

# define well parameters
losstype = 'THIEM'
pumploc = 0
qlimit = 0
ppflag = 1
pumpcap = 0
rw = 0.25

# extract WTE from heads, which will be used to define top of well screen
head[head <= mf.bas6.hnoflo] = np.nan
wte = pp.get_water_table(head, nodata=mf.bas6.hnoflo)

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
    

# convert to recarray to work with python
node_data = node_data.to_records()

# set up stress period data
stress_period_data = {0: pd.DataFrame([[0, wellid, Qw]],
                                      columns=['per', 'wellid', 'qdes']).to_records()}   
                                               
mnw2 = flopy.modflow.ModflowMnw2(model=mf, mnwmax=1, 
                             node_data=node_data,
                             stress_period_data=stress_period_data,
                             itmp=[1],)

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
print(df_flux)

#### plot results ####
## look at output
# figure out timestep
time = 1

## head output
# Create the headfile object
h = bf.HeadFile(os.path.join(mf.model_ws, modelname+'.hds'), text='head')

# extract data matrix
head = h.get_data(totim=time)
head[head <= mf.bas6.hnoflo] = np.nan

# calculate WTE and DDN
wte_pumped = pp.get_water_table(head, nodata=mf.bas6.hnoflo)
ddn = wte - wte_pumped

## make plot
plt.imshow(ddn)
plt.colorbar()