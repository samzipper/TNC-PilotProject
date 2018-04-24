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
import platform

# set up your model
modelname = 'Navarro-SteadyState'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
stream_BC = 'RIV'  # options: 'RIV' 'SFR'

# where is your MODFLOW-2005 executable?
if (modflow_v=='mf2005'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'
    else:
        path2mf = modflow_v
elif (modflow_v=='mfnwt'):
    if platform.system() == 'Windows':
        path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MODFLOW-NWT_1.1.3/bin/MODFLOW-NWT.exe'
    else:
        path2mf = modflow_v

# check if model workspace exists; create if not
model_prefix = 'mf'
model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v, model_prefix+'0')
if not os.path.isdir(model_ws):
    os.makedirs(model_ws)

# Assign name and create modflow model object
mf = flopy.modflow.Modflow.load(modelname+'.nam', 
        exe_name=path2mf, version=modflow_v, model_ws=os.path.join('modflow', modelname, stream_BC, modflow_v))

# update model workspace
mf.change_model_ws(model_ws)

# print info about stream BC package
print('Using ', stream_BC, ' for stream features')
                
## write inputs for no-pumping scenario
mf.write_input()

## copy launch script
shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'), 
        os.path.join(model_ws, 'launch_thisRun.sh'))

### Now: scroll through wells and create WEL input files
## load wel boundary condition info
iwel = pd.read_table(os.path.join('modflow', 'input', 'iwel.txt'), delimiter=' ')

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

# load WTE from SS simulation, which will be used to define top of well screen
wte = np.array(pd.read_csv(os.path.join('modflow', modelname, stream_BC, modflow_v, 'wte.csv'),
                           header=None), dtype=np.float64)

for w in range(0,9): #iwel.shape[0]):
    WellNum = iwel['WellNum'][w]
    wellid = 'Well'+str(WellNum)

    # create output folder
    w_model_ws = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modflow_v, model_prefix+str(WellNum))
    if not os.path.isdir(w_model_ws):
        os.makedirs(w_model_ws)
        mf.model_ws = w_model_ws

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

    # write input (NAM and WEL only)
    mf.write_input(SelPackList=['MNW2'])
        
    # copy launch script
    shutil.copy2(os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, 'launch_thisRun_'+modflow_v+'.sh'),
            os.path.join(w_model_ws, 'launch_thisRun.sh'))


    # copy namefile from template, which points to input package files in mf0 (no pumping) directory
    path_nam_template = os.path.join('modflow', 'HTC', 'Navarro', 'SteadyState', stream_BC, modelname+'_Template_'+modflow_v+'.nam')
    shutil.copy2(path_nam_template, os.path.join(w_model_ws, modelname+'.nam'))

