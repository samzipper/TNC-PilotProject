## Navarro-SteadyState.py
# This script creates a steady-state groundwater flow model for the 
# Navarro River Watershed in California.
#
# Using default units of ITMUNI=4 (days) and LENUNI=2 (meters)

import os
import numpy as np
import flopy
import shapefile as shp
import matplotlib.pyplot as plt
import pandas as pd


#

# where is your MODFLOW-2005 executable?
path2mf = 'C:/Users/Sam/Dropbox/Work/Models/MODFLOW/MF2005.1_12/bin/mf2005.exe'

# Assign name and create modflow model object
modelname = 'Navarro-SteadyState'
mf = flopy.modflow.Modflow(modelname, exe_name=path2mf)

## Set up BAS
# read in text output from R for ibound
ibound = np.array(pd.read_csv(os.path.join('modflow', 'input', 'ibound.txt'),
                              header=None, delim_whitespace=True), dtype=np.int32)
                              
# set up constants - these should be the same as in your R script
nlay = 1
nrow = ibound.shape[0]
ncol = ibound.shape[1]
delr = 1000
delc = 1000

# define starting head (set to 0 everywhere for now)
strt = np.zeros((nrow,ncol), dtype='float')

# set up bas
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
