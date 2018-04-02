import os
import numpy as np
import flopy.utils.binaryfile as bf
import pandas as pd

## set things up
modelname = 'Navarro-SteadyState'
model_ws = os.path.join('modflow', modelname)
time=1

## load RIV input
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')
iriv_ReachData = pd.read_table(os.path.join('modflow', 'input', 'iriv_ReachData.txt'), delimiter=' ')

## load RIV output
rivout = bf.CellBudgetFile(os.path.join(model_ws, modelname+'.riv.out'), verbose=True)
rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
rivout.close()

## join leakage to reach_data
iriv_ReachData = pd.merge(iriv_ReachData[['SegNum', 'row', 'col']], iriv[['row', 'col', 'leakage']],
                          how='left', on=['row','col'])
                      
## summarize by segment number
iriv_out = iriv_ReachData.groupby('SegNum', as_index=False).agg({'leakage': 'sum'})

## save RIV output
np.savetxt(os.path.join(model_ws, 'RIV_leakage.txt'), iriv_out)