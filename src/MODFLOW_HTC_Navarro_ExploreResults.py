## MODFLOW_HTC_Navarro_ExploreResults.py

import os
import numpy as np
import flopy
import pandas as pd
import flopy.utils.binaryfile as bf
import platform
import copy
import matplotlib.pyplot as plt

# estimated area of navarro
area_domain = 222511*100*100  # [m2]

# what stream BC and modflow version?
stream_BC = 'SFR'  # 'RIV' or 'SFR'
modflow_v = 'mfnwt'  # 'mfnwt' or 'mf2005'
timeType = 'Transient'

# define the folder and the prefix
dir_runs = os.path.join('modflow', 'HTC', 'Navarro', timeType, stream_BC, modflow_v)
prefix_runs = 'mf'
modelname = 'Navarro-'+timeType

sp = np.arange(0,24)

## load budgets
# 0 = no pumping
mfl_NoPump = flopy.utils.MfListBudget(os.path.join(dir_runs, prefix_runs+str(0), modelname+".list"))
df_flux_NoPump, df_vol_NoPump = mfl_NoPump.get_dataframes()

plt.plot(np.arange(0,df_flux_NoPump.shape[0]), df_flux_NoPump['IN-OUT'], '-o')

plt.plot(sp, df_flux_NoPump['IN-OUT'], '-o')

plt.plot(sp, df_flux_NoPump['RECHARGE_IN']/area_domain, '-o')
1000*np.mean(df_flux_NoPump['RECHARGE_IN'])*365/area_domain  # [mm/yr]

riv_net_NoPump = (df_flux_NoPump['RIVER_LEAKAGE_IN']-df_flux_NoPump['RIVER_LEAKAGE_OUT'])
plt.plot(sp, 1000*riv_net_NoPump/area_domain)
np.mean(1000*riv_net_NoPump/area_domain)*365  # [mm/yr]

mnw_net_NoPump = (df_flux_NoPump['MNW2_IN']-df_flux_NoPump['MNW2_OUT'])
plt.plot(sp, mnw_net_NoPump)

plt.plot(sp, df_vol_NoPump['RECHARGE_IN']/area_domain, '-o')

# pumping
w_Pump = 151  # what well to look at
mfl_Pump = flopy.utils.MfListBudget(os.path.join(dir_runs, prefix_runs+str(w_Pump), modelname+".list"))
df_flux_Pump, df_vol_Pump = mfl_Pump.get_dataframes()

plt.plot(sp, 1000*df_flux_Pump['RECHARGE_IN']/area_domain)
1000*np.mean(df_flux_Pump['RECHARGE_IN'])*365/area_domain  # [mm/yr]

riv_net_Pump = (df_flux_Pump['RIVER_LEAKAGE_IN']-df_flux_Pump['RIVER_LEAKAGE_OUT'])
plt.plot(sp, 1000*riv_net_Pump/area_domain)
np.mean(1000*riv_net_Pump/area_domain)*365  # [mm/yr]

mnw_net_Pump = (df_flux_Pump['MNW2_IN']-df_flux_Pump['MNW2_OUT'])
plt.plot(sp, mnw_net_Pump)

plt.plot(sp, df_flux_Pump['IN-OUT'], '-o')


# differences
riv_net_diff = riv_net_NoPump - riv_net_Pump
plt.plot(sp, riv_net_diff, '-o')

mnw_net_diff = mnw_net_NoPump - mnw_net_Pump
plt.plot(sp, mnw_net_diff, '-o')

CHB_net_diff = (df_flux_NoPump['CONSTANT_HEAD_IN']-df_flux_NoPump['CONSTANT_HEAD_OUT']) - (df_flux_Pump['CONSTANT_HEAD_IN']-df_flux_Pump['CONSTANT_HEAD_OUT'])
storage_net_diff = (df_flux_NoPump['STORAGE_IN']-df_flux_NoPump['STORAGE_OUT']) - (df_flux_Pump['STORAGE_IN']-df_flux_Pump['STORAGE_OUT'])
recharge_net_diff = (df_flux_NoPump['RECHARGE_IN']-df_flux_NoPump['RECHARGE_OUT']) - (df_flux_Pump['RECHARGE_IN']-df_flux_Pump['RECHARGE_OUT'])
error_net_diff = df_flux_NoPump['IN-OUT'] - df_flux_Pump['IN-OUT']

capture_prc = riv_net_diff/mnw_net_diff

plt.plot(sp, capture_prc, '-o')

## look at riv
iriv = pd.read_table(os.path.join('modflow', 'input', 'iriv.txt'), delimiter=' ')
iriv_ReachData = pd.read_table(os.path.join('modflow', 'input', 'iriv_ReachData.txt'), delimiter=' ')

time=630

# well 0
rivout_NoPump = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(0), modelname+'.riv.out'), verbose=False)
mnwout_NoPump = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(0), modelname+'.mnw2.out'), verbose=False)

iriv_copy = iriv.copy()

rivout_NoPump_3D = rivout_NoPump.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv_copy['leakage'] = rivout_NoPump_3D[0][iriv['lay'],iriv['row'],iriv['col']]
iriv_NoPump_merge = pd.merge(iriv_ReachData[['SegNum', 'row', 'col','seg_proportion']], iriv_copy[['row', 'col', 'leakage']],
                    how='left', on=['row','col'])
            
iriv_NoPump_merge['leakage'] = iriv_NoPump_merge['leakage']*iriv_NoPump_merge['seg_proportion']        
iriv_NoPump_out = iriv_NoPump_merge.groupby('SegNum', as_index=False).agg({'leakage': 'sum'})
            
mnwout_NoPump_data = mnwout_NoPump.get_data(totim=time, text='MNW2', full3D=False)
iriv_NoPump_out['MNW_net'] = sum(mnwout_NoPump_data[0]['q'])

rivout_NoPump.close()
mnwout_NoPump.close()

# well 226
rivout_Pump = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(226), modelname+'.riv.out'), verbose=False)
mnwout_Pump = bf.CellBudgetFile(os.path.join(dir_runs, prefix_runs+str(226), modelname+'.mnw2.out'), verbose=False)

iriv_copy = iriv.copy()
            
rivout_Pump_3D = rivout_Pump.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv_copy['leakage'] = rivout_Pump_3D[0][iriv['lay'],iriv['row'],iriv['col']]
iriv_Pump_merge = pd.merge(iriv_ReachData[['SegNum', 'row', 'col','seg_proportion']], iriv_copy[['row', 'col', 'leakage']],
                    how='left', on=['row','col'])
            
iriv_Pump_merge['leakage'] = iriv_Pump_merge['leakage']*iriv_Pump_merge['seg_proportion']        
iriv_Pump_out = iriv_Pump_merge.groupby('SegNum', as_index=False).agg({'leakage': 'sum'})
            
mnwout_Pump_data = mnwout_Pump.get_data(totim=time, text='MNW2', full3D=False)
iriv_Pump_out['MNW_net'] = sum(mnwout_Pump_data[0]['q'])

rivout_Pump.close()
mnwout_Pump.close()

# difference
leakage_diff = iriv_NoPump_out['leakage'] - iriv_Pump_out['leakage']
MNW_diff = iriv_Pump_out['MNW_net'] - iriv_NoPump_out['MNW_net']
depletion_prc = leakage_diff/MNW_diff

plt.hist(depletion_prc)
plt.hist(leakage_diff)

sum(leakage_diff)
sum(MNW_diff)
sum(depletion_prc)
max(depletion_prc)