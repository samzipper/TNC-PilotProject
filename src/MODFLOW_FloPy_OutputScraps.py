## MODFLOW_FloPyOutputScraps.py
# Random commands for processing, plotting FloPy output

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
sfr_df.to_csv(os.path.join(model_ws, 'output', 'sfr.csv'), index=False)

## process gaging station info
Qoutlet = sfr_df.loc[(sfr_df['segment']==gage_data[0][0]) & (sfr_df['reach']==gage_data[0][1])].Qout

# convert m3/d to mm/yr
Qoutlet_mm = 365*1000*Qoutlet/(gage_data_in.TotDASqKM[0]*1000*1000)


# put data frame into 
Qaquifer = np.zeros((nrow, ncol))
Qaquifer[sfr_df.i, sfr_df.j] = sfr_df.Qout
plt.imshow(np.log10(Qaquifer)); plt.colorbar();

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
sfrleak = cbbobj.get_data(text='  STREAM LEAKAGE', full3D=True)[0]
sfrleak.q[sfrleak.q == 0] = np.nan # remove zero values

# plot of leakage, plan view
im = plt.imshow(sfrleak[0,:,:], interpolation='none', cmap='coolwarm')
cb = plt.colorbar(im, label='SFR Leakage, in cubic meters per day');

# plot total streamflow
sfrQ = sfrleak[0].copy()
sfrQ[sfrQ == 0] = np.nan
sfrQ[sfr_df.row.values-1, sfr_df.column.values-1] = sfr_df[['Qin', 'Qout']].mean(axis=1).values
im = plt.imshow(sfrQ, interpolation='none')
plt.colorbar(im, label='Streamflow, in cubic meters per day');

## load RIV output
rivout = bf.CellBudgetFile(os.path.join(model_ws, modelname+'.riv.out'), verbose=True)
rivout_3D = rivout.get_data(totim=time, text='RIVER LEAKAGE', full3D=True)
iriv['leakage'] = rivout_3D[0][iriv['lay'],iriv['row'],iriv['col']]
rivout.close()

# net leakage, calculated as (Infiltration - Discharge)
# units are [m3/d]
# leakage convention: negative value = gaining stream
leakage_mm_yr = 1000*365*sum(iriv.loc[iriv['Navarro'], 'leakage'])/area_navarro

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
np.savetxt(os.path.join(model_ws, 'head.txt'), head[0,:,:])
np.savetxt(os.path.join(model_ws, 'wtd.txt'), wtd)
            
## look at input
# plot the top of the model
mf.dis.plot()

mf.lpf.hk.plot(colorbar=True)
mf.rch.rech.plot(colorbar=True)