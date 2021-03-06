#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
infile01 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
infile21 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601_2/ctm.01.20160116.nc')
infile02 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602/ctm.01.20160201.nc')
infile22 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602_2/ctm.01.20160216.nc')
infile03 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603/ctm.01.20160301.nc')
infile23 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603_2/ctm.01.20160316.nc')
infile04 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604/ctm.01.20160401.nc')
infile24 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604_2/ctm.01.20160416.nc')
mopfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.10.20160101.nc')
mopfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601_2/ctm.10.20160116.nc')
mopfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602/ctm.10.20160201.nc')
mopfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602_2/ctm.10.20160216.nc')
mopfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603/ctm.10.20160301.nc')
mopfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603_2/ctm.10.20160316.nc')
mopfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604/ctm.10.20160401.nc')
mopfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604_2/ctm.10.20160416.nc')
omhfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/ctm.10.20160101.nc')
omhfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/ctm.10.20160101.nc')
omhfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/ctm.10.20160116.nc')
omhfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/ctm.10.20160201.nc')
omhfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/ctm.10.20160216.nc')
omhfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/ctm.10.20160301.nc')
omhfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/ctm.10.20160316.nc')
omhfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/ctm.10.20160401.nc')
omhfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/ctm.10.20160416.nc')
iasifile01 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1601/ctm.15.20160101.nc')
iasifile21 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1601_2/ctm.15.20160116.nc')
iasifile02 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1602/ctm.15.20160201.nc')
iasifile22 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1602_2/ctm.15.20160216.nc')
iasifile03 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1603/ctm.15.20160301.nc')
iasifile23 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1603_2/ctm.15.20160316.nc')
iasifile04 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1604/ctm.15.20160401.nc')
iasifile24 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1604_2/ctm.15.20160416.nc')
omifile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601/ctm.10.20160101.nc')
omifile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601_2/ctm.10.20160116.nc')
omifile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602/ctm.10.20160201.nc')
omifile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602_2/ctm.10.20160216.nc')
omifile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603/ctm.10.20160301.nc')
omifile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603_2/ctm.10.20160316.nc')
omifile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604/ctm.10.20160401.nc')
omifile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604_2/ctm.10.20160416.nc')
allfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_06/ctm.20.20160101.nc')
allfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_26/ctm.15.20160116.nc')
allfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_06/ctm.20.20160201.nc')
allfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_26/ctm.20.20160216.nc')
allfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_06/ctm.15.20160301.nc')
allfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_26/ctm.20.20160316.nc')
allfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_06/ctm.15.20160401.nc')
allfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_26/ctm.15.20160416.nc')
infile_am = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1608_AM/ctm.00.20160801.nc')

#XOX-BIOB__XCH2O or anTHSRCE__XCH2O or XOX-an-S__XCH2O or xOX-LI-S__xCH2O
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][0,:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len3 = shape(omhfile01.variables["CHEM-L_S__OH"][:,:,:])
len4 = 6
len5 = 8
len6 = len(infile01.variables["CHEM-L_S__OH"][:,0,0])
OH = np.zeros((len4,len5,len6,len2,len3))
AM = infile_am.variables["BXHGHT-S__AD"][0:38,:,:]

OH[0,0,:,:,:] = infile01.variables["CHEM-L_S__OH"][:,:,:]
OH[0,1,:,:,:] = infile21.variables["CHEM-L_S__OH"][:,:,:]
OH[0,2,:,:,:] = infile02.variables["CHEM-L_S__OH"][:,:,:]
OH[0,3,:,:,:] = infile21.variables["CHEM-L_S__OH"][:,:,:]
OH[0,4,:,:,:] = infile03.variables["CHEM-L_S__OH"][:,:,:]
OH[0,5,:,:,:] = infile23.variables["CHEM-L_S__OH"][:,:,:]
OH[0,6,:,:,:] = infile04.variables["CHEM-L_S__OH"][:,:,:]
OH[0,7,:,:,:] = infile24.variables["CHEM-L_S__OH"][:,:,:]

OH[1,0,:,:,:] = omhfile01.variables["CHEM-L_S__OH"][:,:,:]
OH[1,1,:,:,:] = omhfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[1,2,:,:,:] = omhfile02.variables["CHEM-L_S__OH"][:,:,:]
OH[1,3,:,:,:] = omhfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[1,4,:,:,:] = omhfile03.variables["CHEM-L_S__OH"][:,:,:]
OH[1,5,:,:,:] = omhfile23.variables["CHEM-L_S__OH"][:,:,:]
OH[1,6,:,:,:] = omhfile04.variables["CHEM-L_S__OH"][:,:,:]
OH[1,7,:,:,:] = omhfile24.variables["CHEM-L_S__OH"][:,:,:]

OH[2,0,:,:,:] = iasifile01.variables["CHEM-L_S__OH"][:,:,:]
OH[2,1,:,:,:] = iasifile21.variables["CHEM-L_S__OH"][:,:,:]
OH[2,2,:,:,:] = iasifile02.variables["CHEM-L_S__OH"][:,:,:]
OH[2,3,:,:,:] = iasifile21.variables["CHEM-L_S__OH"][:,:,:]
OH[2,4,:,:,:] = iasifile03.variables["CHEM-L_S__OH"][:,:,:]
OH[2,5,:,:,:] = iasifile23.variables["CHEM-L_S__OH"][:,:,:]
OH[2,6,:,:,:] = iasifile04.variables["CHEM-L_S__OH"][:,:,:]
OH[2,7,:,:,:] = iasifile24.variables["CHEM-L_S__OH"][:,:,:]

OH[3,0,:,:,:] = mopfile01.variables["CHEM-L_S__OH"][:,:,:]
OH[3,1,:,:,:] = mopfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[3,2,:,:,:] = mopfile02.variables["CHEM-L_S__OH"][:,:,:]
OH[3,3,:,:,:] = mopfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[3,4,:,:,:] = mopfile03.variables["CHEM-L_S__OH"][:,:,:]
OH[3,5,:,:,:] = mopfile23.variables["CHEM-L_S__OH"][:,:,:]
OH[3,6,:,:,:] = mopfile04.variables["CHEM-L_S__OH"][:,:,:]
OH[3,7,:,:,:] = mopfile24.variables["CHEM-L_S__OH"][:,:,:]

OH[4,0,:,:,:] = omifile01.variables["CHEM-L_S__OH"][:,:,:]
OH[4,1,:,:,:] = omifile21.variables["CHEM-L_S__OH"][:,:,:]
OH[4,2,:,:,:] = omifile02.variables["CHEM-L_S__OH"][:,:,:]
OH[4,3,:,:,:] = omifile21.variables["CHEM-L_S__OH"][:,:,:]
OH[4,4,:,:,:] = omifile03.variables["CHEM-L_S__OH"][:,:,:]
OH[4,5,:,:,:] = omifile23.variables["CHEM-L_S__OH"][:,:,:]
OH[4,6,:,:,:] = omifile04.variables["CHEM-L_S__OH"][:,:,:]
OH[4,7,:,:,:] = omifile24.variables["CHEM-L_S__OH"][:,:,:]

OH[5,0,:,:,:] = allfile01.variables["CHEM-L_S__OH"][:,:,:]
OH[5,1,:,:,:] = allfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[5,2,:,:,:] = allfile02.variables["CHEM-L_S__OH"][:,:,:]
OH[5,3,:,:,:] = allfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[5,4,:,:,:] = allfile03.variables["CHEM-L_S__OH"][:,:,:]
OH[5,5,:,:,:] = allfile23.variables["CHEM-L_S__OH"][:,:,:]
OH[5,6,:,:,:] = allfile04.variables["CHEM-L_S__OH"][:,:,:]
OH[5,7,:,:,:] = allfile24.variables["CHEM-L_S__OH"][:,:,:]

pre_mt = 500
OH_avg = np.mean(OH,axis=1)
OH_mavg = np.sum(OH_avg[:,:,:,:]*AM,axis=1)/np.sum(AM,axis=0)

figure(1,figsize=(14,7.5))

ax = subplot(231)

title(r'OH CTRL: Jan-April 500 hPa (1e5molec/cm3)')
plot1 = OH_mavg[0,:,:]
plot_data, lon = addcyclic(plot1,lon)
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data/1e5,vmin=1,vmax=15)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(232)
title(r'$\Delta$OH: OMI HCHO only')
plot_data = (OH_mavg[1,:,:]-OH_mavg[0,:,:])*100/OH_mavg[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(233)
title(r'$\Delta$OH: IASI only')
plot_data = (OH_mavg[2,:,:] - plot1)*100/OH_mavg[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(234)

title(r'$\Delta$OH: MOPITT only')
plot2 = (OH_mavg[3,:,:]-OH_mavg[0,:,:])*100/OH_mavg[0,:,:]
plot_data= plot2
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(235)
title(r'$\Delta$OH: OMI NO2 only')
plot_data = (OH_mavg[4,:,:]-OH_mavg[0,:,:])*100/OH_mavg[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(236)
title(r'$\Delta$OH: All instrument')
plot_data = (OH_mavg[5,:,:]-OH_mavg[0,:,:])*100/OH_mavg[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
show()
