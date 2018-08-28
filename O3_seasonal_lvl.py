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
iasifile01 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1601/ctm.15.20160101.nc')
iasifile21 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1601_2/ctm.15.20160116.nc')
iasifile02 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1602/ctm.15.20160201.nc')
iasifile22 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1602_2/ctm.15.20160216.nc')
iasifile03 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1603/ctm.15.20160301.nc')
iasifile23 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1603_2/ctm.15.20160316.nc')
iasifile04 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1604/ctm.15.20160401.nc')
iasifile24 = Dataset('/users/jk/16/xzhang/gcadj_std_M_SOB/runs/v8-02-01/geos5_iasio3_1604_2/ctm.15.20160416.nc')
allfile01 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1601_05/ctm.20.20160101.nc')
allfile21 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1601_25/ctm.15.20160116.nc')
allfile02 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1602_05/ctm.15.20160201.nc')
allfile22 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1602_25/ctm.20.20160216.nc')
allfile03 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1603_05/ctm.15.20160301.nc')
allfile23 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1603_25/ctm.20.20160316.nc')
allfile04 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1604_05/ctm.20.20160401.nc')
allfile24 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1604_25/ctm.20.20160416.nc')


#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][0,:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len3 = shape(iasifile01.variables["IJ-AVG-S__Ox"][:,:,:])
len4 = 3
len5 = 8

O3 = np.zeros((len4,len5,len1,len2,len3))

O3[0,0,:,:,:] = infile01.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,1,:,:,:] = infile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,2,:,:,:] = infile02.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,3,:,:,:] = infile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,4,:,:,:] = infile03.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,5,:,:,:] = infile23.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,6,:,:,:] = infile04.variables["IJ-AVG-S__Ox"][:,:,:]
O3[0,7,:,:,:] = infile24.variables["IJ-AVG-S__Ox"][:,:,:]

O3[1,0,:,:,:] = iasifile01.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,1,:,:,:] = iasifile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,2,:,:,:] = iasifile02.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,3,:,:,:] = iasifile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,4,:,:,:] = iasifile03.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,5,:,:,:] = iasifile23.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,6,:,:,:] = iasifile04.variables["IJ-AVG-S__Ox"][:,:,:]
O3[1,7,:,:,:] = iasifile24.variables["IJ-AVG-S__Ox"][:,:,:]

O3[2,0,:,:,:] = allfile01.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,1,:,:,:] = allfile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,2,:,:,:] = allfile02.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,3,:,:,:] = allfile21.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,4,:,:,:] = allfile03.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,5,:,:,:] = allfile23.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,6,:,:,:] = allfile04.variables["IJ-AVG-S__Ox"][:,:,:]
O3[2,7,:,:,:] = allfile24.variables["IJ-AVG-S__Ox"][:,:,:]

pre_mt = 700
O3_avg = np.mean(O3,axis=1)
O3_mt = np.zeros((len4,len2,len3))
for i in range(len4):
    O3_mt[i,:,:] = geos_interpolate(len1,lat,lon,O3_avg[i,:,:,:],ps,pre_mt)
pre_utls = 900
O3_utls = np.zeros((len4,len2,len3))
for i in range(len4):
    O3_utls[i,:,:] = geos_interpolate(len1,lat,lon,O3_avg[i,:,:,:],ps,pre_utls)

figure(1,figsize=(14,7.5))

ax = subplot(231)

title(r'$O_3$ concentrations CTRL: Jan-April 700 hPa[ppbv]')

plot_data, lon = addcyclic(O3_mt[0,:,:],lon)
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=80)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(232)
title(r'$\Delta$$O_3$: IASI only')
plot_data = O3_mt[1,:,:]-O3_mt[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(233)
title(r'$\Delta$$O_3$: All instrument')
plot_data = O3_mt[2,:,:]-O3_mt[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(234)

title(r'$O_3$ concentrations CTRL: Jan-April 900 hPa[ppbv]')

plot_data= O3_utls[0,:,:]
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=120)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(235)
title(r'$\Delta$$O_3$: IASI only')
plot_data = O3_utls[1,:,:]-O3_utls[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(236)
title(r'$\Delta$$O_3$: All instrument')
plot_data = O3_utls[2,:,:]-O3_utls[0,:,:]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For O3 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
show()
