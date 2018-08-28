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
omhfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/ctm.10.20160101.nc')
omhfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/ctm.10.20160116.nc')
omhfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/ctm.10.20160201.nc')
omhfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/ctm.10.20160216.nc')
omhfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/ctm.10.20160301.nc')
omhfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/ctm.10.20160316.nc')
omhfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/ctm.10.20160401.nc')
omhfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/ctm.10.20160416.nc')
allfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_06/ctm.20.20160101.nc')
allfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_26/ctm.15.20160116.nc')
allfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_06/ctm.20.20160201.nc')
allfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_26/ctm.20.20160216.nc')
allfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_06/ctm.15.20160301.nc')
allfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_26/ctm.20.20160316.nc')
allfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_06/ctm.15.20160401.nc')
allfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_26/ctm.15.20160416.nc')


#XOX-BIOB__XCH2O or anTHSRCE__XCH2O or XOX-an-S__XCH2O or xOX-LI-S__xCH2O
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][0,:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len3 = shape(omhfile01.variables["IJ-AVG-S__CH2O"][:,:,:])
len4 = 3
len5 = 8

CH2O = np.zeros((len4,len5,len1,len2,len3))

CH2O[0,0,:,:,:] = infile01.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,1,:,:,:] = infile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,2,:,:,:] = infile02.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,3,:,:,:] = infile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,4,:,:,:] = infile03.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,5,:,:,:] = infile23.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,6,:,:,:] = infile04.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[0,7,:,:,:] = infile24.variables["IJ-AVG-S__CH2O"][:,:,:]

CH2O[1,0,:,:,:] = omhfile01.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,1,:,:,:] = omhfile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,2,:,:,:] = omhfile02.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,3,:,:,:] = omhfile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,4,:,:,:] = omhfile03.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,5,:,:,:] = omhfile23.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,6,:,:,:] = omhfile04.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[1,7,:,:,:] = omhfile24.variables["IJ-AVG-S__CH2O"][:,:,:]

CH2O[2,0,:,:,:] = allfile01.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,1,:,:,:] = allfile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,2,:,:,:] = allfile02.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,3,:,:,:] = allfile21.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,4,:,:,:] = allfile03.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,5,:,:,:] = allfile23.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,6,:,:,:] = allfile04.variables["IJ-AVG-S__CH2O"][:,:,:]
CH2O[2,7,:,:,:] = allfile24.variables["IJ-AVG-S__CH2O"][:,:,:]

len6 = len(infile01.variables["CHEM-L_S__NO2"][:,0,0])
NO2_trop = np.zeros_like(CH2O)
NO2 = np.zeros((len4,len5,len6,len2,len3))

NO2[0,0,:,:,:] = infile01.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,1,:,:,:] = infile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,2,:,:,:] = infile02.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,3,:,:,:] = infile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,4,:,:,:] = infile03.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,5,:,:,:] = infile23.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,6,:,:,:] = infile04.variables["CHEM-L_S__NO2"][:,:,:]
NO2[0,7,:,:,:] = infile24.variables["CHEM-L_S__NO2"][:,:,:]

NO2[1,0,:,:,:] = omhfile01.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,1,:,:,:] = omhfile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,2,:,:,:] = omhfile02.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,3,:,:,:] = omhfile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,4,:,:,:] = omhfile03.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,5,:,:,:] = omhfile23.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,6,:,:,:] = omhfile04.variables["CHEM-L_S__NO2"][:,:,:]
NO2[1,7,:,:,:] = omhfile24.variables["CHEM-L_S__NO2"][:,:,:]

NO2[2,0,:,:,:] = allfile01.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,1,:,:,:] = allfile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,2,:,:,:] = allfile02.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,3,:,:,:] = allfile21.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,4,:,:,:] = allfile03.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,5,:,:,:] = allfile23.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,6,:,:,:] = allfile04.variables["CHEM-L_S__NO2"][:,:,:]
NO2[2,7,:,:,:] = allfile24.variables["CHEM-L_S__NO2"][:,:,:]

NO2_trop[:,:,0:38,:,:] = NO2
pre_mt = 500
CH2O_avg = np.mean(CH2O,axis=1)
CH2O_mt = np.zeros((len4,len2,len3))
NO2_avg = 1e9*np.mean(NO2_trop,axis=1)
NO2_mt = np.zeros((len4,len2,len3))
for i in range(len4):
    CH2O_mt[i,:,:] = geos_interpolate(len1,lat,lon,CH2O_avg[i,:,:,:],ps,pre_mt)
    NO2_mt[i,:,:] = geos_interpolate(len1,lat,lon,NO2_avg[i,:,:,:],ps,pre_mt)
pre_utls = 300
CH2O_utls = np.zeros((len4,len2,len3))
NO2_utls = np.zeros((len4,len2,len3))
for i in range(len4):
    CH2O_utls[i,:,:] = geos_interpolate(len1,lat,lon,CH2O_avg[i,:,:,:],ps,pre_utls)
    NO2_utls[i,:,:] = geos_interpolate(len1,lat,lon,NO2_avg[i,:,:,:],ps,pre_utls)

figure(1,figsize=(14,7.5))

ax = subplot(231)

title(r'CH2O/NO2 CTRL: Jan-April 500 hPa')
plot1 = CH2O_mt[0,:,:]/NO2_mt[0,:,:]
plot_data, lon = addcyclic(plot1,lon)
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=50)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(232)
title(r'$\Delta$CH2O/NO2: IASI only')
plot_data = CH2O_mt[1,:,:]/NO2_mt[1,:,:] - plot1
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
title(r'$\Delta$CH2O/NO2: All instrument')
plot_data = CH2O_mt[2,:,:]/NO2_mt[2,:,:] - plot1
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

title(r'CH2O/NO2 CTRL: Jan-April 300 hPa')
plot2 = CH2O_utls[0,:,:]/NO2_utls[0,:,:]
plot_data= plot2
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For CH2O 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=50)
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
bond = 15
ax = subplot(235)
title(r'$\Delta$CH2O/NO2: IASI only')
plot_data = CH2O_utls[1,:,:]/NO2_utls[1,:,:]-plot2
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
title(r'$\Delta$CH2O/NO2: All instrument')
plot_data = CH2O_utls[2,:,:]/CH2O_utls[2,:,:]-plot2
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
