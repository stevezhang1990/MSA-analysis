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
infile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/ctm.01.20160101.nc')
infile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/ctm.01.20160116.nc')
infile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/ctm.01.20160201.nc')
infile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/ctm.01.20160216.nc')
infile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/ctm.01.20160301.nc')
infile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/ctm.01.20160316.nc')
infile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/ctm.01.20160401.nc')
infile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/ctm.01.20160416.nc')
omhfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/gctm.sf.10.20160101.nc')
omhfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/gctm.sf.10.20160116.nc')
omhfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/gctm.sf.10.20160201.nc')
omhfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/gctm.sf.10.20160216.nc')
omhfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/gctm.sf.10.20160301.nc')
omhfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/gctm.sf.10.20160316.nc')
omhfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/gctm.sf.10.20160401.nc')
omhfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/gctm.sf.10.20160416.nc')
allfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_06/gctm.sf.20.20160101.nc')
allfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_26/gctm.sf.15.20160116.nc')
allfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_06/gctm.sf.20.20160201.nc')
allfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_26/gctm.sf.20.20160216.nc')
allfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_06/gctm.sf.15.20160301.nc')
allfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_26/gctm.sf.20.20160316.nc')
allfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_06/gctm.sf.15.20160401.nc')
allfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_26/gctm.sf.15.20160416.nc')


#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2 = shape(omhfile01.variables["IJ-EMS-S__ISOP_an"][:,:])
ISOP_an = np.zeros((8,len1,len2))
NOX_bb = np.zeros((8,len1,len2))

SF_omh_an = np.zeros((8,len1,len2))
SF_omh_bb = np.zeros((8,len1,len2))

SF_all_an = np.zeros((8,len1,len2))
SF_all_bb = np.zeros((8,len1,len2))

ISOP_an[0,:,:] = infile01.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[1,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[2,:,:] = infile02.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[3,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[4,:,:] = infile03.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[5,:,:] = infile23.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[6,:,:] = infile04.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[7,:,:] = infile24.variables["BIOGSRCE__ISOP"][0,:,:]

SF_omh_an[0,:,:] = omhfile01.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[1,:,:] = omhfile21.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[2,:,:] = omhfile02.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[3,:,:] = omhfile21.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[4,:,:] = omhfile03.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[5,:,:] = omhfile23.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[6,:,:] = omhfile04.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[7,:,:] = omhfile24.variables["IJ-EMS-S__ISOP_an"][:,:]


SF_all_an[0,:,:] = allfile01.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[1,:,:] = allfile21.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[2,:,:] = allfile02.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[3,:,:] = allfile21.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[4,:,:] = allfile03.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[5,:,:] = allfile23.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[6,:,:] = allfile04.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[7,:,:] = allfile24.variables["IJ-EMS-S__ISOP_an"][0,:,:]

ISOP_ap = ISOP_an
ISOP_omh = SF_omh_an
ISOP_all = SF_all_an
#ISOP_tes = ISOP_tes_an+ISOP_tes_bb -1

figure(1,figsize=(14,7.5))

ax = subplot(231)

title('ISOP emissions: Jan-April [1e10 molec/cm3/s]')

plot_data, lon = addcyclic(np.mean(ISOP_ap[:,:,:],axis=0),lon)
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For ISOP 
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=200,cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(232)
title('ISOP emissions SF: OMI CH2O only')
plot_data = np.mean(ISOP_omh[:,:,:],axis=0)
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For ISOP 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0.5,vmax=1.5,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(233)
title('ISOP emissions SF: All instrument')
plot_data = np.mean(ISOP_all[:,:,:],axis=0)
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For ISOP 
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0.5,vmax=1.5,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

show()
