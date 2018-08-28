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
infile05 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1605/ctm.01.20160501.nc')
infile25 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1605_2/ctm.01.20160516.nc')
infile06 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1606/ctm.01.20160601.nc')
infile26 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1606_2/ctm.01.20160616.nc')
infile07 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1607/ctm.01.20160701.nc')
infile27 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1607_2/ctm.01.20160716.nc')
infile08 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1608/ctm.01.20160801.nc')
infile28 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1608_2/ctm.01.20160816.nc')
infile09 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1609/ctm.01.20160901.nc')
infile29 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1609_2/ctm.01.20160916.nc')
infile10 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1610/ctm.01.20161001.nc')
infile210= Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1610_2/ctm.01.20161016.nc')
infile11 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1611/ctm.01.20161101.nc')
infile211 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1611_2/ctm.01.20161116.nc')
infile12 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1612/ctm.01.20161201.nc')
infile212 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1612_2/ctm.01.20161216.nc')


#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile01.variables["DXYP__DXYP"][0,:,:] * (100**2) # cm^2
len1,len2 = shape(infile01.variables["BIOGSRCE__ISOP"][0,:,:])
ISOP_an = np.zeros((24,len1,len2))
ISOP_bb = np.zeros((24,len1,len2))
ISOP_an[0,:,:] = infile01.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[1,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[2,:,:] = infile02.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[3,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[4,:,:] = infile03.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[5,:,:] = infile23.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[6,:,:] = infile04.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[7,:,:] = infile24.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[8,:,:] = infile05.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[9,:,:] = infile25.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[10,:,:] = infile06.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[11,:,:] = infile26.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[12,:,:] = infile07.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[13,:,:] = infile27.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[14,:,:] = infile08.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[15,:,:] = infile28.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[16,:,:] = infile09.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[17,:,:] = infile29.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[18,:,:] = infile10.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[19,:,:] = infile210.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[20,:,:] = infile11.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[21,:,:] = infile211.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[22,:,:] = infile12.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[23,:,:] = infile212.variables["BIOGSRCE__ISOP"][0,:,:]


ISOP_ap = ISOP_an
#ISOP_mopitt = ISOP_mopitt_an+ISOP_mopitt_bb-1
#ISOP_omi = ISOP_omi_an+ISOP_omi_bb -1
#ISOP_tes = ISOP_tes_an+ISOP_tes_bb -1

figure(1,figsize=(12,10))

ax = subplot(321)

title('ISOP emissions: annual emissions [1e9 molec/cm3/s]')

plot_data, lon = addcyclic(np.mean(ISOP_ap[:,:,:],axis=0),lon)
#print plot_data[15,15]
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")

x,y = m1(*np.meshgrid(lon,lat))

#For ISOP 50 for an 150 for anth
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=500, cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(322)
title('ISOP emissions: DJF')
plot_data = (2*np.mean(ISOP_ap[0:4,:,:],axis=0)+np.mean(ISOP_ap[22:23,:,:],axis=0))/3
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ISOP 50 for an 150 for anth
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=500, cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(323)
title('ISOP emissions: MAM')
plot_data = np.mean(ISOP_ap[4:10,:,:],axis=0)
m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ISOP 50 for an 150 for anth
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=500, cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')


ax = subplot(324)
title('ISOP emissions: JJA')

plot_data = np.mean(ISOP_ap[10:16,:,:],axis=0)

m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ISOP 50 for an 150 for anth
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=100, cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')

ax = subplot(325)
title('ISOP emissions: SON')

plot_data = np.mean(ISOP_ap[16:22,:,:],axis=0)

m1=Basemap(lon_0=0)
m1.drawcoastlines()
m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ISOP 50 for an 150 for anth
anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=1,vmax=500, cmap=get_cmap("hot_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')


molISOP = 6.022*10**23 / ( 68.12 )
#ISOP is reported to be either 46 or 40.03
s2d = 60*60*24
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
###World
lon_mask_world = (lon >= -180) & (lon <= 180)
lat_mask_world = (lat >= -90) & (lat<= 90)
ems_sum = np.zeros(24)
for i in range(24):
    ems_sum[i] = (ISOP_ap*surface)[i,lat_mask_world,:][:,lon_mask_world[:-1]].sum()/molISOP*s2d*d2y[i]/10**12
ax = subplot(326)
x_series = arange(24)
plot(x_series,ems_sum)
xticks([0,2,4,6,8,10,12,14,16,18,20,22],['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'])
title('ISOP emissions: time series')

show()
