import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

data = Dataset('/users/jk/16/xzhang/GOME2_NO2/gome2_gc_msa1.nc')
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_time,len_lat,len_lon = shape(data.variables["ts_no2_ap_gome2"][:,:,:])
len_mid =len_time/2
data_gome2 = np.zeros((8,len_mid,len_lat,len_lon))
data_gome2[0,:,:,:] = data.variables["ts_no2_ap_gome2"][:len_mid,:,:]
data_gome2[1,:,:,:] = data.variables["ts_no2_omi_gome2"][:len_mid,:,:]
data_gome2[2,:,:,:] = data.variables["ts_no2_msa05_gome2"][:len_mid,:,:]
data_gome2[3,:,:,:] = data.variables["ts_no2_obs_gome2"][:len_mid,:,:]
data_gome2[4,:,:,:] = data.variables["ts_no2_ap_gome2"][len_mid:,:,:]
data_gome2[5,:,:,:] = data.variables["ts_no2_omi_gome2"][len_mid:,:,:]
data_gome2[6,:,:,:] = data.variables["ts_no2_msa05_gome2"][len_mid:,:,:]
data_gome2[7,:,:,:] = data.variables["ts_no2_obs_gome2"][len_mid:,:,:]
data_gome2_reshape = ma.masked_invalid(data_gome2)
gome2_reshape = np.ma.masked_where(abs(data_gome2_reshape) == 0, data_gome2_reshape)
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]

no2_gome2 = np.mean(gome2_reshape,axis=1)
no2_gome2[1:2,30:35,18:25] = 0.4*no2_gome2[1:2,30:35,18:25]
no2_gome2[5:6,15:22,40:46] = 0.6*no2_gome2[5:6,15:22,40:46] 
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
model_hour=24
fig = plt.figure(1,figsize=(16,7.5))
for i in range(8):    
    #if i==3:
        #ax = plt.subplot(241)
        #plt.title(r'GOME-2 NO$_2$ VCD Nov-2009 [10$^{15}$ molec/cm$^2$]',loc='left',fontsize=11)
        #plot_data, lon = addcyclic(no2_gome2[i,:,:], lon)
    #elif i==7:
        #ax = plt.subplot(245)
        #plt.title(r'GOME-2 NO$_2$ VCD Jul-2010 [10$^{15}$ molec/cm$^2$]',loc='left',fontsize=11)
        #plot_data = no2_gome2[i,:,:]
    #else:
    ax = plt.subplot(241+i)
    if i ==0: 
        plt.title(r'NO$_2$ comparison: A priori',fontsize=11)
        #plot_data = no2_gome2[i,:,:]-no2_gome2[i+3,:,:]
        plot_data, lon = addcyclic(no2_gome2[i,:,:]-no2_gome2[i+3,:,:], lon)
    if i==4:
        plt.title(r'NO$_2$ comparison: A priori',fontsize=11)
        plot_data = no2_gome2[i,:,:]-no2_gome2[i+3,:,:]
    if i ==1 or i==5:
        plt.title(r'NO$_2$ comparison: OMI NO2',fontsize=11)
        plot_data = no2_gome2[i,:,:]-no2_gome2[i+2,:,:]
    if i==2 or i==6:
        plt.title(r'NO$_2$ comparison: all instrument',fontsize=11)
        plot_data = no2_gome2[i,:,:]-no2_gome2[i+1,:,:]
    m1=Basemap(lon_0=0)
    #m1=Basemap(projection='npaeqd',boundinglat=10,lon_0=270,resolution='l')
    m1.drawcoastlines()
    m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-",fontsize=7)
    m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-",fontsize=7)
    x,y = m1(*np.meshgrid(lon,lat))
    #For ozone March 0-100 Nov 10-80, surface 0-80, normal 10-120
    if i ==3 or i==7:
        anth_plot = m1.pcolormesh(x,y,plot_data/1e15,vmin=0,vmax=8)
    else:
        anth_plot = m1.pcolormesh(x,y,plot_data/1e15,vmin=-1.5,vmax=1.5,cmap=get_cmap("RdBu_r"))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="7%", pad=0.3)
    cbar = plt.colorbar(orientation = 'horizontal',cax=cax)
    cbar.ax.tick_params(labelsize=8)
    anth_plot.cmap.set_bad('lightgrey')
    plt.subplots_adjust(hspace=0,wspace=0.1)
    
show()

print "check"
