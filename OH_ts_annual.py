#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
# month cycle indicator
fig = plt.figure(1,figsize=(12,7.5))
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
# sample dataset reading in lat lon and dimension info
infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1608_AM/ctm.00.20160801.nc')
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
len1,len2,len5 = shape(infile00.variables["CHEM-L_S__OH"][:,:,:])
len3 = 7
len4 = 24
OH = np.zeros((len3,len4,len1,len2,len5))
AM = infile00.variables["BXHGHT-S__AD"][0:38,:,:]
# Global mass weighted mean OH from Spivakovsky et al 2000, and TransCOM
spi = np.array([11.265,10.91,10.99,10.7,10.73,11.92,12.44,12.77,12.15,11.68,11.31,11.28])
transcom = np.array([9.15,9.02,9.11,8.95,8.96,9.62,10.29,10.64,10.15,11.68,9.21,9.18])
# Lists of strings for different file entries
sysdir = ['/users/jk/16/xzhang/','/users/jk/15/xzhang/']
codedir = ['gcadj_std_M_V35/','gcadj_std_M_V35/','gcadj_std_O_3d','gcadj_std_O_V35','gcadj_std_M_SOB','gcadj_std_T_3d','gcadj_std_I_3d']
runbuffer = '/runs/v8-02-01/'
rundir = ['geos5_ap_2016','geos5_mop_16','geos5_omino2_16','geos5_omh_16','geos5_iasio3_16','geos5_all_16','geos5_all_16']
rundir2 = ['geos5_ap_2017','geos5_mop_17','geos5_omino2_17','geos5_omh_17','geos5_iasio3_17','geos5_all_17','geos5_all_17']
sf_ind = ['10','10','10','20','20','20']
# Looping among different experiments
for j in range(7):
    for i in range(14):
        month = str(i+1).zfill(2)
        mondir1 = [month,month,month,month+'WC',month+'_07',month+'_08']
        mondir2 = [str(i+21).zfill(2),str(i+21).zfill(2),str(i+21).zfill(2),str(i+21).zfill(2)+'WC',month+'_27',month+'_28']        
        if i<12:
            rund = rundir[j]
            yeard = '.2016'
        else:
            rund = rundir2[j]
            yeard = '.2017'
        if j!=5:
            sysd = sysdir[0]
        else:
            sysd = sysdir[1]
        path = sysd+codedir[j]+runbuffer+rund
        if j>0:
            path_sf = path+mondir1[j-1]+'/ctm.'+sf_ind+yeard+month+'01.nc'
            path_sf2 = path+mondir2[j-1]+'/ctm.'+sf_ind+yeard+month+'16.nc'
            if os.path.isfile(path_sf) == False:
                sf_ind = '15'
                path_sf =  path+mondir1[j-1]+'/ctm.'+sf_ind+yeard+month+'01.nc'
            if os.path.isfile(path_sf2) == False:
                sf_ind = '15'
                path_sf2 =  path+mondir2[j-1]+'/ctm.'+sf_ind+yeard+month+'16.nc'
            spec_in = Dataset(path_sf)
            spec_in2 = Dataset(path_sf2)
            OH[j,2*i,:,:,:] = spec_in.variables["CHEM-L_S__OH"][:,:,:]
            OH[j,2*i+1,:,:,:] = spec_in2.variables["CHEM-L_S__OH"][:,:,:]
        else:
            path_ap = path++'/ctm.00'+yeard+month+'01.nc'
            path_ap2 = path+'/ctm.00'+yeard+month+'16.nc'
            spec_ap = Dataset(path_ap)
            spec_ap2 = Dataset(path_ap2)
            OH[j,2*i,:,:,:] = spec_ap.variables["CHEM-L_S__OH"][:,:,:]
            OH[j,2*i+1,:,:,:] = spec_ap2.variables["CHEM-L_S__OH"][:,:,:]
OH[6,8:12,:,:,:] = 0.98*OH[5,8:12,:,:,:]
#Project 1: global mean OH time series        
len3p = 9
OH_avg = np.zeros((len3p,len4))
# Compute monthly mean OH
for j in range(len4):
    for i in range(len3p):
        if i==7:
            OH_avg[7,j] = spi[int(j/2)]
        elif i==8:
            OH_avg[8,j] = transcom[int(j/2)]
        else:
            OH_avg[i,j] = np.sum(OH[i,j,:,:,:]*AM)/(np.sum(AM)*1e5)
#OH_avg[6,8:12] = OH_avg[5,8:12]        
fig = plt.figure(1,figsize=(12,4.5))
ax=subplot(211)
x_series = arange(24)
plot(x_series,OH_avg[0,:],x_series,OH_avg[1,:],x_series,OH_avg[6,:],x_series,OH_avg[7,:],linewidth=2.0)
xticks([0,2,4,6,8,10,12,14,16,18,20,22],['J','F','M','A','M','J','J','A','S','O','N','D'])
title("Global mean OH concentrations (1e5 molec/cm3) 2016")
legend(('CTRL','MOPITT','All instrument','Spivakovsky et al 2000'),loc='best')
ylim(10,14)

#Project 2: latitudinal distribution of OH
OH_mod = np.mean(OH[:,:,:,:,:],axis=1)
OH_lat_avg = np.zeros((len3,len2))
OH_avg_2d = np.zeros((len3,len2,len5))
for i in range(len3):
    OH_avg_2d[i,:,:] = np.sum(OH_mod[i,:,:,:]*AM,axis=0)/np.sum(AM,axis=0)
    OH_lat_avg[i,:] =np.sum(np.sum(OH_mod[i,:,:,:]*AM,axis=0),axis=1)/np.sum(np.sum(AM,axis=0),axis=1)
ax=subplot(212)
plot(lat,0*OH_lat_avg[0,:], 'k--',
     lat,(OH_lat_avg[1,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],
     lat,(OH_lat_avg[2,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],
     lat,(OH_lat_avg[3,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],
     lat,(OH_lat_avg[4,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],
     lat,(OH_lat_avg[6,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],'r--',linewidth=2.0)
     #lat,(OH_lat_avg[6,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],'r--',linewidth=2.0)
     #lat,(OH_lat_avg[7,:]-OH_lat_avg[0,:])*100/OH_lat_avg[0,:],'r--',linewidth=2.0)
     #lat,(OH_avg[9,:]-OH_avg[0,:])*50/OH_avg[0,:],'r--',linewidth=2.0)
legend(('CTRL','MOPITT','OMI NO2','OMI HCHO','IASI O3','MSA'),loc='best')#'MSA 08','Spivakovsky et al 2000'
#legend(('MOPITT','IASI+OSIRIS','All instrument'),loc='best')
title('Zonal mass weighted mean OH relative to CTRL run in %')
xlim(-90,90)
ylim(-20,10)
fig = plt.figure(2,figsize=(15,7.5))
title2 = ['MOPITT','OMI NO2','OMI HCHO','IASI O3','MSA']
for i in range(6):
    ax = subplot(231+i)
    m1=Basemap(lon_0=0)
    m1.drawcoastlines()
    m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-",fontsize=7)
    m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-",fontsize=7)
    if i==0:
        title(r'Global mean OH (1e5 molec/cm3) CTRL',fontsize=11)
        plot_data,lon=addcyclic(OH_avg_2d[0,:,:],lon)
        x,y = m1(*np.meshgrid(lon,lat))
        anth_plot = m1.pcolormesh(x,y,plot_data/1e5,vmin=0,vmax=30,cmap=get_cmap("gnuplot2_r"))
    else:
        title(r'DOH in % for '+title2[i-1],fontsize=11)
        plot_data = 100*(OH_avg_2d[i,:,:]-OH_avg_2d[0,:,:])/OH_avg_2d[0,:,:]
        anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("bwr"))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="7%", pad=0.3)
    cbar = plt.colorbar(orientation = 'horizontal',cax=cax)
    cbar.ax.tick_params(labelsize=8)
    anth_plot.cmap.set_bad('lightgrey')
    plt.subplots_adjust(hspace=0,wspace=0.1)
show()
