#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import os
def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1608_AM/ctm.00.20160801.nc')
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
len1,len2,len3 = shape(infile00.variables["IJ-AVG-S__Ox"][:,:,:])
len4 = 4
len5 = 24
O3 = np.zeros((len4,len5,len1,len2,len3))
O3_crt = np.zeros((len4,len5,5,len2,len3))
pre_crt = np.array([800,500,300])
ps = infile00.variables["DAO-FLDS__PS-PBL"][0,:,:]
for j in range(4):
    for i in range(12):
        month = str(i+1).zfill(2)
        runbuffer = '/runs/v8-02-01/'
        if j==0:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_M_V35/'
            rundir = 'geos5_ap_20'
            mondir1 = '16'
            mondir2 = '16'
        elif j==1:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_M_SOB/'
            rundir = 'geos5_iasio3_16'
            mondir1 = month+'_WC'
            mondir2 = str(i+21).zfill(2)+'_WC'
            sf_ind = '20'
        elif j==2:
            sysdir = '/users/jk/15/xzhang/'
            codedir = 'gcadj_std_T_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_07'
            mondir2 = month+'_27'
            sf_ind = '20'
        elif j==3:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_I_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_08'
            mondir2 = month+'_28'
            sf_ind = '20'
        path = sysdir+codedir+runbuffer+rundir
        if j>0:
            path_sf = path+mondir1+'/ctm.'+sf_ind+'.2016'+month+'01.nc'
            path_sf2 = path+mondir2+'/ctm.'+sf_ind+'.2016'+month+'16.nc'
            if os.path.isfile(path_sf) == False:
                sf_ind = '15'
                path_sf =  path+mondir1+'/ctm.'+sf_ind+'.2016'+month+'01.nc'
            if os.path.isfile(path_sf2) == False:
                sf_ind = '15'
                path_sf2 =  path+mondir2+'/ctm.'+sf_ind+'.2016'+month+'16.nc'
            spec_in = Dataset(path_sf)
            spec_in2 = Dataset(path_sf2)
            O3[j,2*i,:,:,:] = spec_in.variables["IJ-AVG-S__Ox"][:,:,:]
            O3[j,2*i+1,:,:,:] = spec_in2.variables["IJ-AVG-S__Ox"][:,:,:]
        else:
            path_ap = path+mondir1+'/ctm.00.2016'+month+'01.nc'
            path_ap2 = path+mondir2+'/ctm.00.2016'+month+'16.nc'
            spec_ap = Dataset(path_ap)
            spec_ap2 = Dataset(path_ap2)
            O3[j,2*i,:,:,:] = spec_ap.variables["IJ-AVG-S__Ox"][:,:,:]
            O3[j,2*i+1,:,:,:] = spec_ap2.variables["IJ-AVG-S__Ox"][:,:,:]
        for m in range(3):
            O3_crt[j,2*i,m,:,:] = geos_interpolate(len1,lat,lon,O3[j,2*i,:,:,:],ps,pre_crt[m])
            O3_crt[j,2*i+1,m,:,:] = geos_interpolate(len1,lat,lon,O3[j,2*i+1,:,:,:],ps,pre_crt[m])
            #O3_crt_mon[j,i,m,:,:] = np.mean(O3_crt[j,2*i:2*i+2,m,:,:],axis=0)
O3_crt_sea = np.zeros((len4,4,5,len2,len3))
O3_crt_sea[:,0,:,:,:] = (2*np.mean(O3_crt[:,0:4,:,:,:],axis=1)+np.mean(O3_crt[:,22:24,:,:,:],axis=1))/3
O3_crt_sea[:,1,:,:,:] = np.mean(O3_crt[:,4:10,:,:,:],axis=1)
O3_crt_sea[:,2,:,:,:] = np.mean(O3_crt[:,10:16,:,:,:],axis=1)
O3_crt_sea[:,3,:,:,:] = np.mean(O3_crt[:,16:22,:,:,:],axis=1)
O3_crt_ts = np.zeros((len4,24,5,3))
O3_crt_ts[:,:,:,0] = np.mean(np.mean(O3_crt[:,:,:,4:15,:],axis=3),axis=3)
O3_crt_ts[:,:,:,1] = np.mean(np.mean(O3_crt[:,:,:,15:31,:],axis=3),axis=3)
O3_crt_ts[:,:,:,2] = np.mean(np.mean(O3_crt[:,:,:,31:42,:],axis=3),axis=3)
title2 = ['DJF ','MAM ','JJA ','SON ']
title3 = ['CTRL','IASI','All instrument']
for j in range(3):
    for m in range(3):
        figure(3*j+m+1,figsize=(12,7.5))
        for k in range(4):
            ax = subplot(221+k)
            m1=Basemap(lon_0=0)
            m1.drawcoastlines()
            m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
            m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
#For CO
            if j==0:
                title1 = r'$O_3$ concentration: '
                plot_data, lon = addcyclic(O3_crt_sea[j,k,m,:,:],lon)
                if k==0 and m==0:
                    x,y = m1(*np.meshgrid(lon,lat))
                if m==2:
                    anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=130)
                else:
                    anth_plot = m1.pcolormesh(x,y,plot_data,vmin=1,vmax=80)
            else:
                title1 = r'$\Delta$$O_3$ concentration: '
                plot_data = O3_crt_sea[j,k,m,:,:]-O3_crt_sea[0,k,m,:,:]
                anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-15,vmax=15,cmap=get_cmap("RdBu_r"))
            if k==0:
                title(title1+title2[k]+str(pre_crt[m])+'hPa '+title3[j])
            else:
                title(title1+title2[k])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="7%", pad=0.3)
            plt.colorbar(orientation = 'horizontal',cax=cax)
            anth_plot.cmap.set_bad('lightgrey')
figure(10,figsize=(16,12))
x_series = arange(24)
title1=['800hPa','500hPa','300hPa']
title2=['75-30S','30S-30N','30-75N']
for i in range(3):
    for j in range(3):
        ax = subplot(331+3*i+j)
        plot(x_series,O3_crt_ts[0,:,i,j],x_series,O3_crt_ts[1,:,i,j],x_series,O3_crt_ts[2,:,i,j],linewidth=2.0)
        xticks([0,2,4,6,8,10,12,14,16,18,20,22],['J','F','M','A','M','J','J','A','S','O','N','D'])
        title(r'O3 distribution'+title1[i]+title2[j])
        if 3*i*j==0:
            legend(('CTRL','MOPITT','All instrument'),loc='best')
    #ylim(10,14)
show()
