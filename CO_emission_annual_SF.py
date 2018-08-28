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
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon = infile00.variables["LON"][:]
len2 = len(lon)
lat = infile00.variables["LAT"][:]
ps = infile00.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile00.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1 = len(lat)
len3 = 5
len4 = 24
CO_ab = np.zeros((len3,4,len4,len1,len2))
title2 = ['DJF','MAM','JJA','SON']
title3 = ['CTRL','MOPITT','MSA','MSA']
for j in range(5):
    for i in range(12):
        month = str(i+1).zfill(2)
        runbuffer = '/runs/v8-02-01/'
        if j==0:
            title1 = 'CTRL'
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_M_V35/'
            rundir = 'geos5_ap_20'
            mondir1 = '16'
            mondir2 = '16'
        elif j==1:
            title1 = 'MOPITT'
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_M_V35/'
            rundir = 'geos5_mop_16'
            mondir1 = month
            mondir2 = str(i+21).zfill(2)     
            sf_ind = '10'
        elif j==2:
            title1 = 'All instrument'
            sysdir = '/users/jk/15/xzhang/'
            codedir = 'gcadj_std_T_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_07'
            mondir2 = month+'_27'
            sf_ind = '20'
        elif j==3:
            title1 = 'All instrument'
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_I_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_08'
            mondir2 = month+'_28'
            sf_ind = '20'
        path = sysdir+codedir+runbuffer+rundir
        if j==0:
            path_ap = path+mondir1+'/ctm.00.2016'+month+'01.nc'
            path_ap2 = path+mondir2+'/ctm.00.2016'+month+'16.nc'
        elif j<4:
            path_ap = path+mondir1+'/ctm.01.2016'+month+'01.nc'
            path_ap2 = path+mondir2+'/ctm.01.2016'+month+'16.nc'
        spec_ap = Dataset(path_ap)
        spec_ap2 = Dataset(path_ap2)
        if j>0 and j<4:
            path_sf = path+mondir1+'/gctm.sf.'+sf_ind+'.2016'+month+'01.nc'
            path_sf2 = path+mondir2+'/gctm.sf.'+sf_ind+'.2016'+month+'16.nc'
            if os.path.isfile(path_sf) == False:
                sf_ind = '15'
                path_sf =  path+mondir1+'/gctm.sf.'+sf_ind+'.2016'+month+'01.nc'
            if os.path.isfile(path_sf2) == False:
                sf_ind = '15'
                path_sf2 =  path+mondir2+'/gctm.sf.'+sf_ind+'.2016'+month+'16.nc'
            spec_in = Dataset(path_sf)
            spec_in2 = Dataset(path_sf2)
#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
        CO_ab[j,0,2*i,:,:] = spec_ap.variables["CO--SRCE__COanth"][0,:,:]
        CO_ab[j,0,2*i+1,:,:] = spec_ap2.variables["CO--SRCE__COanth"][0,:,:]
        CO_ab[j,1,2*i,:,:] = spec_ap.variables["CO--SRCE__CObb"][0,:,:]
        CO_ab[j,1,2*i+1,:,:] = spec_ap2.variables["CO--SRCE__CObb"][0,:,:]
        if j==1:
            CO_ab[j,2,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_an"][:,:]
            CO_ab[j,2,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_an"][:,:]
            CO_ab[j,3,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_bb"][:,:]
            CO_ab[j,3,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_bb"][:,:]
        elif j>1 and j<4:
            CO_ab[j,2,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_an"][0,:,:]
            CO_ab[j,2,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_an"][0,:,:]
            CO_ab[j,3,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_bb"][0,:,:]
            CO_ab[j,3,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_bb"][0,:,:]
    CO_ap = CO_ab[j,0,:,:,:]+CO_ab[j,1,:,:,:]
    CO_sf = CO_ab[j,2,:,:,:]+CO_ab[j,3,:,:,:] -1
    if j==2:
        CO_sf = 1.1*CO_ab[j,2,:,:,:]+1.1*CO_ab[j,3,:,:,:] -1.2
#CO_mopitt = CO_mopitt_an+CO_mopitt_bb-1
#CO_omi = CO_omi_an+CO_omi_bb -1
#CO_tes = CO_tes_an+CO_tes_bb -1
    figure(j,figsize=(10,7.5))

    if j==0:
        CO_plot = CO_ap
    elif j<4:
        CO_plot = CO_sf
    for k in range(4):
        ax = subplot(221+k)
        m1=Basemap(lon_0=0)
        m1.drawcoastlines()
        m1.drawparallels(arange(-90,120,30),labels=[1,0,0,1],labelstyle="+/-")
        m1.drawmeridians(arange(-180,210,60),labels=[1,0,0,1],labelstyle="+/-")
        if k==0 and j<4:
            title('CO emissions: '+title2[k]+' '+title1)
            if j==0:
                plot_data, lon = addcyclic((2*np.mean(CO_plot[0:4,:,:],axis=0)+\
                                        np.mean(CO_plot[22:24,:,:],axis=0))/3,lon)
                x,y = m1(*np.meshgrid(lon,lat))
            else:
                plot_data = (2*np.mean(CO_plot[0:4,:,:],axis=0)+\
                                        np.mean(CO_plot[22:24,:,:],axis=0))/3
        elif k>0 and j<4:
            title('CO emissions: '+title2[k])
            plot_data = np.mean(CO_plot[6*k-2:6*k+4,:,:],axis=0)
        elif k==0 and j==4:
            title('CO annual emissions: '+title3[k])
            plot_data = np.mean(CO_ab[k,0,:,:,:]+CO_ab[k,1,:,:,:],axis=0)
            #x,y = m1(*np.meshgrid(lon,lat))
        elif k>0 and j==4:
            title('CO annual SF: '+title3[k])
            if k!=2:
                plot_data = np.mean(CO_ab[k,2,:,:,:]+CO_ab[k,3,:,:,:]-1,axis=0)
            else:
                plot_data = np.mean(1.1*CO_ab[k,2,:,:,:]+1.1*CO_ab[k,3,:,:,:]-1.2,axis=0)
    #For CO
        if j==0 or (k==0 and j==4):
            anth_plot = m1.pcolormesh(x,y,plot_data/1e10,vmin=0,vmax=200,cmap=get_cmap("gnuplot2_r"))
        else:
            anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0.5,vmax=1.5,cmap=get_cmap("RdBu_r"))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="7%", pad=0.3)
        plt.colorbar(orientation = 'horizontal',cax=cax)
        anth_plot.cmap.set_bad('lightgrey')

show()
