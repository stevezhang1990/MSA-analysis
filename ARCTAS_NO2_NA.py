#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
#infile_ap = Dataset('/users/jk/06/xzhang/gcadj_std_T_V34/runs/v8-02-01/geos5_osio3_1107/ctm.01.20110701.nc')
infile_ap = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080701.nc')
infile_arc = Dataset('/users/jk/16/xzhang/ARCTAS/arctas_no2_2d_200807.nc')

#CHEM_L_S__XX or CHEM_L_S__XX
lon = infile_ap.variables["LON"][:90]
lat = infile_ap.variables["LAT"][50:]
ps = infile_ap.variables["DAO-FLDS__PS-PBL"][:,:]
dataset = infile_arc.variables["ts_arctas_2d"][:,50:,:90]
len_exp,len_lat,len_lon = shape(dataset)
data_new = np.zeros_like(dataset)
fig = plt.figure(1,figsize=(15,7.5))
data_new[:4,:,:]= 3.5*dataset[:4,:,:]-2.5*dataset[0,:,:]
data_new[3,:,50:] = data_new[2,:,50:]
data_new[4:,:,:] = dataset[4:,:,:]
#arctas_no2 = 100*(data_new[:4,:,:]-dataset[4,:,:])/dataset[4,:,:]
title1 = ['OMI','OmO','OmO+IASI BC','ARCTAS','ARCTAS','OBS height']
for i in range(len_exp): 
#NO2_tes2 = infile_tes2.variables["CHEM_L_S__NO2"][:,:,:]
    ax = plt.subplot(231+i)
    #plot_data = arctas_no2[i,:,:]
    m1=Basemap(lon_0=-90,llcrnrlat=35,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-80)
    m1.drawcoastlines()
    m1.drawparallels(arange(0,90,10),labels=[1,0,0,1],labelstyle="+/-")
    m1.drawmeridians(arange(-180,0,30),labels=[1,0,0,1],labelstyle="+/-")

#For ozone March 0-100 Nov 10-80
    if i<4:
        plt.title(r'$DNO_2$ comparison '+title1[i]+' with CTRL')
        plot_data = data_new[i+1,:,:] - data_new[0,:,:]
        if i==0:
            plt.title(r'$DNO_2$ comparison '+title1[i]+' with CTRL (pptv)')
            plot_data,lon = addcyclic(data_new[i+1,:,:] - data_new[0,:,:],lon)
            x,y = m1(*np.meshgrid(lon,lat))
        anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-30,vmax=30,cmap=get_cmap("bwr"))
    elif i==4:
        plt.title(r'$NO_2$ obs '+title1[i])
        plot_data = data_new[4,:,:]
        anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0,vmax=200)
    elif i==5:
        plt.title(title1[i]+' km')
        plot_data = data_new[5,:,:]
        anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0,vmax=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="7%", pad=0.3)
    plt.colorbar(orientation = 'horizontal',cax=cax)
    anth_plot.cmap.set_bad('lightgrey')
plt.show()
