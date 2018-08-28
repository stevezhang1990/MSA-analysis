#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic

path = '/users/jk/16/xzhang/gcadj_std_I_WC/runs/v8-02-01/'
exp = ['1','2','3','4','5','6','7','8','9']
exp2 = ['01','04','07','10','13','16','19','22','25']
infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1608_AM/ctm.00.20160801.nc')
lat_ini = 25
lat_fin = 41
lon_ini = 8
lon_fin = 25
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
lon_gc = infile00.variables["LON"][lon_ini:lon_fin]
lat_gc = infile00.variables["LAT"][lat_ini:lat_fin]
len1,len2,len3 = shape(infile00.variables["IJ-AVG-S__Ox"][:,:,:])
O3_gc = np.zeros((9,2,len1,len2,len3))
O3_crt = np.zeros((9,2,3,len2,len3))
pre_crt = np.array([800,500,300])
#pres = np.zeros((len1,len2,len3))
ps = infile00.variables["DAO-FLDS__PS-PBL"][0,:,:]
for i in range(len(exp)):
    infile = Dataset(path+'geos5_iasio3_0807_'+exp[i]+'/ctm.01.200807'+exp2[i]+'.nc')
    infile2 = Dataset(path+'geos5_iasio3_0807_'+exp[i]+'/ctm.15.200807'+exp2[i]+'.nc')
    O3_gc[i,0,:,:,:] = infile.variables["IJ-AVG-S__Ox"][:,:,:]
    O3_gc[i,1,:,:,:] = infile2.variables["IJ-AVG-S__Ox"][:,:,:]
    for j in range(2):
        for m in range(3):
            pres, O3_crt[i,j,m,:,:] = geos_interpolate(len1,lat,lon,O3_gc[i,j,:,:,:],ps,pre_crt[m],pre_output=True)

O3_avg = np.mean(O3_gc[:,:,:,lat_ini:lat_fin,lon_ini:lon_fin],axis=0)
O3_avg_crt = np.mean(O3_crt[:,:,:,lat_ini:lat_fin,lon_ini:lon_fin],axis=0)
O3_diff_crt = O3_avg_crt[1,:,:,:] - O3_avg_crt[0,:,:,:]
O3_diff_avg = O3_avg[1,:,:,:] - O3_avg[0,:,:,:]
O3_diff_zonal = np.mean(O3_diff_avg,axis=2)
fig = plt.figure(1,figsize=(12,7.5))
bond =5
for i in range(3):
    ax = plt.subplot(221+i)
    m1=Basemap(lon_0=-90,llcrnrlat=10,urcrnrlat=60,llcrnrlon=-140,urcrnrlon=-60)
    m1.drawcoastlines()
    m1.drawparallels(arange(0,90,10),labels=[1,0,0,1],labelstyle="+/-")
    m1.drawmeridians(arange(-180,0,30),labels=[1,0,0,1],labelstyle="+/-")
    title(r'O$_3$ July-2008 IASI O3: '+str(pre_crt[i])+'hPa')
    if i==0:
        plot_data,lon_gc = addcyclic(O3_diff_crt[0,:,:],lon_gc)
        x,y = m1(*np.meshgrid(lon_gc,lat_gc))
    else:
        plot_data = O3_diff_crt[i,:,:]
    anth_plot = m1.pcolormesh(x,y,plot_data,vmin=-bond,vmax=bond,cmap=get_cmap("RdBu_r"))      
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="7%", pad=0.3)
    plt.colorbar(orientation = 'horizontal',cax=cax)
    anth_plot.cmap.set_bad('lightgrey')
ax = plt.subplot(224)
z_len = 47
y = np.mean(np.mean(pres[:,lat_ini:lat_fin,lon_ini:lon_fin],axis=1),axis=1)
delta = 1
clevs = np.arange(-5,5,delta) 
CS = contourf(lat_gc,y,O3_diff_zonal,clevs,cmap=plt.cm.RdBu_r,animated=True,extend='both')
plt.title('DO3 altitude vs latitude IASI O3')
plt.colorbar(mappable=None, cax=None, ax=None)
plt.gca().invert_yaxis()
print "done"
show()
