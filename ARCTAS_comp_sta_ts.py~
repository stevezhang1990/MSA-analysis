#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
#to plot the spaghetti, make sure you prepare a reshape file first.
atom_file = ['atom1','atom2']
file_i = 0
data_ap = np.load('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/'+atom_file[file_i]+'_gc_reshape.npy')
len_spec,len_lat,len_alt = shape(data_ap)
len_exp = 5
dataset = np.zeros((len_exp,len_spec,len_lat,len_alt))
dataset[0,:,:,:] = data_ap
dataset[1,:,:,:] = np.load('/users/jk/16/xzhang/TS_annual/geos5_iasio3_2016/'+atom_file[file_i]+'_gc_reshape.npy')
dataset[3,:,:,:] = np.load('/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/'+atom_file[file_i]+'_gc_reshape.npy')
dataset[2,:,:,:] = np.load('/users/jk/16/xzhang/TS_annual/geos5_all_08_2016/'+atom_file[file_i]+'_gc_reshape.npy')
dataset[4,:,:,:] = np.load('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_obs_reshape.npy')
for i in range(131-45):
    for j in range(40):
        if dataset[2,1,45+i,j] > 0:
            dataset[2,1,45+i,j] = dataset[2,1,45+i,j] + 5
def get_chart(data_in):
    len_exp1,len_spec1,len_alt1,len_lon1 = shape(data_in)
    lvl_crt = int(len_alt1/2)
    data_in_chart = np.zeros((len_exp1,len_spec1,13))
    data_in_chart[:,:,0] = np.mean(np.mean(data_in[:,:,:,:],axis=2),axis=2)
    data_in_chart[:,:,1] = np.mean(np.mean(data_in[:,:,0:lvl_crt,:],axis=2),axis=2)
    data_in_chart[:,:,2] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,:],axis=2),axis=2)
    data_in_chart[:,:,3] = np.mean(np.mean(data_in[:,:,0:lvl_crt,150:180],axis=2),axis=2)
    data_in_chart[:,:,4] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,150:180],axis=2),axis=2)
    data_in_chart[:,:,5] = np.mean(np.mean(data_in[:,:,0:lvl_crt,120:150],axis=2),axis=2)
    data_in_chart[:,:,6] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,120:150],axis=2),axis=2)
    data_in_chart[:,:,7] = np.mean(np.mean(data_in[:,:,0:lvl_crt,60:120],axis=2),axis=2)
    data_in_chart[:,:,8] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,60:120],axis=2),axis=2)
    data_in_chart[:,:,9] = np.mean(np.mean(data_in[:,:,0:lvl_crt,30:60],axis=2),axis=2)
    data_in_chart[:,:,10] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,30:60],axis=2),axis=2)
    data_in_chart[:,:,11] = np.mean(np.mean(data_in[:,:,0:lvl_crt,0:30],axis=2),axis=2)
    data_in_chart[:,:,12] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,0:30],axis=2),axis=2)
    return data_in_chart
#data2 = np.ma.masked_invalid(dataset)
#data = np.ma.masked_where(abs(data2)==0, data2)
#data = np.ma.masked_where(abs(dataset)==5, dataset)
data = dataset
o3_gc = get_chart(data)

diff = np.zeros_like(data)
for i in range(len_exp):
    diff_raw = data[4,:,:,:] - data[i,:,:,:] 
    #print np.max(data[i,:,:,:]), np.min(data[i,:,:,:])
    #diff_filt = np.ma.masked_where(abs(diff_raw) > 199, diff_raw)
    #diff[i,:,:,:] = np.ma.masked_where(abs(diff_raw) == 0, diff_raw)
    diff[i,:,:,:] = diff_raw
diff = ma.masked_invalid(diff)
diff_abs = diff**2
delta = 5
latitude = np.arange(-90,91,1)
altitude = np.arange(0,14000,200)
#chart_mean = get_chart(diff)
#chart_rmse = (get_chart(diff_abs))**0.5
title1 = ['Atlantic','Pacific']
title2 = ['CTRL','IASI O3','All instrument']
clevs = np.arange(-20,20,delta)
# Draw Contours
fig = plt.figure(1,figsize=(14,9))
for i in range(len_exp-2):
    for j in range(len_spec):
        ax = subplot(321+i*2+j)
        CS = contourf(latitude,altitude,np.transpose(diff[i,j,:,:]),clevs,cmap=plt.cm.RdBu_r,extend='both',\
                      extendfrac='auto')
        plt.title('CO diff between '+atom_file[file_i]+' and '+title2[i]+' '+title1[j])
        plt.colorbar(mappable=None, cax=None, ax=None)
show()
