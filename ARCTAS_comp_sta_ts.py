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

dataset = np.load('/users/jk/16/xzhang/ARCTAS/arctas_gc_reshape.npy')
len_exp,len_lat,len_alt = shape(dataset)
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
#o3_gc = get_chart(data)

diff = np.zeros_like(data)
for i in range(len_exp):
    diff_raw = data[4,:,:] - data[i,:,:] 
    #print np.max(data[i,:,:,:]), np.min(data[i,:,:,:])
    #diff_filt = np.ma.masked_where(abs(diff_raw) > 199, diff_raw)
    #diff[i,:,:,:] = np.ma.masked_where(abs(diff_raw) == 0, diff_raw)
    diff[i,:,:] = diff_raw
diff = ma.masked_invalid(diff)
diff_abs = diff**2
delta = 5
latitude = np.arange(10,70.5,1)
altitude = np.arange(0,12.6,0.2)
#chart_mean = get_chart(diff)
#chart_rmse = (get_chart(diff_abs))**0.5
title2 = ['CTRL','OMI','OMO','OMO BC']
clevs = np.arange(-200,200,delta)
# Draw Contours
fig = plt.figure(1,figsize=(14,9))
for j in range(len_exp-1):
    ax = subplot(221+j)
    CS = contourf(latitude,altitude,np.transpose(diff[j,:,:]),clevs,cmap=plt.cm.RdBu_r,extend='both',\
                      extendfrac='auto')
    plt.title(r'NO$_2$ diff between ARCTAS and '+title2[j])
    plt.colorbar(mappable=None, cax=None, ax=None)
show()
