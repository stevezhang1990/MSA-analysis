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
title3 = ['ATom-1','ATom-2']
file_i = 1
dataset = np.load('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_merged_reshape.npy')
#ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_pres_msa1.nc')
#airden_gc = ts_pres.variables["ts_air_density"][0,:,:,:]

data2 = np.ma.masked_invalid(dataset)
data = np.ma.masked_where(abs(data2)==0, data2)
data[:7,:,:,:,0] = 0.5*data[:7,:,:,:,0]
data[:,:,:,:,0] = 1000*data[:,:,:,:,0]
len1,len2,len3,len4,len5 = shape(data)
spec_mean = np.zeros((len1,len2,3,len4,len5))
spec_mean[:,:,0,:,:] = np.mean(data[:,:,15:60,:,:],axis=2)
spec_mean[:,:,1,:,:] = np.mean(data[:,:,60:120,:,:],axis=2)
spec_mean[:,:,2,:,:] = np.mean(data[:,:,120:165,:,:],axis=2)
spec_mean_trop = np.ma.masked_where(abs(spec_mean[:,:,:,5:50,:])==0, spec_mean[:,:,:,5:50,:])
spec_mean_filt = np.around(np.mean(np.ma.masked_where(abs(spec_mean_trop)>0.8*np.max(spec_mean_trop[7,:,:,:,:]),\
                                                      spec_mean_trop),axis=3),decimals=1)                                    
spec_plot = np.zeros((len1,len2,13,len5))
for i in range(13):
    if i<12:
        spec_plot[:,:,i,:] = np.mean(np.mean(data[:,:,:,5*i:5*(i+1),:],axis=3),axis=2)
    else:
        spec_plot[:,:,i,:] = np.mean(np.mean(data[:,:,:,60:,:],axis=3),axis=2)
props = dict(boxstyle='square', facecolor='white', alpha=0.0)
title1 = ['Atlantic ','Pacific ']
title2 = ['OH',r'$O_3$','CO','HCHO',r'$NO_2$']
title4 = ['MOP','IASI','MOP','OMI','OMI']
altitude = np.arange(0,13,1)
#marker2 = ['v','o','^']
marker1 = ['k','b','g','r']
exp_ind = np.array([[1,1,5,7],[1,2,6,7],[2,1,5,7],[1,4,5,7],[1,3,5,7]])
fig = plt.figure(1,figsize=(12,6.5))
for j in range(len2):
    for k in arange(len5-1):
        ax = subplot(241+k+4*j)
        title(title3[file_i]+' '+title2[k]+' '+title1[j])
        #for i in range(3):
        for l in range(4):
            plot(spec_plot[int(exp_ind[k,l]),j,:,k],altitude,marker1[l])
        textstr = '\n'.join((
            '75-30N, '+'30N-30S, '+'30S-75S',
            'CTRL='+str(spec_mean_filt[int(exp_ind[k,0]),j,2,k])+', '+str(spec_mean_filt[int(exp_ind[k,0]),j,1,k])\
            +', '+str(spec_mean_filt[int(exp_ind[k,0]),j,0,k]),
            title4[k]+'='+str(spec_mean_filt[int(exp_ind[k,1]),j,2,k])+', '+str(spec_mean_filt[int(exp_ind[k,1]),j,1,k])\
            +', '+str(spec_mean_filt[int(exp_ind[k,1]),j,0,k]),
            'MSA='+str(spec_mean_filt[int(exp_ind[k,2]),j,2,k])+', '+str(spec_mean_filt[int(exp_ind[k,2]),j,1,k])\
            +', '+str(spec_mean_filt[int(exp_ind[k,2]),j,0,k]),
            'OBS='+str(spec_mean_filt[int(exp_ind[k,3]),j,2,k])+', '+str(spec_mean_filt[int(exp_ind[k,3]),j,1,k])\
            +', '+str(spec_mean_filt[int(exp_ind[k,3]),j,0,k])))
        ax.text(0.98, 0.05, textstr, transform=ax.transAxes, fontsize=8,horizontalalignment='right',
                 verticalalignment='bottom', bbox=props)
show()
