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

#len1, len2 = np.shape(file_atom)
atom_file = ['atom1','atom2']
file_i = 0
file_atom1 = loadtxt('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_1.csv')
file_atom2 = loadtxt('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_2.csv')
file_atom = np.vstack((file_atom1,file_atom2))
len3, len2 = np.shape(file_atom)
file_post = Dataset('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_gc.nc')
len1,len3,len5 = shape(file_post.variables["atom_gc_obs"][:,:,:])
atom_spec = np.zeros((len1+1,len3,len5))
atom_spec[:7,:,:] = file_post.variables["atom_gc_obs"][:,:,:]
 
year, month, day, time, lon, lat, alt, oh_atom, o3_noy, o3_ucats, co_geos, co_ucats, ch2o_atom, no2_atom = \
    file_atom[:,0], file_atom[:,1], file_atom[:,2], file_atom[:,3], file_atom[:,4], file_atom[:,5],\
    file_atom[:,6], file_atom[:,7], file_atom[:,8], file_atom[:,9], file_atom[:,10], file_atom[:,11], \
    file_atom[:,12], file_atom[:,13] 
#OH,O3,CO,CH2O,NO2
atom_spec[7,:,0] = oh_atom
atom_spec[7,:,1] = np.mean(np.ma.masked_where(np.ma.masked_invalid(file_atom[:,8:10])>299,\
                                              np.ma.masked_invalid(file_atom[:,8:10])),axis=1)
atom_spec[7,:,2] = np.mean(np.ma.masked_where(np.ma.masked_invalid(file_atom[:,10:12])>299,\
                                              np.ma.masked_invalid(file_atom[:,10:12])),axis=1)
atom_spec[7,:,3] = np.ma.masked_where(np.ma.masked_invalid(file_atom[:,12])>299,\
                                              np.ma.masked_invalid(file_atom[:,12]))
atom_spec[7,:,4] = np.ma.masked_where(np.ma.masked_invalid(file_atom[:,13])>299,\
                                              np.ma.masked_invalid(file_atom[:,13]))

latitude = np.arange(-90,91,1)
altitude = np.arange(0,14000,200)
spec_obs = np.zeros((8,3,2,len(latitude),len(altitude),len5))
for i in range(len3):
    for j in range((len(latitude))-1):
        for k in range((len(altitude))-1):
            if lat[i]>latitude[j] and abs(lat[i])<=90:
                if lat[i]<latitude[j+1]:
                    if alt[i]>altitude[k]:
                        if alt[i]<altitude[k+1] and abs(lon[i])<=180:
                               #print lon[i]
                            if lon[i] > -75 and lon[i]<20:
                                crt_lon = 0
                            else:
                                crt_lon = 1
                            for spec_i in range(len5):
                                if atom_spec[7,i,spec_i] < 299 and atom_spec[7,i,spec_i]>0:
                                    spec_obs[:,0,crt_lon,j,k,spec_i] = spec_obs[:,0,crt_lon,j,k,spec_i]+atom_spec[:,i,spec_i]
                                    spec_obs[:,1,crt_lon,j,k,spec_i] = spec_obs[:,1,crt_lon,j,k,spec_i]+1
                                   #if o3_ucats[i] < 299 and o3_ucats[i]>5:
                                       #o3_obs[j,k] = o3_obs[j,k]+oh_atom[i]#100*(o3_diff2[i]/oh_atom[i])
                                       #o3_obs[crt_lon,j,k] = o3_obs[crt_lon,j,k]+o3_gc[i]
                                       #o3_count[crt_lon,j,k] = o3_count[crt_lon,j,k]+1
spec_obs[:,2,:,:,:,:] = spec_obs[:,0,:,:,:,:]/spec_obs[:,1,:,:,:,:]
#np.save("/users/jk/16/xzhang/TS_annual/geos5_ap_2016/"+atom_file[file_i]+"_gc_reshape.npy",o3_post)
np.save("/users/jk/15/xzhang/ATom/"+atom_file[file_i]+"_merged_reshape.npy",spec_obs[:,2,:,:,:,:])
print "done"                           

