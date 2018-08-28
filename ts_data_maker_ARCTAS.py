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

path = '/users/jk/16/xzhang/ARCTAS/'
infile = Dataset(path+'arctas_no2_200807.nc')
lat = infile.variables["lat_orbit"][:]
lon = infile.variables["lon_orbit"][:]
alt = infile.variables["alt_orbit"][:]
no2_arctas = np.zeros((6,len(lat)))
no2_arctas[0,:] = infile.variables["NO2_ap_NA"][:]
no2_arctas[1,:] = infile.variables["NO2_omi_NA"][:]
no2_arctas[2,:] = infile.variables["NO2_omo_NA"][:]
no2_arctas[3,:] = infile.variables["NO2_omo_BC"][:]
no2_arctas[4,:] = infile.variables["NO2_NCAR_orbit"][:]
no2_arctas[5,:] = infile.variables["NO2_UCB_orbit"][:]
latitude = np.arange(10,70.5,1)
altitude = np.arange(0,12.6,0.2)
no2_obs = np.zeros((5,len(latitude),len(altitude)))
no2_count = np.zeros((5,len(latitude),len(altitude)))
no2_post = np.zeros((5,len(latitude),len(altitude)))
for i in range((len(lat))):
    for j in range((len(latitude))-1):
        for k in range((len(altitude))-1):
            if lat[i]>latitude[j] and abs(lat[i])<=70 and lat[i]>=10:
                if lat[i]<=latitude[j+1]:
                   if alt[i]>altitude[k]:
                       if alt[i]<=altitude[k+1] and lon[i]<=-40 and lon[i]>=-140:
                           for m in range(5):
                               if no2_arctas[4,i] <1e4 and no2_arctas[4,i]>5:
                                   no2_obs[m,j,k] = no2_obs[m,j,k]+no2_arctas[m,i]
                                   no2_count[m,j,k] = no2_count[m,j,k]+1
                               if no2_arctas[5,i] < 1e4 and no2_arctas[5,i]>5:
                                   if m<5:
                                       no2_obs[m,j,k] = no2_obs[m,j,k]+no2_arctas[m,i]
                                   else:
                                       no2_obs[4,j,k] = no2_obs[4,j,k]+no2_arctas[5,i]
                                   no2_count[m,j,k] = no2_count[m,j,k]+1
for j in range((len(latitude))):
    for k in range((len(altitude))):
        for i in range(5):
            if no2_obs[i,j,k] != 0:
                no2_post[i,j,k] = no2_obs[i,j,k]/no2_count[i,j,k]
np.save(path+'arctas_gc_reshape.npy',no2_post)
#np.save("/users/jk/15/xzhang/ATom/"+atom_file[file_i]+"_obs_reshape.npy",no2_post)
print "done"                           

