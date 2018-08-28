import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

data = Dataset('/users/jk/16/xzhang/SHADOZ/shadoz_obs_gc_2016.nc')
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_pres_2016.nc')
len_time,len_lvl,len_lat,len_lon = shape(data.variables["shadoz_o3_obs"][:,:,:,:])
len_time1 = len_time#-6*24
gc_shadoz = np.zeros((5,len_time1,len_lvl,len_lat,len_lon))
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
gc_shadoz[0,:,:,:,:] = data.variables["shadoz_o3_ap"][0:len_time1,:,:,:]
gc_shadoz[1,:,:,:,:] = data.variables["shadoz_o3_iasio3"][0:len_time1,:,:,:]
gc_shadoz[2,:,:,:,:] = data.variables["shadoz_o3_msa07"][0:len_time1,:,:,:]
gc_shadoz[3,:,:,:,:] = data.variables["shadoz_o3_msa08"][0:len_time1,:,:,:]
gc_shadoz[4,:,:,:,:] = data.variables["shadoz_o3_obs"][0:len_time1,:,:,:]
ts_p_center = np.mean(np.ma.masked_invalid(ts_pres.variables["ts_p_center"][:,:,:,:]),axis=0)
o3_shadoz = np.zeros((30,5,len_time1,len_lvl))
lat_obs = np.zeros(30)
lon_obs = np.zeros(30)
diff_obs = np.zeros(len_lvl)
plot_counter = 0
#gc_shadoz = np.ma.masked_invalid(gc_shadoz)
for i in range(len_lat):
   for j in range(len_lon):
      if np.sum(gc_shadoz[1,:,:,i,j])!=0:
         o3_shadoz[plot_counter,:,:,:] = gc_shadoz[:,:,:,i,j]
         lat_obs[plot_counter] = lat[i]
         lon_obs[plot_counter] = lon[j]
         plot_counter = plot_counter+1
dataset = Dataset('/users/jk/16/xzhang/SHADOZ/shadoz_obs_2016.nc','w',format='NETCDF4_CLASSIC')
len_station = dataset.createDimension('len_station',plot_counter)
len_level = dataset.createDimension('len_lvl',len_lvl)
len_time = dataset.createDimension('len_time',len_time)
lat_shadoz = dataset.createVariable('lat_obs',np.float32,('len_station'))
lat_shadoz[:] = lat_obs[:plot_counter] 
lon_shadoz = dataset.createVariable('lon_obs',np.float32,('len_station'))
lon_shadoz[:] = lon_obs[:plot_counter]
O3_shadoz = dataset.createVariable('shadoz_o3_obs',np.float32,('len_station','len_time','len_lvl'))
O3_shadoz[:,:,:] = o3_shadoz[:plot_counter,4,:,:]
O3_ap_sha = dataset.createVariable('shadoz_o3_ap',np.float32,('len_station','len_time','len_lvl'))
O3_ap_sha[:,:,:] = o3_shadoz[:plot_counter,0,:,:]
O3_iasio3_sha = dataset.createVariable('shadoz_o3_iasio3',np.float32,('len_station','len_time','len_lvl'))
O3_iasio3_sha[:,:,:] = o3_shadoz[:plot_counter,1,:,:]
O3_msa07_sha = dataset.createVariable('shadoz_o3_msa07',np.float32,('len_station','len_time','len_lvl'))
O3_msa07_sha[:,:,:] = o3_shadoz[:plot_counter,2,:,:]
O3_msa08_sha = dataset.createVariable('shadoz_o3_msa08',np.float32,('len_station','len_time','len_lvl'))
O3_msa08_sha[:,:,:] = o3_shadoz[:plot_counter,3,:,:]

dataset.close()
print "done"         
print "check"
show()

