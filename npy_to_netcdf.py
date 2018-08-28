import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

def npy_2_netcdf(path,filename):
   ts = np.load(path+filename+'.npy')
   len_time,len_spec,len_trop,len_lat,len_lon = shape(ts)
   dataset = Dataset(path+filename+'.nc','w',format='NETCDF4_CLASSIC')
   level = dataset.createDimension('level',len_trop)
   lat = dataset.createDimension('lat',len_lat)
   lon = dataset.createDimension('lon',len_lon)
   time = dataset.createDimension('time',len_time)
   ts_co = dataset.createVariable('ts_co',np.float32,('time','level','lat','lon'))
   ts_co[:,:,:,:] = ts[:,0,:,:,:]
   ts_ch2o = dataset.createVariable('ts_ch2o',np.float32,('time','level','lat','lon'))
   ts_ch2o[:,:,:,:] = ts[:,1,:,:,:]
   ts_no2 = dataset.createVariable('ts_no2',np.float32,('time','level','lat','lon'))
   ts_no2[:,:,:,:] = ts[:,2,:,:,:]
   ts_o3 = dataset.createVariable('ts_o3',np.float32,('time','level','lat','lon'))
   ts_o3[:,:,:,:] = ts[:,3,:,:,:]
   dataset.close()
   print "npy to netCDF conversion complete"
   return 

path = '/users/jk/15/xzhang/TS_annual/geos5_ap_2016/'
path1 = '/users/jk/15/xzhang/TS_annual/geos5_mop_2016/'
path2 = '/users/jk/15/xzhang/TS_annual/geos5_all_05_2016/'
path3 = '/users/jk/15/xzhang/TS_annual/geos5_all_06_2016/'
filename = 'ts_ap_2016'
filename1 = 'ts_mop_2016'
filename2 = 'ts_all_05_2016'
filename3 = 'ts_all_06_2016'

npy_2_netcdf(path3,filename3)
