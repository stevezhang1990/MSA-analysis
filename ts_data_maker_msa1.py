#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os

file_count = 0
infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_M_V34/runs/v8-02-01/geos5_mop_0911/ctm.01.20091102.nc')
len_lvl, len_lat, len_lon = shape(infile00.variables["IJ-AVG-S__CO"][:,:,:])
#len_time = 24*366
ts_day_fin = 26*24
len_spec = 3
data_in = np.zeros((ts_day_fin,len_spec,len_lvl,len_lat,len_lon))
#mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
year = np.array([2009,2010])
month = np.array([11,7])

len_atime = 26
day0 = np.array([2,2-13])
path1 = '/users/jk/16/xzhang/TS_annual/geos5_omino2_msa1/'
print "check"

#MM = str(month[month_range]).zfill(2)
#month = str(i+1).zfill(2)
#index1 = str(5).zfill(2)
#month2 = str(i+21).zfill(2)
#day_counter = np.sum(mon_counter[0:i])
ts_day_ini = 0
#spec_ap = np.load(path+'ts_ap_16'+month+'.npy')
#print "ap"
spec_in = np.load(path1+'ts_omino2_0911.npy')
spec_in2 = np.load(path1+'ts_omino2_1007.npy')
#len_time_ap= len(spec_ap[:,0,0,0,0])
len_time_in = len(spec_in[:,0,0,0,0])
len_time_in2 = len(spec_in2[:,0,0,0,0])
#data_ap[ts_day_ini:ts_day_fin,:,:,:,:] = spec_ap
data_in[0:ts_day_ini+len_time_in,:,:,:,:] = spec_in
data_in[ts_day_fin-len_time_in2:ts_day_fin,:,:,:,:] = spec_in2

dataset = Dataset(path1+'ts_omino2_msa1.nc','w',format='NETCDF4_CLASSIC')
level = dataset.createDimension('level',len_lvl)
lat = dataset.createDimension('lat',len_lat)
lon = dataset.createDimension('lon',len_lon)
time = dataset.createDimension('time',ts_day_fin)
ts_co = dataset.createVariable('ts_co',np.float32,('time','level','lat','lon'))
ts_co[:,:,:,:] = data_in[:,0,:,:,:]
ts_no2 = dataset.createVariable('ts_no2',np.float32,('time','level','lat','lon'))
ts_no2[:,:,:,:] = data_in[:,1,:,:,:]
ts_o3 = dataset.createVariable('ts_o3',np.float32,('time','level','lat','lon'))
ts_o3[:,:,:,:] = data_in[:,2,:,:,:]
dataset.close()
print "npy to netCDF conversion complete"    

# for the annual file, the shape is time,species#,ver_lvl,lat,lon
# species #: 0:CO 2:O3 1:NO2
print "done"
show()
