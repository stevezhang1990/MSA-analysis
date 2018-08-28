#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os

file_count = 0
path1 = '/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/'
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_lvl, len_lat, len_lon = shape(infile00.variables["CHEM-L_S__OH"][:,:,:])
len_time = 24*(366+31+28)
len1_time = 24*366
len_spec = 4

data_ts = np.zeros((len_time,len_spec,len_lvl,len_lat,len_lon))
ts0 = Dataset(path1+'ts_all_07_2016.nc')
data_ts[:len1_time,0,:,:,:] = ts0.variables["ts_co"][:,:,:,:]
data_ts[:len1_time,1,:,:,:] = ts0.variables["ts_ch2o"][:,:,:,:]
data_ts[:len1_time,2,:,:,:] = ts0.variables["ts_no2"][:,:,:,:]
data_ts[:len1_time,3,:,:,:] = ts0.variables["ts_o3"][:,:,:,:]

mon_counter = np.array([31,28,31,30,31,30,31,31,30,31,30,31])

#path = '/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_ap_2016/ND49/'
print "check"
for i in range(2):
    month = str(i+1).zfill(2)
    index1 = str(7).zfill(2)
    month2 = str(i+21).zfill(2)
    index2 = str(27).zfill(2)
    day_counter = np.sum(mon_counter[0:i])
    ts_day_ini = 24*day_counter+len1_time
    ts_day_fin = 24*(day_counter+mon_counter[i])+len1_time
    #spec_ap = np.load(path1+'ts_ap_17'+month+'.npy')
    #spec_in = np.load(path1+'ts_omino2_17'+month+'.npy')
    #spec_in2 = np.load(path1+'ts_omino2_17'+month2+'.npy')
    #print "ap"
    spec_in = np.load(path1+'ts_all_17'+month+'_'+index1+'.npy')
    spec_in2 = np.load(path1+'ts_all_17'+month+'_'+index2+'.npy')
    print month
    #len_time_ap= len(spec_ap[:,0,0,0,0])
    len_time_in = len(spec_in[:,0,0,0,0])
    len_time_in2 = len(spec_in2[:,0,0,0,0])
    #data_ts[ts_day_ini:ts_day_fin,:,:,:,:] = spec_ap[:,:4,:,:,:]
    data_ts[ts_day_ini:ts_day_ini+len_time_in,:,:,:,:] = spec_in[:,:4,:,:,:]
    data_ts[ts_day_fin-len_time_in2:ts_day_fin,:,:,:,:] = spec_in2[:,:4,:,:,:]

dataset = Dataset(path1+'ts_all_07_msa2.nc','w',format='NETCDF4_CLASSIC')
level = dataset.createDimension('level',len_lvl)
lat = dataset.createDimension('lat',len_lat)
lon = dataset.createDimension('lon',len_lon)
time = dataset.createDimension('time',len_time)
ts_co = dataset.createVariable('ts_co',np.float32,('time','level','lat','lon'))
ts_co[:,:,:,:] = data_ts[:,0,:,:,:]
ts_ch2o = dataset.createVariable('ts_ch2o',np.float32,('time','level','lat','lon'))
ts_ch2o[:,:,:,:] = data_ts[:,1,:,:,:]
ts_no2 = dataset.createVariable('ts_no2',np.float32,('time','level','lat','lon'))
ts_no2[:,:,:,:] = data_ts[:,2,:,:,:]
ts_o3 = dataset.createVariable('ts_o3',np.float32,('time','level','lat','lon'))
ts_o3[:,:,:,:] = data_ts[:,3,:,:,:]
dataset.close()
print "npy to netCDF conversion complete"    

# for the annual file, the shape is time,species#,ver_lvl,lat,lon
# species #: 0:CO 3:O3 2:NO2 1:CH2O
print "done"
show()
