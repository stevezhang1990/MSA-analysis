#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os

file_count = 0
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
path = '/users/jk/16/xzhang/AQS/'
aqs_data = np.load(path+'aqs_2010.npy')

#len_time,len_lat,len_lon = shape(ts_ap.variables["ts_o3"][:,0,:,:])
#ts_species = np.zeros((4,3,len_time,len_lat,len_lon))
#ts_species[0,1,:,:,:] = ts_ap.variables["ts_co"][:,0,:,:]
#ts_species[1,1,:,:,:] = ts_mop.variables["ts_co"][:,0,:,:]
#ts_species[2,1,:,:,:] = ts_msa05.variables["ts_co"][:,0,:,:]
#ts_species[3,1,:,:,:] = ts_msa06.variables["ts_co"][:,0,:,:]
#ts_species[0,0,:,:,:] = ts_ap.variables["ts_o3"][:,0,:,:]
#ts_species[1,0,:,:,:] = ts_mop.variables["ts_o3"][:,0,:,:]
#ts_species[2,0,:,:,:] = ts_msa05.variables["ts_o3"][:,0,:,:]
#ts_species[3,0,:,:,:] = ts_msa06.variables["ts_o3"][:,0,:,:]
#ts_species[0,2,:,:,:] = ts_ap.variables["ts_no2"][:,0,:,:]
#ts_species[1,2,:,:,:] = ts_mop.variables["ts_no2"][:,0,:,:]
#ts_species[2,2,:,:,:] = ts_msa05.variables["ts_no2"][:,0,:,:]
#ts_species[3,2,:,:,:] = ts_msa06.variables["ts_no2"][:,0,:,:]
len_lat, len_lon = shape(infile00.variables["DXYP__DXYP"][0,:,:])
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
len_lat_NA = 16
len_lon_NA = 17
# 10-70N, 140-60 W
time_orb,lat_orb,lon_orb = shape(aqs_data)
aqs_gc = np.zeros((time_orb,len_lat_NA,len_lon_NA))
#for counter in range(spe_orb):
for i in range(len_lat_NA):
    for j in range(len_lon_NA):
        for n in range(time_orb):
            aqs_gc_acc = 0
            aqs_gc_count = 0
            for k in range(4):
                for l in range(5):
                    if 4*i+k < lat_orb and 5*j+l < lon_orb:
                        if aqs_data[n,4*i+k,5*j+l] > 0:
                            aqs_gc_acc = aqs_gc_acc + aqs_data[n,4*i+k,5*j+l] 
                            aqs_gc_count = aqs_gc_count + 1
            if aqs_gc_count > 0:
                aqs_gc[n,i,j] = aqs_gc_acc/aqs_gc_count
o3_aqs = aqs_gc[:,:,:]
#no2_aqs = aqs_gc[:,:,:,2]
#co_aqs = aqs_gc[:,:,:,0]

dataset = Dataset(path+'aqs_gc_2010.nc','w',format='NETCDF4_CLASSIC')
lat_nc = dataset.createDimension('lat_NA',len_lat_NA)
lon_nc = dataset.createDimension('lon_NA',len_lon_NA)
time_nc = dataset.createDimension('time',time_orb)
#ts_co_aqs = dataset.createVariable('ts_co_aqs',np.float32,('time','lat_NA','lon_NA'))
#ts_co_aqs[:,:,:] = co_aqs[:,:,:]
ts_o3_aqs = dataset.createVariable('ts_o3_aqs',np.float32,('time','lat_NA','lon_NA'))
ts_o3_aqs[:,:,:] = o3_aqs[:,:,:]
#ts_no2_aqs = dataset.createVariable('ts_no2_aqs',np.float32,('time','lat_NA','lon_NA'))
#ts_no2_aqs[:,:,:] = no2_aqs[:,:,:]
dataset.close()


#data_filt = np.zeros((len_time,len_lat,len_lon,len_spec))
#mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
#path = '/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_ap_2016/ND49/'
#for filename in glob.glob(os.path.join(path, '*.npy')):
    #infile = np.load(filename)
    #len1,len2,len3,len4,len5 = shape(infile)
    #date_num0 = np.sum(mon_counter[0:file_count])
    #j = 24*date_num0
    #for i in range(len1):
        #data_filt[i+j,:,:,0] = infile[i,0,0,:,:]
        #data_filt[i+j,:,:,1] = infile[i,3,0,:,:]
        #data_filt[i+j,:,:,2] = infile[i,2,0,:,:]
    #file_count = file_count + 1

#np.save('/users/jk/15/xzhang/AQS/aqs_gc_2016.npy',data_filt)
print "done"
show()
