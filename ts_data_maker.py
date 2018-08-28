#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset

infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_T_V34/runs/v8-02-01/geos5_all_0911_05/ctm.01.20091102.nc')
len1,len2,len3 = shape(infile00.variables["IJ-AVG-S__CO"][:38,:,:])
strat_co = infile00.variables["IJ-AVG-S__CO"][38:,:,:]
strat_no2 = infile00.variables["IJ-AVG-S__NOx"][38:,:,:]
strat_o3 = infile00.variables["IJ-AVG-S__O3"][38:,:,:]
len_hour = 24
#change this
len_day = 16
#change this
day0 = 16
#change this
month0 = 1#7
year0 = 2017#2016
len4 = len_hour*len_day
len5 = 5
ts = np.zeros((len4,len5,len1,len2,len3))
#change this
#home_dir = '/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1607_28/ND49/'
home_dir = '/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1721/ND49/'
len_trop = 38
for i in arange(len_day):
    for j in range(len_hour):
        day = str(i+day0).zfill(2)
        hour = str(j).zfill(2)
        month = str(month0).zfill(2)
        year = str(year0)
        filename = 'ts.'+year+month+day+'.'+hour+'0000.nc'
        dir_name = home_dir+filename
        infile = Dataset(dir_name)
        len4_count = 24*i+j
        ts[len4_count,0,:,:,:] = infile.variables["IJ-AVG-S__CO"][:,:,:]
        ts[len4_count,1,:,:,:] = infile.variables["IJ-AVG-S__CH2O"][:,:,:]
        ts[len4_count,2,:,:,:] = infile.variables["TIME-SER__NO2"][:,:,:]
        ts[len4_count,3,:,:,:] = infile.variables["IJ-AVG-S__O3"][:,:,:]
        #ts[len4_count,4,:,:,:] = infile.variables["TIME-SER__OH"][:,:,:]
        #ts[len4_count,0,len_trop:,:,:] = strat_co
        #ts[len4_count,1,:,:,:] = infile.variables["IJ-AVG-S__CH2O"][:,:,:]
        #ts[len4_count,1,len_trop:,:,:] = strat_no2
        #ts[len4_count,2,len_trop:,:,:] = strat_o3
        #ts[len4_count,0,:,:,:] = infile.variables["TIME-SER__AIRDEN"][:,:,:]
        #ts[len4_count,1,:,:,:] = infile.variables["BXHGHT-S__BXHEIGHT"][:,:,:]
#change this
outfile = 'ts_omino2_1721.npy'
#the ts npy file includes matrix in shape of time,species#,ver_lvl,lat,lon
#species#: 0:CO,1:CH2O,2:NO2,3:O3
dir_out = home_dir+outfile
np.save(dir_out,ts)
print "done"
show()
