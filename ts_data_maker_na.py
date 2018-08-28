#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os

infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_1/ctm.05.20080701.nc')
len1,len2,len3 = shape(infile00.variables["IJ-AVG-S__CO"][:,:,:])
len_hour = 24
#change this
len_day = 27
#change this
day0 = 1
#change this
month0 = 7#7
year0 = 2008#2016
len4 = len_hour*len_day
len5 = 5
ts = np.zeros((len4,len1,len2,len3))
path = '/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/'
i=0
#for i in range(9):
if i==0:
    run_ind = str(i+1)
    #home_dir = path+'geos5_omo_0807_BC_'+run_ind+'/ND49/'
    home_dir = path+'geos5_ap_0807_NA/ND49/'
    for filename in glob.glob(os.path.join(home_dir,'*.nc')):
        infile = Dataset(filename)
        filename_str = filename.split('/')
        file_str = filename_str[-1].split('.')
        day = float(file_str[1])-year0*10000-month0*100
        hour = float(file_str[2])/10000
        len4_count = int(24*(day-1)+hour)
        ts[len4_count,:,:,:] = infile.variables["TIME-SER__NO2"][:,:,:]
#change this
outfile = 'ts_ap_0807_NA.npy'
#the ts npy file includes matrix in shape of time,species#,ver_lvl,lat,lon
#species#: 0:CO,1:CH2O,2:NO2,3:O3
out_dir = '/users/jk/16/xzhang/ARCTAS/'
dir_out = out_dir+outfile
np.save(dir_out,ts)
print "done"
show()
