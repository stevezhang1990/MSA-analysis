import os
import numpy as np
from get_intmap import *
import glob
from netCDF4 import Dataset
from scipy import *
from pylab import *

infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_lat,len_lon = shape(infile00.variables["CHEM-L_S__OH"][0,:,:])
path = '/data/ctm/GEOS_4x5/GEOS5/'
year = np.array([2009,2010])
month = np.array([11,7])
month_counter=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
ps = np.zeros((24*26,len_lat,len_lon))
day0 = np.array([2,-13])
for month_range in range(2):
   for day in arange(2,15):
      YYYY = str(year[month_range]).zfill(4)
      MM = str(month[month_range]).zfill(2)
      DD = str(day).zfill(2)
      filename = path+YYYY+'/'+MM+'/'+'GEOS_5.'+YYYY+MM+DD+'.I3.4x5.nc'
      #print filename
      infile = Dataset(filename)
      for k in range(24):
         I3_i = int(k/3)
         ts_i = k + 24*(day-day0[month_range])
         #print ts_i
         ps[ts_i,:,:] = infile.variables["PS"][I3_i,:,:]
      
np.save('/users/jk/15/xzhang/TS_annual/geos5_ap_msa1/ts_ps_msa1.npy',ps)
print "done"
