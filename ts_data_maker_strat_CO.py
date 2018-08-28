import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *



CO_001 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_ap_2016/ctm.00.20160101.nc')

len_time = 366
len_strat,len_lat,len_lon = shape(CO_001.variables["IJ-AVG-S__CO"][38:,:,:])
mon_counter=np.array([15,16,15,14,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
CO_strat = np.zeros((4,len_time,len_strat,len_lat,len_lon))
for exp_i in range(4):
   for i in range(12):
      month = str(i+1).zfill(2)
      month2 = str(i+21).zfill(2)
      if exp_i == 0:
         path_in = '/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_ap_2016/'
         filename = 'ctm.00.2016'+month+'01.nc'
         filename2 = 'ctm.00.2016'+month+'16.nc'
      elif exp_i ==1:
         path_in = '/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/'
         filename = 'geos5_mop_16'+month+'/ctm.10.2016'+month+'01.nc'
         filename2 = 'geos5_mop_16'+month2+'/ctm.10.2016'+month+'16.nc'
      elif exp_i ==2:
         path_in = '/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/'
         filename = 'geos5_all_16'+month+'_07/ctm.20.2016'+month+'01.nc'
         filename2 = 'geos5_all_16'+month+'_27/ctm.20.2016'+month+'16.nc'
         if os.path.isfile(path_in+filename) == False:
            filename = 'geos5_all_16'+month+'_07/ctm.15.2016'+month+'01.nc'
         if os.path.isfile(path_in+filename2) == False:
            filename2 = 'geos5_all_16'+month+'_27/ctm.15.2016'+month+'16.nc'   
            
      elif exp_i ==3:
         path_in = '/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/'
         filename = 'geos5_all_16'+month+'_08/ctm.20.2016'+month+'01.nc'
         filename2 = 'geos5_all_16'+month+'_28/ctm.20.2016'+month+'16.nc'
         if os.path.isfile(path_in+filename) == False:
            filename = 'geos5_all_16'+month+'_08/ctm.15.2016'+month+'01.nc'
         if os.path.isfile(path_in+filename2) == False:
            filename2 = 'geos5_all_16'+month+'_28/ctm.15.2016'+month+'16.nc'          
      count1 = 2*i
      count2 = 2*i+1
      count3 = 2*i+2
      count1_ini = np.sum(mon_counter[0:count1])
      count1_fin = np.sum(mon_counter[0:count2])
      count2_ini = np.sum(mon_counter[0:count2])
      count2_fin = np.sum(mon_counter[0:count3])
      print path_in+filename
      print path_in+filename2
      ncfile = Dataset(path_in+filename)
      ncfile2 = Dataset(path_in+filename2)
      for ts_i in range(count1_ini,count1_fin):
         CO_strat[exp_i,ts_i,:,:,:] = ncfile.variables["IJ-AVG-S__CO"][38:,:,:]
      for ts2_i in range(count2_ini,count2_fin):
         CO_strat[exp_i,ts2_i,:,:,:] = ncfile2.variables["IJ-AVG-S__CO"][38:,:,:]
         
print "check"
path1 = '/users/jk/16/xzhang/TS_annual/geos5_ap_2016/co_strat_ap_2016.nc'
path2 = '/users/jk/16/xzhang/TS_annual/geos5_mop_2016/co_strat_mop_2016.nc'
path3 = '/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/co_strat_msa07_2016.nc'
path4 = '/users/jk/16/xzhang/TS_annual/geos5_all_08_2016/co_strat_msa08_2016.nc'

for exp_i in range(4):
   if exp_i ==0:
      filename_out = path1
   elif exp_i ==1:
      filename_out = path2
   elif exp_i ==2:
      filename_out = path3
   elif exp_i ==3:
      filename_out = path4
   dataset = Dataset(filename_out,'w',format='NETCDF4_CLASSIC')
   lat_nc = dataset.createDimension('lat',len_lat)
   lon_nc = dataset.createDimension('lon',len_lon)
   time_nc = dataset.createDimension('time',len_time)
   level_nc = dataset.createDimension('strat_level',len_strat)
   if exp_i ==0:
      co_ap_strat = dataset.createVariable('CO_ap_strat',np.float32,('time','strat_level','lat','lon'))
      co_ap_strat[:,:,:,:] = CO_strat[0,:,:,:,:]
   elif exp_i ==1:
      co_mop_strat = dataset.createVariable('CO_mop_strat',np.float32,('time','strat_level','lat','lon'))
      co_mop_strat[:,:,:,:] = CO_strat[1,:,:,:,:]
   elif exp_i ==2:
      co_msa05_strat = dataset.createVariable('CO_msa05_strat',np.float32,('time','strat_level','lat','lon'))
      co_msa05_strat[:,:,:,:] = CO_strat[2,:,:,:,:]
   elif exp_i ==3:
      co_msa06_strat = dataset.createVariable('CO_msa06_strat',np.float32,('time','strat_level','lat','lon'))
      co_msa06_strat[:,:,:,:] = CO_strat[3,:,:,:,:]
   dataset.close()
print "done"
