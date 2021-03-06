import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
atom_file = ['atom1','atom2']
day0 = np.array([16*24,0])
file_i = 0
file_atom1 = loadtxt('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_1.csv')
file_atom2 = loadtxt('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_2.csv')
file_atom = np.vstack((file_atom1,file_atom2))
len1_atom,len2_atom = shape(file_atom)
atom_gc = np.zeros((7,len1_atom,5))
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_pres_msa1.nc')
XNUMOLAIR = 6.022e23/(28.9644e-3)
#OH,O3,CO,CH2O,NO2
#len1_atom1,len2_atom1 = shape(file_atom1)
#len2_atom1,len2_atom2 = shape(file_atom2)
#file_atom = np.zeros((len1_atom1+len2_atom1,len2_atom1))
#file_atom[:len1_atom1,:] = file_atom1
#file_atom[len1_atom1:,:] = file_atom2
airden_gc = ts_pres.variables["ts_air_density"][:,:,:,:] 
year, month, day, time, lon, lat, alt, \
      OH_atom, O3_noy, O3_ucats, co_geos,\
      co_ucats, ch2o_atom, no2_atom = \
      file_atom[:,0], file_atom[:,1], \
      file_atom[:,2], file_atom[:,3], file_atom[:,4],\
      file_atom[:,5], file_atom[:,6], file_atom[:,7],\
      file_atom[:,8], file_atom[:,9], file_atom[:,10],\
      file_atom[:,11], file_atom[:,12], file_atom[:,13]
day_tot = np.array([31,27])
day1_tot = np.array([15,14])
len_trop=38
len_time = 24*day_tot[file_i]
len_time1 = 24*day1_tot[file_i]
len_spec = np.array([4,5])
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
date1 = ['1608','1702']
date2 = ['1628','1722']
len_lvl,len_lat,len_lon = shape(infile00.variables["CHEM-L_S__OH"][:,:,:])
ts_atom = np.zeros((7,len_time,len_spec[file_i],len_lvl,len_lat,len_lon))
rundir = '/users/jk/16/xzhang/TS_annual/geos5_'
expdir = ['ap','mop','iasio3','omino2','omh','all_07','all_08']
expdir2 = [['05','08','25','28'],['07','08','27','28']]
for j in range(7):
    if j==0 and file_i == 0:
        filename = rundir+expdir[j]+'_2016/ts_'+expdir[j]+'_'+date1[file_i]+'.npy'
    elif j<5:
        filename1 = rundir+expdir[j]+'_2016/ts_'+expdir[j]+'_'+date1[file_i]+'.npy'
        filename2 = rundir+expdir[j]+'_2016/ts_'+expdir[j]+'_'+date2[file_i]+'.npy'
    else:
        filename1 = rundir+expdir[j]+'_2016/ts_all_'+date1[file_i]+'_'+expdir2[file_i][j-5]+'.npy'
        filename2 = rundir+expdir[j]+'_2016/ts_all_'+date1[file_i]+'_'+expdir2[file_i][j-3]+'.npy'
    if j==0 and file_i ==0:
        ts_atom[j,:,:,:,:,:] = np.load(filename)
    else:
        ts_atom[j,:len_time1,:,:,:,:] = np.load(filename1)
        ts_atom[j,len_time1:,:,:,:,:] = np.load(filename2)
   
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_pres_msa1.nc')
print "check"
for i in range(len1_atom):
   day_atom = int(file_atom[i,2])
   hour_atom = int(file_atom[i,3]/100)
   ts_i = 24*(day_atom-1)+hour_atom
   lon_atom = file_atom[i,4]
   lat_atom = file_atom[i,5]
   alt_atom = file_atom[i,6]*1e-3
   OH_atom = file_atom[i,7]
   O3_noy = file_atom[i,8]
   O3_ucats = file_atom[i,9]
   co_geos = file_atom[i,10]
   co_ucats = file_atom[i,11]
   ch2o_atom = file_atom[i,12]
   no2_atom = file_atom[i,13]
   lat_i,lon_j = get_gc_ij(lat_atom,lon_atom)
   lvl_k = get_gc_lvl(alt_atom)
   hour_i = day0[file_i]+hour_atom-1
   #print OH_atom,O3_noy,O3_ucats,co_geos,co_ucats,ch2o_atom,no2_atom
   #print lvl_k
   for j in range(7):
       if lat_atom > -90 and lon_atom > -180:
          if O3_noy>0 or O3_ucats>0:
             atom_gc[j,i,1] = ts_atom[j,ts_i,3,lvl_k,lat_i,lon_j]
          if OH_atom>0 and file_i == 1:
             atom_gc[j,i,0] = ts_atom[j,ts_i,4,lvl_k,lat_i,lon_j]*2e12/airden_gc[hour_i,lvl_k,lat_i,lon_j]
          if co_geos>0 or co_ucats>0:
             atom_gc[j,i,2] = ts_atom[j,ts_i,0,lvl_k,lat_i,lon_j]
          if ch2o_atom>0:
             atom_gc[j,i,3] = ts_atom[j,ts_i,1,lvl_k,lat_i,lon_j]*1e3
          if no2_atom>0:
             atom_gc[j,i,4] = ts_atom[j,ts_i,2,lvl_k,lat_i,lon_j]*1e3
          #print atom_gc[j,i,:]
print "check"
dataset = Dataset('/users/jk/15/xzhang/ATom/'+atom_file[file_i]+'_gc.nc','w',format='NETCDF4_CLASSIC')
len_orbit = dataset.createDimension('len_orbit',len1_atom)
len_run = dataset.createDimension('len_run',7)
len_species = dataset.createDimension('len_species',5)            
atom_gc_obs = dataset.createVariable('atom_gc_obs',np.float32,('len_run','len_orbit','len_species'))
atom_gc_obs[:,:,:] = atom_gc[:,:,:]
dataset.close()
show()
