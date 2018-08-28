import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
ts_sp = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_pres_2016.nc')
ts_bxh = np.load('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_bxh_0911.npy')
ts_bxh1 = np.load('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_bxh_1007.npy')
len_time,len_lvl,len_lat,len_lon = shape(ts_sp.variables["ts_p_center"])
ts_pcenter = ts_sp.variables["ts_p_center"]
ts_pedge = ts_sp.variables["ts_p_edge"]
len_time = 26*24
len_mid = 26*12
len_time_b = 48
mon_counter=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
len_nov_ini = np.sum(mon_counter[:10])+1
len_jul_ini = np.sum(mon_counter[:6])+1
ini11 = 24*len_nov_ini
fin11 = 24*(len_nov_ini+13)
ini07 = 24*len_jul_ini
fin07 = 24*(len_jul_ini+13)
ts_pcenter_msa1 = np.zeros((len_time,len_lvl,len_lat,len_lon))
ts_bxh_msa1 = np.zeros((len_time_b,2,len_lvl,len_lat,len_lon))
ts_pedge_msa1 = np.zeros((len_time,len_lvl+1,len_lat,len_lon))
ts_pcenter_msa1[:len_mid,:,:,:] = ts_pcenter[ini11:fin11,:,:,:]
ts_pcenter_msa1[len_mid:,:,:,:] = ts_pcenter[ini07:fin07,:,:,:]
ts_pedge_msa1[:len_mid,:,:,:] = ts_pedge[ini11:fin11,:,:,:]
ts_pedge_msa1[len_mid:,:,:,:] = ts_pedge[ini07:fin07,:,:,:]
ts_bxh_msa1[:24,:,:,:,:] = ts_bxh[:,:,:,:,:]
ts_bxh_msa1[24:,:,:,:,:] = ts_bxh1[:,:,:,:,:]
##model_hour=24
##len_lvl = 47
##p_edge = np.zeros((len_time,len_lvl+1,len_lat,len_lon))
##p_center = np.zeros((len_time,len_lvl,len_lat,len_lon))
##
##for i in range(len_time):
##   p_edge[i,:,:,:] = geos5_pressure(lat_gc,lon_gc,ts_sp[i,:,:])
##   for j in range(len_lvl):
##      p_center[i,j,:,:] = 0.5*(p_edge[i,j,:,:]+p_edge[i,j+1,:,:])
##print "check"
path = '/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/'
filename = 'ts_pres_msa1'
dataset = Dataset(path+filename+'.nc','w',format='NETCDF4_CLASSIC')
level_edge = dataset.createDimension('level_edge',len_lvl+1)
level_center = dataset.createDimension('level_center',len_lvl)
lat = dataset.createDimension('lat',len_lat)
lon = dataset.createDimension('lon',len_lon)
time = dataset.createDimension('time',len_time)
time_season = dataset.createDimension('time_season',len_time_b)
ts_pc = dataset.createVariable('ts_p_center',np.float32,('time','level_center','lat','lon'))
ts_pc[:,:,:,:] = ts_pcenter_msa1[:,:,:,:]
ts_edge = dataset.createVariable('ts_p_edge',np.float32,('time','level_edge','lat','lon'))
ts_edge[:,:,:,:] = ts_pedge_msa1[:,:,:,:]
ad = dataset.createVariable('ts_air_density',np.float32,('time_season','level_center','lat','lon'))
ad[:,:,:,:] = ts_bxh_msa1[:,0,:,:,:]
bxh = dataset.createVariable('ts_bxheight',np.float32,('time_season','level_center','lat','lon'))
bxh[:,:,:,:] = ts_bxh_msa1[:,1,:,:,:]
dataset.close()
print "done"
