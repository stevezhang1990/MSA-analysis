#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
from get_intmap import *
    
file_count = 0
ap_corr = 1.0672
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
ts_sp = np.load('/users/jk/15/xzhang/TS_annual/geos5_ap_2016/ts_ps_2016.npy')
ts_ap = Dataset('/users/jk/15/xzhang/TS_annual/geos5_ap_2016/ts_ap_2016.nc')
ts_pres = Dataset('/users/jk/15/xzhang/TS_annual/geos5_ap_2016/ts_pres_2016.nc')
ts_mop = Dataset('/users/jk/15/xzhang/TS_annual/geos5_mop_2016/ts_mop_2016.nc')
ts_msa05 = Dataset('/users/jk/15/xzhang/TS_annual/geos5_all_05_2016/ts_all_05_2016.nc')
ts_msa06 = Dataset('/users/jk/15/xzhang/TS_annual/geos5_all_06_2016/ts_all_06_2016.nc')
len_lat, len_lon = shape(infile00.variables["DXYP__DXYP"][0,:,:])
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
filename = '/users/jk/15/xzhang/TCCON/tccon_2016.nc'
tccon_data = Dataset(filename)
lat_obs = tccon_data.variables["station_lat"][:]
lon_obs = tccon_data.variables["station_lon"][:]
sza = tccon_data.variables["sta_sza"][:,:]
co_col_obs = tccon_data.variables["sta_xco"][:,:]
co_ak = tccon_data.variables["sta_ak_co"][:,:,:]
block
sza_ak = tccon_data.variables["sta_ak_zenith"][:,:]
prior_height = tccon_data.variables["sta_prior_height"][:,:]
p_ak = tccon_data.variables["sta_ak_pres"][:,:]
co_ap = tccon_data.variables["sta_prior_co"][:,:,:]
pres_ap = tccon_data.variables["sta_prior_pres"][:,:,:]
pwt = tccon_data.variables["prior_pwt"][:]
len_site,len_time = shape(sza)
len_akz = len(sza_ak[0,:])
len_time_ap, len_p = shape(co_ap[0,:,:])
len_lvl = 47
len_time, len_trop, len_lat, len_lon = shape(ts_ap.variables["ts_co"][:,:,:,:])
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
len_time1 = 24*np.sum(mon_counter[0:4])
len_time1_ap = np.sum(mon_counter[0:4])
CO_gc_annual = np.zeros((4,len_time,len_lvl,len_lat,len_lon))
CO_gc_annual[0,:,0:len_trop,:,:] = ts_ap.variables["ts_co"][:,:,:,:]
CO_gc_annual[1,:,0:len_trop,:,:] = ts_mop.variables["ts_co"][:,:,:,:]
CO_gc_annual[2,:,0:len_trop,:,:] = ts_msa05.variables["ts_co"][:,:,:,:]
CO_gc_annual[3,:,0:len_trop,:,:] = ts_msa06.variables["ts_co"][:,:,:,:]
ts_co_col = np.zeros((5,len_time,len_lat,len_lon))
gc_p_edge = ts_pres.variables["ts_p_edge"][:,0:len_lvl,:,:]
gc_p_center = ts_pres.variables["ts_p_center"][:,:,:,:]
print "check"

for i in range(len_site):
    lat_i,lon_j = get_gc_ij(lat_obs[i],lon_obs[i])
    for j in range(len_time1):
        if co_col_obs[i,j] > 0:
            jj = int(j/24)
            ts_co_col[4,j,lat_i,lon_j] = co_col_obs[i,j]
            n = int(round((sza[i,j]-sza_ak[i,0])/(sza_ak[i,1]-sza_ak[i,0])))
            #print sza_ak[i,n]
            #print sza[i,j]
            mod_p0 = gc_p_edge[j,0,lat_i,lon_j]
            mod_pc = gc_p_center[j,:,lat_i,lon_j]
            obs_pc = pres_ap[i,jj,:]
            obs_p0 = pres_ap[i,jj,0]
            block
            for exp_i in range(4):
                gc_co = np.zeros(len_p)
                co_apc = np.zeros(len_p)
                co_pert = np.zeros(len_p)
                co_hat = np.zeros(len_p)
                gc_co_native = CO_gc_annual[exp_i,j,:,lat_i,lon_j]
                MAP = get_intmap(len_lvl,mod_pc,mod_p0,len_p,obs_pc,obs_p0)
                for k in range(len_p):
                    gc_co[k] = 0
                    for kk in range(len_lvl):
                        gc_co[k] = gc_co[k] + MAP[kk,k]*gc_co_native[kk]  
                for k in range(len_p):
                    co_apc[k] = ap_corr*co_ap[i,jj,k]
                    co_pert[k] = gc_co[k] - co_apc[k]
                    co_hat[k] =  pwt[k]*co_apc[k] + pwt[k]*co_ak[i,k,n]*co_pert[k]
                ts_co_col[exp_i,j,lat_i,lon_j] = np.sum(co_hat)
                #print ts_co_col[:,j,lat_i,lon_j]
dataset = Dataset('/users/jk/15/xzhang/TCCON/tccon_gc_col1.nc','w',format='NETCDF4_CLASSIC')
lat_nc = dataset.createDimension('lat',len_lat)
lon_nc = dataset.createDimension('lon',len_lon)
time_nc = dataset.createDimension('time',len_time1)
ts_co_ap_tccon = dataset.createVariable('ts_co_ap_tccon',np.float32,('time','lat','lon'))
ts_co_ap_tccon[:,:,:] = ts_co_col[0,0:len_time1,:,:]
ts_co_mop_tccon = dataset.createVariable('ts_co_mop_tccon',np.float32,('time','lat','lon'))
ts_co_mop_tccon[:,:,:] = ts_co_col[1,0:len_time1,:,:]
ts_co_msa05_tccon = dataset.createVariable('ts_co_msa05_tccon',np.float32,('time','lat','lon'))
ts_co_msa05_tccon[:,:,:] = ts_co_col[2,0:len_time1,:,:]
ts_co_msa06_tccon = dataset.createVariable('ts_co_msa06_tccon',np.float32,('time','lat','lon'))
ts_co_msa06_tccon[:,:,:] = ts_co_col[3,0:len_time1,:,:]
ts_co_obs_tccon = dataset.createVariable('ts_co_obs_tccon',np.float32,('time','lat','lon'))
ts_co_obs_tccon[:,:,:] = ts_co_col[4,0:len_time1,:,:]
dataset.close()
print "done"
show()
