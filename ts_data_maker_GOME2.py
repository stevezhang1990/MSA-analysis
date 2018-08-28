#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
from get_intmap import *
    
file_count = 0
#xnumolair = 6.022d23/(28.964d-3)
infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_M_V34/runs/v8-02-01/geos5_mop_0911/ctm.01.20091102.nc')
ts_ap = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_ap_2016.nc')
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_pres_msa1.nc')
ts_pres1 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_pres_2016.nc')
ts_omi = Dataset('/users/jk/16/xzhang/TS_annual/geos5_omino2_2016/ts_omino2_2016.nc')
ts_msa07 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/ts_all_07_2016.nc')
gome2_data = Dataset('/users/jk/16/xzhang/GOME2_NO2/gome2_msa2.nc')
ts_msa08 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_08_2016/ts_all_08_2016.nc')
len_time, len_lvl, len_lat, len_lon = shape(ts_ap.variables["ts_co"][:,:,:,:])
#len_lat, len_lon = shape(infile00.variables["DXYP__DXYP"][0,:,:])
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
amf_total_gc = gome2_data.variables["amf_total_gc"][:,:,:]
amf_trop_gc = gome2_data.variables["amf_trop_gc"][:,:,:]
trop_level_gc = gome2_data.variables["trop_level_gc"][:,:,:]
ak_total_gc = gome2_data.variables["ak_total_gc"][:,:,:,:]
pres_gc = gome2_data.variables["pres_gc"][:,:,:,:]
pres_pc = 0.5*(pres_gc[:,:,:,:-1]+pres_gc[:,:,:,1:])
len_p = len(pres_pc[0,0,0,:])
ts_no2_col = np.zeros((5,len_time,len_lat,len_lon))
vcdt_gc = gome2_data.variables["vcdt_gc"][:,:,:]
vcd_gc = gome2_data.variables["vcd_gc"][:,:,:]
#pcenter_gc = ts_pres1.variables["ts_p_center"][:,:,:,:]
#pedge_gc = ts_pres1.variables["ts_p_edge"][:,:,:,:]
bxh_gc = ts_pres.variables["ts_bxheight"][:,:,:,:]
airden_gc = ts_pres.variables["ts_air_density"][:,:,:,:]

#mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
#len_time1 = 24*np.sum(mon_counter[0:4])
#len_time1_ap = np.sum(mon_counter[0:4])
NO2_gc_annual = np.zeros((4,len_time,len_lvl,len_lat,len_lon))
NO2_gc_annual[0,:,:,:,:] = ts_ap.variables["ts_no2"][:,:,:,:]
NO2_gc_annual[1,:,:,:,:] = ts_omi.variables["ts_no2"][:,:,:,:]
NO2_gc_annual[2,:,:,:,:] = ts_msa07.variables["ts_no2"][:,:,:,:]
NO2_gc_annual[3,:,:,:,:] = ts_msa08.variables["ts_no2"][:,:,:,:]


gc_p_edge = ts_pres1.variables["ts_p_edge"][:,0:len_lvl+1,:,:]
gc_p_center = ts_pres1.variables["ts_p_center"][:,:,:,:]

print "check"

for i in range(len_time):
    for j in range(len_lat):
        for k in range(len_lon):
            if vcdt_gc[i,j,k] > 0:
                ts_no2_col[4,i,j,k] = vcdt_gc[i,j,k]*1e15
                #print vcdt_gc[i,j,k]*1e15
                mod_p0 = gc_p_edge[i,0,j,k]
                mod_pedge = gc_p_edge[i,:,j,k]
                obs_pedge = pres_gc[i,j,k,:]
                obs_p0 = pres_gc[i,j,k,0]
                crt_lvl = trop_level_gc[i,j,k]
                for exp_i in range(4):
                    gc_no2 = np.zeros(len_p)
                    no2_hat = np.zeros(len_p)
                    ak_trop = np.zeros(len_p)
                    gc_no2_native = NO2_gc_annual[exp_i,i,:,j,k]*1e-9*airden_gc[i%24,:38,j,k]*bxh_gc[i%24,:38,j,k]*100
                    #print np.sum(gc_no2_native[:38])
                    MAP = get_intmap(len_lvl,mod_pedge,mod_p0,len_p,obs_pedge,obs_p0,mode='PC')
                    for h in range(len_p):
                        gc_no2[h] = 0
                        for hh in range(len_lvl):
                            gc_no2[h] = gc_no2[h] + MAP[hh,h]*gc_no2_native[hh]
                #for h in arange(1,len_p-1):
                    #if gc_no2[h] == 0:
                        #if gc_no2[h-1] > 0:
                            #lvl_crt1 = h-1
                        #if gc_no2[h+1] > 0:
                            #lvl_crt2 = h+2
                        #if lvl_crt1*lvl_crt2>0:
                            #for hh in arange(lvl_crt1,lvl_crt2):
                                #gc_no2[hh] = (gc_no2[lvl_crt1]+gc_no2[lvl_crt2])/np.sum()
                    #print amf_trop_gc[i,j,k],ak_total_gc[i,j,k,:],amf_total_gc[i,j,k]
                    for h in range(len_p):
                        ak_trop[h] = ak_total_gc[i,j,k,h]*amf_total_gc[i,j,k]/amf_trop_gc[i,j,k]
                        no2_hat[h] = ak_trop[h]*gc_no2[h]
                    ts_no2_col[exp_i,i,j,k] = np.sum(no2_hat[:crt_lvl])
                #print ts_no2_col[:,i,j,k]

dataset = Dataset('/users/jk/16/xzhang/GOME2_NO2/gome2_gc_msa2.nc','w',format='NETCDF4_CLASSIC')
lat_nc = dataset.createDimension('lat',len_lat)
lon_nc = dataset.createDimension('lon',len_lon)
time_nc = dataset.createDimension('time',len_time)
ts_no2_ap_tccon = dataset.createVariable('ts_no2_ap_gome2',np.float32,('time','lat','lon'))
ts_no2_ap_tccon[:,:,:] = ts_no2_col[0,:,:,:]
ts_no2_omi_tccon = dataset.createVariable('ts_no2_omi_gome2',np.float32,('time','lat','lon'))
ts_no2_omi_tccon[:,:,:] = ts_no2_col[1,:,:,:]
ts_no2_msa07_tccon = dataset.createVariable('ts_no2_msa07_gome2',np.float32,('time','lat','lon'))
ts_no2_msa07_tccon[:,:,:] = ts_no2_col[2,:,:,:]
ts_no2_msa08_tccon = dataset.createVariable('ts_no2_msa08_gome2',np.float32,('time','lat','lon'))
ts_no2_msa08_tccon[:,:,:] = ts_no2_col[3,:,:,:]
ts_no2_obs_tccon = dataset.createVariable('ts_no2_obs_gome2',np.float32,('time','lat','lon'))
ts_no2_obs_tccon[:,:,:] = ts_no2_col[4,:,:,:]
dataset.close()
print "done"
show()
