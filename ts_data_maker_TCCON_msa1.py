import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

infile00 = Dataset('/users/jk/15/xzhang/gcadj_std_M_V34/runs/v8-02-01/geos5_mop_0911/ctm.01.20091102.nc')
ts_ap = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_ap_msa1.nc')
ts_mop = Dataset('/users/jk/16/xzhang/TS_annual/geos5_mop_msa1/ts_mop_msa1.nc')
ts_msa05 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_05_msa1/ts_all_05_msa1.nc')
#ts_msa06 = Dataset('/users/jk/15/xzhang/TS_annual/geos5_all_06_2016/ts_all_06_2016.nc')
#strat_ap = Dataset('/users/jk/15/xzhang/TS_annual/geos5_ap_2016/co_strat_ap_2016.nc')
#strat_mop = Dataset('/users/jk/15/xzhang/TS_annual/geos5_mop_2016/co_strat_mop_2016.nc')
#strat_msa05 = Dataset('/users/jk/15/xzhang/TS_annual/geos5_all_05_2016/co_strat_msa05_2016.nc')
#strat_msa06 = Dataset('/users/jk/15/xzhang/TS_annual/geos5_all_06_2016/co_strat_msa06_2016.nc')
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_msa1/ts_pres_msa1.nc')
print "check"

len_time,len_trop,len_lat,len_lon = shape(ts_ap.variables["ts_co"][:,:,:,:])
len_lvl = 47
ts_co = np.zeros((3,len_time,len_lvl,len_lat,len_lon))
ts_co[0,:,:,:,:] = ts_ap.variables["ts_co"][:,:,:,:]
ts_co[1,:,:,:,:] = ts_mop.variables["ts_co"][:,:,:,:]
ts_co[2,:,:,:,:] = ts_msa05.variables["ts_co"][:,:,:,:]
#ts_co[3,:,0:len_trop,:,:] = ts_msa06.variables["ts_co"][:,:,:,:]
#ts_co[0,0:366,len_trop:,:,:] = strat_ap.variables["CO_ap_strat"][:,:,:,:]
#ts_co[1,0:366,len_trop:,:,:] = strat_mop.variables["CO_mop_strat"][:,:,:,:]
#ts_co[2,0:366,len_trop:,:,:] = strat_msa05.variables["CO_msa05_strat"][:,:,:,:]
#ts_co[3,0:366,len_trop:,:,:] = strat_msa06.variables["CO_msa06_strat"][:,:,:,:]
ts_co_col = np.zeros((4,len_time,len_lat,len_lon))
ts_edge = ts_pres.variables["ts_p_edge"][:,0:len_lvl,:,:]
ts_p_center = ts_pres.variables["ts_p_center"][:,:,:,:]
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]

model_hour=24
obs_lvl=71
obs_max=200
mon_counter=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
pwt=np.ndarray(shape=(obs_max,obs_lvl), dtype=np.float64, order='C')
tccon_ak_ts=np.ndarray(shape=(obs_max,obs_lvl), dtype=np.float64, order='C')
tccon_co_ap=np.ndarray(shape=(obs_max,obs_lvl), dtype=np.float64, order='C')
tccon_pressure=np.ndarray(shape=(obs_max,obs_lvl), dtype=np.float64, order='C')
lev=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
tccon_hour=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
tccon_minute=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
XCO=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
tccon_lon=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
tccon_lat=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
station_id=np.ndarray(shape=(obs_max), dtype='|S2', order='C')
tccon_psurf=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
tccon_p_center= np.ndarray(shape=(obs_lvl), dtype=np.float64, order='C')
year = np.array([2009,2010])
month = np.array([11,7])
day0 = np.array([2,2-13])
print "check"
#year = str(2016)
#SPECIFY YEAR, MONTH, AND DAY
#for i in range(12):
   #for j in range(mon_counter[i]):
for month_range in range(2):
   for j in arange(2,15):
      YYYY = str(year[month_range]).zfill(4)
      MM = str(month[month_range]).zfill(2)
      #define iy, i correctly, or you need check the availability of this file using isfile
   #month = str(i+1).zfill(2)
      DD = str(j).zfill(2)
   #day_counter = np.sum(mon_counter[0:i])+j
   #ts_day = 24*day_counter
      ts_day = j-day0[month_range] 
   # READ TCCON FILE
      FILE_NAME='/users/jk/16/dengf/obs/TCCON_2016_CO/tcco_xco_'+YYYY+MM+DD+'.dat'
      #print FILE_NAME
      file_counter=0
      if os.path.isfile(FILE_NAME):
         with open(FILE_NAME) as f:
            while True:
               line=f.readline()
               if ("" == line):               
                  f.close()
                  break;
               list1=line.split()
               tccon_hour[file_counter]=int(list1[2])
               tccon_minute[file_counter]=int(list1[3])
               XCO[file_counter]=float(list1[4])*1.0e+9
               tccon_lon[file_counter]=float(list1[5])
               tccon_lat[file_counter]=float(list1[6])
               tccon_psurf[file_counter]=float(list1[11])
               station_id[file_counter]=list1[12]
               line=f.readline()
               for k in range(71):
                  pwt[file_counter,k]=float(line.split()[k])
               line=f.readline()
               for k in range(71):
                  tccon_ak_ts[file_counter,k]=float(line.split()[k])
               line=f.readline()
               for k in range(71):
                  tccon_co_ap[file_counter,k]=float(line.split()[k])*1.0e+9
               line=f.readline()
               for k in range(71):
                  tccon_pressure[file_counter,k]=float(line.split()[k])
               file_counter=file_counter+1
               # LOCATE TCCON LOCATIONS IN GC MODEL
         #print FILE_NAME, file_counter, station_id
         #block
         for file_i in range(file_counter):
            lat_i,lon_j = get_gc_ij(tccon_lat[file_i],tccon_lon[file_i])
         # LOCATE TCCON TIME IN GC MODEL   
            hour=int(tccon_hour[file_i])
            ts_i = hour+24*ts_day
         # LOCATE TCCON OBS: CO MIXING RATIO AND PRESSURE IN GC MODEL
            ts_co_col[3,ts_i,lat_i,lon_j] = XCO[file_i]
            for exp_i in range(3):
               GC_CO_NATIVE = np.zeros(len_lvl)
               MAP = np.zeros((len_lvl,obs_lvl))
               GC_CO = np.zeros(obs_lvl)
               CO_ap = np.zeros(obs_lvl)
               CO_pert = np.zeros(obs_lvl)
               CO_hat = np.zeros(obs_lvl)
               GC_CO_NATIVE[:]=ts_co[exp_i,ts_i,:,lat_i,lon_j]
            #GC_CO_NATIVE[len_trop:]=ts_co[exp_i,day_counter,len_trop:,lat_i,lon_j]
               
            # LOCATE GC SURFACE PRESSURE  
               gc_psurf =ts_edge[ts_i,0, lat_i, lon_j]
            # MAPPING FROM GC LVLS TO TCCON LVLS, NEED #GC_LVLS,PCENTER ON TCCON GRID, GC_PSURF, #TCCON LVLS, TCCON PRESSURE, TCCON PSURF
               MAP[:, :]=get_intmap(len_lvl,ts_p_center[ts_i,:,lat_i,lon_j], gc_psurf, obs_lvl,tccon_pressure[file_i,:], tccon_psurf[file_i])
            # COMPUTE from GC_CO_NATIVE TO GC_CO  
               for k0 in range(obs_lvl):
                  GC_CO[k0] = 0.0
                  for kk0 in range(len_lvl):
                     GC_CO[k0] = GC_CO[k0] + MAP[kk0,k0]*GC_CO_NATIVE[kk0]            
               for k0 in range (obs_lvl):
                  CO_ap[k0]=tccon_co_ap[file_i, k0]*1.0672
                  CO_pert[k0] = GC_CO[k0] - CO_ap[k0]
                  #print GC_CO[k0], CO_ap[k0]
                  CO_hat[k0] = pwt[file_i,k0]*CO_ap[k0]+pwt[file_i,k0]*tccon_ak_ts[file_i,k0]*CO_pert[k0]
               ts_co_col[exp_i,ts_i,lat_i,lon_j] =  np.sum(CO_hat)#*np.sum(tccon_ak_ts[file_i,:])*0.8/np.sum(tccon_ak_ts[file_i,0:k_trop])
            #print ts_co_col[:,ts_i,lat_i,lon_j]
len_time1 = 24*26#24*np.sum(mon_counter[0:4])
len_time1_ap = 26#np.sum(mon_counter[0:4])
dataset = Dataset('/users/jk/16/xzhang/TCCON/tccon_gc_msa1_feng.nc','w',format='NETCDF4_CLASSIC')
lat_nc = dataset.createDimension('lat',len_lat)
lon_nc = dataset.createDimension('lon',len_lon)
time_nc = dataset.createDimension('time',len_time1)
ts_co_ap_tccon = dataset.createVariable('ts_co_ap_tccon',np.float32,('time','lat','lon'))
ts_co_ap_tccon[:,:,:] = ts_co_col[0,0:len_time1,:,:]
ts_co_mop_tccon = dataset.createVariable('ts_co_mop_tccon',np.float32,('time','lat','lon'))
ts_co_mop_tccon[:,:,:] = ts_co_col[1,0:len_time1,:,:]
ts_co_msa05_tccon = dataset.createVariable('ts_co_msa05_tccon',np.float32,('time','lat','lon'))
ts_co_msa05_tccon[:,:,:] = ts_co_col[2,0:len_time1,:,:]
#ts_co_msa06_tccon = dataset.createVariable('ts_co_msa06_tccon',np.float32,('time','lat','lon'))
#ts_co_msa06_tccon[:,:,:] = ts_co_col[3,0:len_time1,:,:]
ts_co_obs_tccon = dataset.createVariable('ts_co_obs_tccon',np.float32,('time','lat','lon'))
ts_co_obs_tccon[:,:,:] = ts_co_col[3,0:len_time1,:,:]
dataset.close()
print "done"
