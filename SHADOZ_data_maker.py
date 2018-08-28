#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
import csv
from get_intmap import *

path = '/users/jk/16/xzhang/SHADOZ/'
year = 2016
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
ts_ap = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_ap_2016.nc')
ts_mop = Dataset('/users/jk/16/xzhang/TS_annual/geos5_iasio3_2016/ts_iasio3_2016.nc')
ts_msa07 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/ts_all_07_2016.nc')
ts_msa08 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_08_2016/ts_all_08_2016.nc')
strat_ap = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/o3_strat_ap_2016.nc')
strat_mop = Dataset('/users/jk/16/xzhang/TS_annual/geos5_iasio3_2016/o3_strat_iasio3_2016.nc')
strat_msa07 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_07_2016/o3_strat_msa07_2016.nc')
strat_msa08 = Dataset('/users/jk/16/xzhang/TS_annual/geos5_all_08_2016/o3_strat_msa08_2016.nc')
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
len_time,len_trop,len_lat,len_lon = shape(ts_ap.variables["ts_o3"][:,:,:,:])
len_lvl = 47
#len_lvl,len_lat,len_lon = shape(infile00.variables["IJ-AVG-S__Ox"][:,:,:])
ts_o3 = np.zeros((4,len_time,len_lvl,len_lat,len_lon))
ts_o3[0,:,0:len_trop,:,:] = ts_ap.variables["ts_o3"][:,:,:,:]
ts_o3[1,:,0:len_trop,:,:] = ts_mop.variables["ts_o3"][:,:,:,:]
ts_o3[2,:,0:len_trop,:,:] = ts_msa07.variables["ts_o3"][:,:,:,:]
ts_o3[3,:,0:len_trop,:,:] = ts_msa08.variables["ts_o3"][:,:,:,:]
ts_o3[0,0:366,len_trop:,:,:] = strat_ap.variables["O3_ap_strat"][:,:,:,:]
ts_o3[1,0:366,len_trop:,:,:] = strat_mop.variables["O3_mop_strat"][:,:,:,:]
ts_o3[2,0:366,len_trop:,:,:] = strat_msa07.variables["O3_msa05_strat"][:,:,:,:]
ts_o3[3,0:366,len_trop:,:,:] = strat_msa08.variables["O3_msa06_strat"][:,:,:,:]

len_dtime = 366*24
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
data_raw = np.zeros((len_dtime,len_lvl,len_lat,len_lon,6))

for filename in glob.glob(os.path.join(path, '*.dat')):
    with open(filename, 'rU') as f:
        print filename
        line_i = 0
        reader = csv.reader(f,delimiter =' ',quoting=csv.QUOTE_NONE)
        for line in reader:
            if line_i == 0:
                line_invalid = int(line[0])
            elif line_i == 7:
                for col_j in range(len(line)):
                    if line[col_j] == ':':
                        lat = float(line[col_j+1])
                        #print lat
            elif line_i == 8:
                for col_j in range(len(line)):
                    if line[col_j] == ':':
                        lon = float(line[col_j+1])
                        #print lon
                        lat_i,lon_j = get_gc_ij(lat,lon)
                        #print lon_gc[lon_j]
            elif line_i == 10:
                for col_j in range(len(line)):
                    if line[col_j] == ':':
                        date = int(line[col_j+1])
            elif line_i == 11:
                for col_j in range(len(line)):
                    if line[col_j] == ':':
                        line1 = line[col_j+1]
                        hh,mm,ss = line1.split(":")
                        hour = int(hh)
            elif line_i >= line_invalid:
                line_data = [c for c in line if c !='']
                #print line_data
                time_tau = float(line_data[0])
                pres = float(line_data[1])                    
                alt = float(line_data[2])
                O3 = float(line_data[6])*1e3
                #print O3
                #print time_tau,pres,alt,O3
                month = int((date-2016*1e4)/1e2)
                day = date-100*month - 2016*1e4
                day_i = int(np.sum(mon_counter[0:(month-1)])+day-1)
                time_i = int(24*day_i+hour)
                #print day, time_i                
                lvl_l = get_gc_lvl(alt)
                if O3 > 0 and O3 < 1e5:
                    data_raw[time_i,lvl_l,lat_i,lon_j,4] = data_raw[time_i,lvl_l,lat_i,lon_j,4] + 1
                    data_raw[time_i,lvl_l,lat_i,lon_j,5] = data_raw[time_i,lvl_l,lat_i,lon_j,5] + O3
                    for k in range(4):
                        if lvl_l >= len_trop:
                            data_raw[time_i,lvl_l,lat_i,lon_j,k] = ts_o3[k,day_i,lvl_l,lat_i,lon_j]
                        else:
                            data_raw[time_i,lvl_l,lat_i,lon_j,k] = ts_o3[k,time_i,lvl_l,lat_i,lon_j]
                    #print data_raw[time_i,lvl_l,lat_i,lon_j,:]
            line_i = line_i + 1
dataset = Dataset(path+'shadoz_obs_gc_2016.nc','w',format='NETCDF4_CLASSIC')
len_time = dataset.createDimension('len_time',len_dtime)
len_level = dataset.createDimension('len_lvl',len_lvl)
len_lat = dataset.createDimension('len_lat',len_lat)
len_lon = dataset.createDimension('len_lon',len_lon)
O3_shadoz = dataset.createVariable('shadoz_o3_obs',np.float32,('len_time','len_lvl','len_lat','len_lon'))
O3_shadoz[:,:,:,:] = data_raw[:,:,:,:,5]/data_raw[:,:,:,:,4]
O3_ap_sha = dataset.createVariable('shadoz_o3_ap',np.float32,('len_time','len_lvl','len_lat','len_lon'))
O3_ap_sha[:,:,:,:] = data_raw[:,:,:,:,0]
O3_iasio3_sha = dataset.createVariable('shadoz_o3_iasio3',np.float32,('len_time','len_lvl','len_lat','len_lon'))
O3_iasio3_sha[:,:,:,:] = data_raw[:,:,:,:,1]
O3_msa07_sha = dataset.createVariable('shadoz_o3_msa07',np.float32,('len_time','len_lvl','len_lat','len_lon'))
O3_msa07_sha[:,:,:,:] = data_raw[:,:,:,:,2]
O3_msa08_sha = dataset.createVariable('shadoz_o3_msa08',np.float32,('len_time','len_lvl','len_lat','len_lon'))
O3_msa08_sha[:,:,:,:] = data_raw[:,:,:,:,3]

dataset.close()
print "done"
