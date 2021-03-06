#!/usr/bin/python

from scipy import *
from pylab import *
import csv
import numpy as np
from get_intmap import *
import icartt
import numpy.ma as ma
from netCDF4 import Dataset

path = '/users/jk/16/xzhang/ARCTAS/'
date_num = ['01','04','05','08','09','10','13']
expdir = ['ap_0807_NA','omi_0807_NA','omo_0807_NA','omo_0807_BC']
year = 2008
month = 7
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080701.nc')
lat_gc = infile00.variables["LAT"][:]
lon_gc = infile00.variables["LON"][:]
len_time = 24*27
len_lvl = 47
len_lat = len(lat_gc)
len_lon = len(lon_gc)
ts_gc = np.zeros((5,len_time,len_lvl,len_lat,len_lon))
ts_arctas = np.zeros((6,len_time,len_lvl,len_lat,len_lon))
for j in range(4):
    filename = path+'ts_'+expdir[j]+'.npy'
    ts_gc[j,:,:,:,:] = 1e3*np.load(filename)
for file_i in range(len(date_num)):
    data_MMS = icartt.Dataset(path+'ARCTAS-mrg60-dc8_merge_200807'+date_num[file_i]+'_R14.ict')
    time_MMS = np.array(data_MMS['UTC'],dtype=float)
    lat_MMS = np.array(data_MMS['LATITUDE'],dtype=float)
    lon_MMS = np.array(data_MMS['LONGITUDE'],dtype=float)-360
    alt_MMS = np.array(data_MMS['GPS_Altitude'],dtype=float)
    NO2_NCAR_MMS = np.array(data_MMS['NO2_NCAR'],dtype=float)
    NO2_UCB_MMS = np.array(data_MMS['NO2_UCB'],dtype=float)
    len2 = 13
    data_med = np.zeros((len(time_MMS),len2))
    hms_m,hms_s = divmod(time_MMS,60)
    hms_h,hms_m = divmod(hms_m,60)
    #with open(path+'NOxyO3_DC8_200807'+date_num[file_i]+'_R3.ict', 'rU') as f:
        #reader = csv.reader(f,delimiter =' ',quoting=csv.QUOTE_NONE)
    #data_NO2 = icartt.Dataset(path+'NOxyO3_DC8_200807'+date_num[file_i]+'_R3.ict')
        #line_i = 0
        #for line in reader:
            #if line_i == 0:
                #line_invalid = int(line[0])
            #elif line_i >= line_invalid:
                #line_data = [c for c in line if c !='']
                #NO2_MMS[line_i-line_invalid,0] = float(line_data[0])
                #NO2_MMS[line_i-line_invalid,1] = float(line_data[5])
            #line_i = line_i + 1
# to know the variable names: type data.varnames
    #time1_NO2 = NO2_MMS[:line_i-line_invalid,0]
    #NO2_CL = NO2_MMS[:line_i-line_invalid,1]
#---0:year---1:month---2:day---3:time---4:lon---5:lat---6:alt(m)--
#---7:NO2
    for i in range(len(time_MMS)):
        if NO2_NCAR_MMS[i] > 0 or NO2_NCAR_MMS[i] > 0:
            data_med[i,7] = np.ma.masked_where(NO2_NCAR_MMS[i]<=0,NO2_NCAR_MMS[i])
            data_med[i,8] = np.ma.masked_where(NO2_UCB_MMS[i]<=0,NO2_UCB_MMS[i])
            data_med[i,0] = int(year)
            data_med[i,1] = int(month)
            hms = 100*hms_h[i]+hms_m[i]+hms_s[i]/60.0
            if hms > 2400:
                data_med[i,2] = int(date_num[file_i])+1
                data_med[i,3] = np.around(hms-2400,decimals=1)
            else:
                data_med[i,2] = int(date_num[file_i])
                data_med[i,3] = np.around(hms,decimals=1)
            data_med[i,4] = np.around(lon_MMS[i],decimals=2)
            data_med[i,5] = np.around(lat_MMS[i],decimals=2)
            data_med[i,6] = np.around(alt_MMS[i],decimals=2)
            ts_i = int(24*(data_med[i,2]-1)+data_med[i,3]/100)
            lat_i,lon_j = get_gc_ij_NA(data_med[i,5],data_med[i,4])
            #print data_med[i,5],data_med[i,4],lat_gc[lat_i],lon_gc[lon_j]
            lvl_k = get_gc_lvl(data_med[i,6])
            if lat_i < len_lat and lon_j< len_lon and lvl_k <47:
                ts_arctas[4,ts_i,lvl_k,lat_i,lon_j] = np.mean(data_med[i,7:9])
                ts_arctas[5,ts_i,lvl_k,lat_i,lon_j] = data_med[i,6]
                for j in range(4):               
                    #data_med[i,9+j] = ts_arctas[j,ts_i,lvl_k,lat_i,lon_j]
                    ts_arctas[j,ts_i,lvl_k,lat_i,lon_j] = ts_gc[j,ts_i,lvl_k,lat_i,lon_j]
    #data_filt = ma.masked_invalid(data_med)
    #data_avg = np.zeros((len3,len2))
    #time_avg = np.zeros((len3,len2))
    #data_pre = np.zeros((len3,len2))
    #data_fin = np.zeros(len2)
    #time_fin = 0
    #if file_i == 0:
        #data_med2 = data_med
    #else:
        #data_med2 = np.vstack((data_med2,data_med))
ts_arctas_filt = ma.masked_invalid(np.ma.masked_where(ts_arctas<=0,ts_arctas))
ts_arctas_2d = np.mean(np.mean(ts_arctas_filt,axis=1),axis=1)
#dataset = Dataset(path+'arctas_no2_200807.nc','w',format='NETCDF4_CLASSIC')
#data_orbit = dataset.createDimension('data_orbit',len(data_med2[:,0]))
#year_orbit = dataset.createVariable('year_orbit',np.float32,('data_orbit'))
#year_orbit[:] = data_med2[:,0]
#month_orbit = dataset.createVariable('month_orbit',np.float32,('data_orbit'))
#month_orbit[:] = data_med2[:,1]
#lat_orbit = dataset.createVariable('lat_orbit',np.float32,('data_orbit'))
#lat_orbit[:] = data_med2[:,5]
#lon_orbit = dataset.createVariable('lon_orbit',np.float32,('data_orbit'))
#lon_orbit[:] = data_med2[:,4]
#alt_orbit = dataset.createVariable('alt_orbit',np.float32,('data_orbit'))
#alt_orbit[:] = data_med2[:,6]
#day_orbit = dataset.createVariable('day_orbit',np.float32,('data_orbit'))
#day_orbit[:] = data_med2[:,2]
#hms_orbit = dataset.createVariable('hms_orbit',np.float32,('data_orbit'))
#hms_orbit[:] = data_med2[:,3]
#NO2_NCAR_orbit = dataset.createVariable('NO2_NCAR_orbit',np.float32,('data_orbit'))
#NO2_NCAR_orbit[:] = data_med2[:,7]
#NO2_UCB_orbit = dataset.createVariable('NO2_UCB_orbit',np.float32,('data_orbit'))
#NO2_UCB_orbit[:] = data_med2[:,8]
#NO2_ap_NA = dataset.createVariable('NO2_ap_NA',np.float32,('data_orbit'))
#NO2_ap_NA[:] = data_med2[:,9]
#NO2_omi_NA = dataset.createVariable('NO2_omi_NA',np.float32,('data_orbit'))
#NO2_omi_NA[:] = data_med2[:,10]
#NO2_omo_NA = dataset.createVariable('NO2_omo_NA',np.float32,('data_orbit'))
#NO2_omo_NA[:] = data_med2[:,11]
#NO2_omo_BC = dataset.createVariable('NO2_omo_BC',np.float32,('data_orbit'))
#NO2_omo_BC[:] = data_med2[:,12]
dataset = Dataset(path+'arctas_no2_2d_200807.nc','w',format='NETCDF4_CLASSIC')
len_exp = dataset.createDimension('len_exp',6) 
len_lat = dataset.createDimension('len_lat',len_lat)
len_lon = dataset.createDimension('len_lon',len_lon)
ts_2d = dataset.createVariable('ts_arctas_2d',np.float32,('len_exp','len_lat','len_lon'))
ts_2d[:,:,:] = ts_arctas_2d[:,:,:]
dataset.close()
print "done"
