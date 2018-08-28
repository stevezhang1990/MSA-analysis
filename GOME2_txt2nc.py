from scipy import *
from pylab import *
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset

line_p = 34
column_p = 6
line2 = 1000000
column1 = 67#16
column2 = 20 #6
column3 = 22 #10
path1 = '/users/jk/16/xzhang/GOME2_NO2/'
year = 2016
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
for month in range(12):
    for day in range(mon_counter[month]):
        YYYY= str(year).zfill(4)
        MM = str(month+1).zfill(2)
        DD = str(day+1).zfill(2)
        i0 = 0
        i1 = 0
        i2 = 0
        i3 = 0
        i = 0
    #open raw data
        filename = path1+'ORBIT_DATA/'+YYYY+'/'+MM+'/'+'no2track'+YYYY+MM+DD+'.txt'
        print "open "+filename
        data0 = np.empty((line_p,column_p),dtype='|S225')
        data1 = np.empty((line2,column1),dtype='|S225')
        data2 = np.empty((line2,column2),dtype='|S225')
        data3 = np.empty((line2,column3),dtype='|S225')

        with open(filename, 'rU') as f:
            reader = csv.reader(f,delimiter =' ',quoting=csv.QUOTE_NONE)
            for line in reader:
            #print line
                i = i + 1
            #str.isalnum()
                if len(line) == column_p: #and str.isdigit(line[0]) == True:
                    for j in range(column_p):
                        data0[i0,j] = line[j]
                    i0 = i0 + 1
                elif len(line) == column1: #and str.isdigit(line[0]) == True:
                    for j in range(column1):
                        data1[i1,j] = line[j]
                    i1 = i1 + 1
                elif len(line) == column2: #and str.isdigit(line[0]) == True:
                    for j in range(column2):
                        data2[i2,j] = line[j]
                    i2 = i2 + 1        
                elif len(line) == column3: #and float(line[0]) != 1: #and str.isdigit(line[0]) == True:
                    for j in range(column3):
                        data3[i3,j] = line[j]
                    i3 = i3 + 1

        data1 = data1[0:i1,:]
        data2 = data2[0:i2,:]
        data3 = data3[1:i3,:]
        AB = np.zeros((line_p,2))
        data_raw = np.zeros((i1,101))
        for i in range(line_p):
            AB[i,0]= float(data0[i,0])
            AB[i,1]= float(data0[i,2])
#date,time,lon,lat,vcd,sigvcd,vcdtrop,sigvcdt,vcdstrat,sgvcds,flagtrop,psurf,sigvcdak (VCDerr for AK)
#sigvcdtak (TVCDerr for AK), AK (43 levels)
# sza, vza, relative azimuth angle, subset counter, loncorn, latcorn
# scd, amf, amftrop, amftrop, amfgeo, scdstr,clfrac,cltpress,albclr,crfrac,level_tropo,ghostcol
        for i in range(i1):
            for counter1 in range(12):
                data_raw[i,counter1] = float(data1[i,2*counter1])
        #data_raw[i,2] = (float(data1[i,4])+180)%360-180
            for j in range(line_p):
                data_raw[i,j+12] = AB[j,0]+data_raw[i,11]*AB[j,1]
                data_raw[i,j+46] = float(data1[i,28+j])
            for counter2 in range(3):
                data_raw[i,80+counter2] = float(data2[i,2*counter2])
            for k in range(4):
                data_raw[i,83+k] = float(data2[i,8+k])
                data_raw[i,87+k] = float(data2[i,13+k])
            for counter3 in range(10):
                data_raw[i,91+counter3] = float(data3[i,2*counter3])

        dataset = Dataset(path1+YYYY+'/'+MM+'/'+'GOME2_NO2_metopa_'+YYYY+MM+DD+'_v2.3.nc','w',format='NETCDF4_CLASSIC')
        pcenter = dataset.createDimension('pcenter',line_p)
        pedge = dataset.createDimension('pedge',line_p+1)
        n_orbit = dataset.createDimension('n_orbit',i1)
        corner_dim = dataset.createDimension('corner_dim',4)
        date = dataset.createVariable('date_orbit',np.int32,('n_orbit'))
        date[:] = data_raw[:,0].astype(int)
        time = dataset.createVariable('time_orbit',np.float32,('n_orbit'))
        time[:] = data_raw[:,1]
        lon = dataset.createVariable('lon_orbit',np.float32,('n_orbit'))
        lon[:] = data_raw[:,2] #(data_raw[:,2]+180)%360-180
        lat = dataset.createVariable('lat_orbit',np.float32,('n_orbit'))
        lat[:] = data_raw[:,3]
        vcd = dataset.createVariable('total_VCD_orbit',np.float32,('n_orbit'))
        vcd[:] = data_raw[:,4]
        sigvcd = dataset.createVariable('total_VCD_error_orbit',np.float32,('n_orbit'))
        sigvcd[:] = data_raw[:,5]
        vcdt = dataset.createVariable('trop_VCD_orbit',np.float32,('n_orbit'))
        vcdt[:] = data_raw[:,6]
        sigvcdt = dataset.createVariable('trop_VCD_error_orbit',np.float32,('n_orbit'))
        sigvcdt[:] = data_raw[:,7]
        vcdstrat = dataset.createVariable('strat_VCD_orbit',np.float32,('n_orbit'))
        vcdstrat[:] = data_raw[:,8]
        sigvcds = dataset.createVariable('strat_VCD_error_orbit',np.float32,('n_orbit'))
        sigvcds[:] = data_raw[:,9]
        flagtrop = dataset.createVariable('trop_flag_orbit',np.int32,('n_orbit'))
        flagtrop[:] = data_raw[:,10].astype(int)
        pres = dataset.createVariable('pressure_levels_orbit',np.float32,('n_orbit','pedge'))
        pres[:,:] = data_raw[:,11:46]
        avk = dataset.createVariable('averaging_kernels_orbit',np.float32,('n_orbit','pcenter'))
        avk[:,:] = data_raw[:,46:80]
        sza = dataset.createVariable('solar_zenith_angle_orbit',np.float32,('n_orbit'))
        sza[:] = data_raw[:,80]
        vza = dataset.createVariable('viewing_zenith_angle_orbit',np.float32,('n_orbit'))
        vza[:] = data_raw[:,81]
        raa = dataset.createVariable('relative_azimuth_angle_orbit',np.float32,('n_orbit'))
        raa[:] = data_raw[:,82]
        lon_corn = dataset.createVariable('corner_lon_orbit',np.float32,('n_orbit','corner_dim'))
        lon_corn[:,:] = data_raw[:,83:87]
        lat_corn = dataset.createVariable('corner_lat_orbit',np.float32,('n_orbit','corner_dim'))
        lat_corn[:,:] = data_raw[:,87:91]
        scd = dataset.createVariable('total_SCD_orbit',np.float32,('n_orbit'))
        scd[:] = data_raw[:,91]
        amf = dataset.createVariable('total_AMF_orbit',np.float32,('n_orbit'))
        amf[:] = data_raw[:,92]
        amftrop = dataset.createVariable('trop_AMF_orbit',np.float32,('n_orbit'))
        amftrop[:] = data_raw[:,93]
        amfgeo = dataset.createVariable('geometrical_AMF_orbit',np.float32,('n_orbit'))
        amfgeo[:] = data_raw[:,94]
        scdstrat = dataset.createVariable('strat_SCD_orbit',np.float32,('n_orbit'))
        scdstrat[:] = data_raw[:,95]
        cloud_frac = dataset.createVariable('cloud_fraction_orbit',np.float32,('n_orbit'))
        cloud_frac[:] = data_raw[:,96]
        cloud_tp =dataset.createVariable('cloud_top_pres_orbit',np.float32,('n_orbit'))
        cloud_tp[:] = data_raw[:,97]
        albedo_frac =dataset.createVariable('albedo_fraction_orbit',np.float32,('n_orbit'))
        albedo_frac[:] = data_raw[:,98]
        crfrac = dataset.createVariable('cloud_radiance_frac_orbit',np.int32,('n_orbit'))
        crfrac[:] = data_raw[:,99].astype(int)
        ltropo = dataset.createVariable('tropopause_level',np.int32,('n_orbit'))
        ltropo[:] = data_raw[:,100].astype(int)
        dataset.close()    

print "done"
