from scipy import *
from pylab import *
import csv, os, glob, h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import h5py
def tai2utc(tai93):
    leapsec = np.array([ 15638399, 47174400, 94608001, 141868802, 189302403, 410227204, 504921605, 615254406, 709862407, 757382408 ])
    min_base = 60
    hour_base = 60*min_base
    day_base = 24*hour_base
    year_base = 365*day_base
    mdays = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    tai = tai93
    sec =0
    for i in range(len(leapsec)):        
        if int(tai93) == leapsec[i]+1:
            sec = 60 + tai93- int(tai93)
        if int(tai93) > leapsec[i]:
            tai = tai -1
    year = 1993 + int(tai/year_base)
    wk = tai-(year-1993)*year_base-int((year-1993)/4)*day_base
    if wk<0:
        year = year-1
        wk = tai - (year-1993)*year_base-int((year-1993)/4)*day_base
    days = int(wk/day_base)
    wk = wk - days*day_base
    month = 1
    for i in range(len(mdays)):
        leap = 0
        if month ==2 and year%4==0:
            leap = 1
        if month ==i+1 and days>=mdays[i]+leap:
            month = month+1
            days = days - mdays[i] - leap
    day = days+1
    hour = int(wk/hour_base)
    wk = wk-hour*hour_base
    imin = int(wk/min_base)
    if sec<60:
        sec = wk - imin*min_base
    return year,month,day,hour,imin,sec
    
    

path1 = '/users/jk/16/xzhang/TES_O3/ORBIT_DATA/2006/'
year = 2006
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
dname = ['' for x in range(11)]

#YYYY= str(year).zfill(4)
#MM = str(month+1).zfill(2)
#DD = str(day+1).zfill(2)
for filename in glob.glob(os.path.join(path1, '*.he5')):
    print filename
    infile = h5py.File(filename, mode='r')
    dname[0]= '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/Altitude'
    dname[1] = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/AveragingKernel'
    dname[2] = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/O3'
    dname[3] = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/ObservationErrorCovariance'
    dname[4] = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/Pressure'
    dname[5] = '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/SpeciesRetrievalQuality'
    dname[7] =  '/HDFEOS/SWATHS/O3NadirSwath/Data Fields/ConstraintVector'
    dname[8] =  '/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields/Latitude'
    dname[9] =  '/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields/Longitude'
    dname[10] =  '/HDFEOS/SWATHS/O3NadirSwath/Geolocation Fields/Time'
    len1,len2,len3 = shape(infile[dname[1]][:,:,:])
    ymd = np.zeros(len1)
    i_crt = np.array([0,len1,len1,len1,len1,len1],dtype=int)
    j=0
    for i in range(len1):
        iyear,imonth,iday,ihour,imin,isec = tai2utc(infile[dname[10]][i])
        if i==0:
            day0 = iday
        elif iday != day0:
            day0 = iday
            j=j+1
            i_crt[j] = int(i)
        ymd[i] = 1e4*iyear+imonth*100+iday+(ihour*3600+imin*60+isec)/(24*3600.0)
    for i in range(4):
        len_obs = len(ymd[i_crt[i]:i_crt[i+1]])
        if len_obs!=0:
            year = int(ymd[i_crt[i]]/1e4)
            YYYY = str(year).zfill(4)
            month = int((ymd[i_crt[i]]-year*1e4)/100)
            MM = str(month).zfill(2)
            day = int(ymd[i_crt[i]]-year*1e4-month*100)
            DD = str(day).zfill(2)
            dataset = Dataset(path1+'tes_aura_nadir_'+YYYY+MM+DD+'_O3_v6.nc','w',format='NETCDF4_CLASSIC')
            targets = dataset.createDimension('targets',len_obs)
            lev = dataset.createDimension('lev',len2)
            lev_var = dataset.createVariable('lev',np.int32,('lev'))
            lev_var[:] = len2
            targets_var = dataset.createVariable('targets',np.int32,('targets'))
            targets_var[:] = len_obs
            altitude = dataset.createVariable('altitude',np.float32,('targets','lev'))
            altitude[:,:] = infile[dname[0]][i_crt[i]:i_crt[i+1],:]
            averagingkernel = dataset.createVariable('averagingkernel',np.float32,('targets','lev','lev'))
            averagingkernel[:,:,:] = infile[dname[1]][i_crt[i]:i_crt[i+1],:,:]
            constraintvector = dataset.createVariable('constraintvector',np.float32,('targets','lev'))
            constraintvector[:,:] = infile[dname[7]][i_crt[i]:i_crt[i+1],:]
            latitude = dataset.createVariable('latitude',np.float32,('targets'))
            latitude[:] = infile[dname[8]][i_crt[i]:i_crt[i+1]]
            longitude = dataset.createVariable('longitude',np.float32,('targets'))
            longitude[:] = infile[dname[9]][i_crt[i]:i_crt[i+1]]
            observationerrorcovariance = dataset.createVariable('observationerrorcovariance',np.float32,('targets','lev','lev'))
            observationerrorcovariance[:,:,:] = infile[dname[3]][i_crt[i]:i_crt[i+1],:,:]
            pressure = dataset.createVariable('pressure',np.float32,('targets','lev'))
            pressure[:,:] = infile[dname[4]][i_crt[i]:i_crt[i+1],:]
            species = dataset.createVariable('species',np.float32,('targets','lev'))
            species[:,:] = infile[dname[2]][i_crt[i]:i_crt[i+1],:]
            speciesretrievalconverged = dataset.createVariable('speciesretrievalconverged',np.int32,('targets'))
            speciesretrievalconverged[:] = infile[dname[5]][i_crt[i]:i_crt[i+1]].astype(int)
            time = dataset.createVariable('time',np.float32,('targets'))
            time[:] = infile[dname[10]][i_crt[i]:i_crt[i+1]]
            yyyymmdd = dataset.createVariable('yyyymmdd',np.float32,('targets'))
            yyyymmdd[:] =ymd[i_crt[i]:i_crt[i+1]]
            dataset.close()

print "done"
