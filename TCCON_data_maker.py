#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
len_file = 0
path = '/users/jk/15/xzhang/TCCON/'
# find all tccon files
for filename in glob.glob(os.path.join(path, '*.public.nc')):
    len_file = len_file + 1
obs_lvl = 71
zen_lvl = 16
len_dtime = 366*24
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
len_atime = 366
grid_info = np.zeros((len_file,2))
height_info = np.zeros((len_file,2,71))
co_info = np.zeros((len_file,3,len_dtime))
co_info_norm = np.zeros((len_file,3,len_dtime))
co_ak_info = np.zeros((len_file,72,16))
co_ap_info = np.zeros((len_file,2,len_atime,71))
pwt = np.array([0.039219, 0.116322, 0.093536, 0.085170, 0.077237, \
                0.069775, 0.062885, 0.056620, 0.050887, 0.045608, \
                0.040725, 0.036187, 0.031941, 0.028035, 0.024422, \
                0.021054, 0.018096, 0.015448, 0.013107, 0.011082, \
                0.009325, 0.007851, 0.006629, 0.005620, 0.004771, \
                0.004060, 0.003473, 0.002971, 0.002538, 0.002175, \
                0.001861, 0.001593, 0.001360, 0.001158, 0.000988, \
                0.000846, 0.000724, 0.000622, 0.000534, 0.000460, \
                0.000397, 0.000343, 0.000296, 0.000257, 0.000223, \
                0.000195, 0.000170, 0.000148, 0.000130, 0.000115, \
                0.000102, 0.000091, 0.000080, 0.000071, 0.000063, \
                0.000056, 0.000050, 0.000044, 0.000039, 0.000034, \
                0.000030, 0.000027, 0.000023, 0.000020, 0.000018, \
                0.000016, 0.000014, 0.000012, 0.000011, 0.000009, 0.000004])

file_i = 0
for filename in glob.glob(os.path.join(path, '*.public.nc')):
    infile = Dataset(filename)
    year = infile.variables["year"][:]
    day = infile.variables["day"][:]
    hour = infile.variables["hour"][:]
    lon = infile.variables["long_deg"][:]
    lat = infile.variables["lat_deg"][:]
    sza = infile.variables["asza_deg"][:]
    co = infile.variables["xco_ppb"][:]
    len_time = len(year)
    prior_co = infile.variables["prior_co"][:,:]
    prior_height = infile.variables["prior_Height"][:]
    ak_zenith = infile.variables["ak_zenith"][:]
    ak_p = infile.variables["ak_P_hPa"][:]
    ak_co = infile.variables["ak_co"][:,:]
    len_pdate, len_pheight = shape(prior_co)
    len_hpa, len_zenith = shape(ak_co)
    prior_year = infile.variables["prior_year"][:]
    prior_month = infile.variables["prior_month"][:]
    prior_day = infile.variables["prior_day"][:]
    prior_pres = infile.variables["prior_Pressure"][:]
    grid_info[file_i,:] = np.array([lon[0],lat[0]])
    height_info[file_i,0,:] = prior_height
    height_info[file_i,1,:] = ak_p
    co_ak_info[file_i,71,:] = ak_zenith
    co_ak_info[file_i,0:71,:] = ak_co    
    for i in range(len_time):
        if year[i] == 2016:
            j = 24*(day[i]-1)+int(hour[i])
            co_info[file_i,0,j] = co_info[file_i,0,j] + sza[i]
            co_info[file_i,1,j] = co_info[file_i,1,j] + co[i]
            co_info[file_i,2,j] = co_info[file_i,2,j] + 1
    for k in range(len_dtime):
        if co_info[file_i,2,k] > 0:
            co_info_norm[file_i,0,k] = co_info[file_i,0,k]/co_info[file_i,2,k]
            co_info_norm[file_i,1,k] = co_info[file_i,1,k]/co_info[file_i,2,k]
    for ii in range(len_pdate):
        if prior_year[ii] == 2016:
            jj = np.sum(mon_counter[0:int(prior_month[ii]-1)])+prior_day[ii]-1
            co_ap_info[file_i,0,jj,:] = prior_co[ii,:]
            co_ap_info[file_i,1,jj,:] = prior_pres[ii,:]
    file_i = file_i + 1
    print file_i

dataset = Dataset(path+'tccon_2016.nc','w',format='NETCDF4_CLASSIC')
len_station = dataset.createDimension('len_station',len_file)
len_lvl = dataset.createDimension('len_lvl',obs_lvl)
len_zen = dataset.createDimension('len_zen',zen_lvl)
len_dtime = dataset.createDimension('len_time',len_dtime)
len_atime = dataset.createDimension('len_time_ap',len_atime)
lat_tc = dataset.createVariable('station_lat',np.float32,('len_station'))
lat_tc[:] = grid_info[:,1]
lon_tc = dataset.createVariable('station_lon',np.float32,('len_station'))
lon_tc[:] = grid_info[:,0]
prior_height_tc = dataset.createVariable('sta_prior_height',np.float32,('len_station','len_lvl'))
prior_height_tc[:,:] = height_info[:,0,:]
ak_pres_tc = dataset.createVariable('sta_ak_pres',np.float32,('len_station','len_lvl'))
ak_pres_tc[:,:] = height_info[:,1,:]
ak_zenith_tc = dataset.createVariable('sta_ak_zenith',np.float32,('len_station','len_zen'))
ak_zenith_tc[:,:] = co_ak_info[:,71,:]
ak_co_tc = dataset.createVariable('sta_ak_co',np.float32,('len_station','len_lvl','len_zen'))
ak_co_tc[:,:,:] = co_ak_info[:,0:71,:]
xco_tc = dataset.createVariable('sta_xco',np.float32,('len_station','len_time'))
xco_tc[:,:] = co_info_norm[:,1,:]
sza_tc = dataset.createVariable('sta_sza',np.float32,('len_station','len_time'))
sza_tc[:,:] = co_info_norm[:,0,:]
prior_co_tc = dataset.createVariable('sta_prior_co',np.float32,('len_station','len_time_ap','len_lvl'))
prior_co_tc[:,:,:] = co_ap_info[:,0,:,:]
prior_pres_tc = dataset.createVariable('sta_prior_pres',np.float32,('len_station','len_time_ap','len_lvl'))
prior_pres_tc[:,:,:] = co_ap_info[:,1,:,:]
prior_pwt_tc = dataset.createVariable('prior_pwt',np.float32,('len_lvl'))
prior_pwt_tc[:] = pwt
dataset.close()
print "done"
