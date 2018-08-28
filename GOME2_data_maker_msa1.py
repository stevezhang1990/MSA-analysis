#!/usr/bin/python
from scipy import *
from pylab import *
from netCDF4 import Dataset
import glob
import os
from get_intmap import *

path = '/users/jk/16/xzhang/GOME2_NO2/'
year = np.array([2009, 2010])
month = np.array([11,7])
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_lat,len_lon = shape(infile00.variables["IJ-AVG-S__CO"][0,:,:])

obs_lvl = 34
len_dtime = 26*24#366*24
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
len_atime = 26#366
data_raw = np.zeros((len_dtime,len_lat,len_lon,79))
day0 = np.array([2,2-13])
#retrieval method Y = Atrop*Xgome2 Atrop=Atotal*AMFtrop/AMFtotal

for month_range in range(2):
    YYYY = str(year[month_range]).zfill(4)
    MM = str(month[month_range]).zfill(2)
    for filename in glob.glob(os.path.join(path+YYYY+'/'+MM+'/', '*.nc')):
        infile = Dataset(filename)
        print filename
        date = infile.variables["date_orbit"][:]
        day = date - year[month_range]*10000 - month[month_range]*100
        hour = infile.variables["time_orbit"][:]/10000
        lon = (infile.variables["lon_orbit"][:]+180)%360-180
        lat = infile.variables["lat_orbit"][:]
        sza = infile.variables["solar_zenith_angle_orbit"][:]
        vza = infile.variables["viewing_zenith_angle_orbit"][:]    
        cloud_frac = infile.variables["cloud_fraction_orbit"][:]
        len_orbit = len(date)
        amf_total = infile.variables["total_AMF_orbit"][:]
        amf_trop = infile.variables["trop_AMF_orbit"][:]
        ak_total = infile.variables["averaging_kernels_orbit"][:,:]
        flag_trop = infile.variables["trop_flag_orbit"][:]
        pres = infile.variables["pressure_levels_orbit"][:,:]
        vcdt = infile.variables["trop_VCD_orbit"][:]
        vcd = infile.variables["total_VCD_orbit"][:]
        sigvcdt = infile.variables["trop_VCD_error_orbit"][:]
        sigvcd = infile.variables["total_VCD_error_orbit"][:]
        trop_lvl = infile.variables["tropopause_level"][:]
        len_pedge = len(pres[0,:])
        len_pcenter = len(ak_total[0,:])  
        for i in range(len_orbit):
            if flag_trop[i] == 0 and sza[i] < 85 and vza[i]< 75 and \
                abs(lat[i]<80) and cloud_frac[i] <= 0.2 and cloud_frac[i] >0 \
                and amf_total[i]>0 and vcdt[i] > 0 and vcd[i] > 0 and \
                vcdt[i]< vcd[i] and sigvcdt[i] > 0 and sigvcd[i] > 0:
                #print "check"
                time_i = 24*(day[i]-day0[month_range])+int(hour[i])
                lat_i,lon_j = get_gc_ij(lat[i],lon[i])
                data_raw[time_i,lat_i,lon_j,0] = data_raw[time_i,lat_i,lon_j,0] + 1
                data_raw[time_i,lat_i,lon_j,1] = data_raw[time_i,lat_i,lon_j,1] + amf_total[i]
                data_raw[time_i,lat_i,lon_j,2] = data_raw[time_i,lat_i,lon_j,2] + amf_trop[i]
                data_raw[time_i,lat_i,lon_j,3:37] = data_raw[time_i,lat_i,lon_j,3:37] + ak_total[i,:]
                data_raw[time_i,lat_i,lon_j,37:72] = data_raw[time_i,lat_i,lon_j,37:72] + pres[i,:]/100
                data_raw[time_i,lat_i,lon_j,72] = data_raw[time_i,lat_i,lon_j,72] + vcdt[i]*(1-sigvcdt[i]/vcdt[i])**2
                data_raw[time_i,lat_i,lon_j,73] = data_raw[time_i,lat_i,lon_j,73] + vcd[i]*(1-sigvcd[i]/vcd[i])**2
                data_raw[time_i,lat_i,lon_j,74] = data_raw[time_i,lat_i,lon_j,74] + sigvcdt[i]
                data_raw[time_i,lat_i,lon_j,75] = data_raw[time_i,lat_i,lon_j,75] + sigvcd[i]
                data_raw[time_i,lat_i,lon_j,76] = data_raw[time_i,lat_i,lon_j,76] + trop_lvl[i]
                data_raw[time_i,lat_i,lon_j,77] = data_raw[time_i,lat_i,lon_j,77] + (1-sigvcdt[i]/vcdt[i])**2
                data_raw[time_i,lat_i,lon_j,78] = data_raw[time_i,lat_i,lon_j,78] + (1-sigvcd[i]/vcd[i])**2
dataset = Dataset(path+'gome2_msa1.nc','w',format='NETCDF4_CLASSIC')
len_time = dataset.createDimension('len_time',len_dtime)
len_pedge = dataset.createDimension('len_pedge',len_pedge)
len_pcenter = dataset.createDimension('len_pcenter',len_pcenter)
len_lat = dataset.createDimension('len_lat',len_lat)
len_lon = dataset.createDimension('len_lon',len_lon)
amf_total_gc = dataset.createVariable('amf_total_gc',np.float32,('len_time','len_lat','len_lon'))
amf_total_gc[:,:,:] = data_raw[:,:,:,1]/data_raw[:,:,:,0]
amf_trop_gc = dataset.createVariable('amf_trop_gc',np.float32,('len_time','len_lat','len_lon'))
amf_trop_gc[:,:,:] = data_raw[:,:,:,2]/data_raw[:,:,:,0]
ak_total_gc = dataset.createVariable('ak_total_gc',np.float32,('len_time','len_lat','len_lon','len_pcenter'))
for k in range(34):
    ak_total_gc[:,:,:,k] = data_raw[:,:,:,3+k]/data_raw[:,:,:,0]
pres_gc = dataset.createVariable('pres_gc',np.float32,('len_time','len_lat','len_lon','len_pedge'))
for k in range(35):
    pres_gc[:,:,:,k] = data_raw[:,:,:,37+k]/data_raw[:,:,:,0]
vcdt_gc = dataset.createVariable('vcdt_gc',np.float32,('len_time','len_lat','len_lon'))
vcdt_gc[:,:,:] = data_raw[:,:,:,72]/data_raw[:,:,:,77]
#vcdt_gc[:,:,:] = data_raw[:,:,:,72]/data_raw[:,:,:,0]
vcd_gc = dataset.createVariable('vcd_gc',np.float32,('len_time','len_lat','len_lon'))
vcd_gc[:,:,:] = data_raw[:,:,:,73]/data_raw[:,:,:,78]
sigvcdt_gc = dataset.createVariable('sigvcdt_gc',np.float32,('len_time','len_lat','len_lon'))
sigvcdt_gc[:,:,:] = data_raw[:,:,:,74]/data_raw[:,:,:,0]
sigvcd_gc = dataset.createVariable('sigvcd_gc',np.float32,('len_time','len_lat','len_lon'))
sigvcd_gc[:,:,:] = data_raw[:,:,:,75]/data_raw[:,:,:,0]
tropl_gc = dataset.createVariable('trop_level_gc',np.int32,('len_time','len_lat','len_lon'))
tropl_gc[:,:,:] = (data_raw[:,:,:,76]/data_raw[:,:,:,0]).astype(int)
dataset.close()
print "done"
