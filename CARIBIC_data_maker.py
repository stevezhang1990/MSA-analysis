#!/usr/bin/python
from scipy import *
from pylab import *
import csv
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


data1 = Dataset('/users/jk/07/xzhang/CARIBIC/CARIBIC_10s_data_20090526_270_FRA_DEN_10s_V19.nc')
data2 = Dataset('/users/jk/07/xzhang/CARIBIC/CARIBIC_10s_data_20090527_271_DEN_FRA_10s_V18.nc')
data3 = Dataset('/users/jk/07/xzhang/CARIBIC/CARIBIC_10s_data_20090527_272_FRA_KIX_10s_V20.nc')
data4 = Dataset('/users/jk/07/xzhang/CARIBIC/CARIBIC_10s_data_20090528_273_KIX_FRA_10s_V18.nc')
time1 = data1.variables["time"][:]
time2 = data2.variables["time"][:]
time3 = data3.variables["time"][:]
time4 = data4.variables["time"][:]
len1 = len(time1)
len2 = len(time2)
len3 = len(time3)
len4 = len(time4)
len_data = len1+len2+len3+len4
time = np.zeros((len_data))
year = np.zeros((len_data))
doy = np.zeros((len_data))
month = np.zeros((len_data))
day = np.zeros((len_data))
lon = np.zeros((len_data))
lat = np.zeros((len_data))
pres = np.zeros((len_data))
O3 = np.zeros((len_data))
NO = np.zeros((len_data))
NO2 = np.zeros((len_data))
NOy = np.zeros((len_data))

for i in range(len_data):
    year[i] = 2009
    month[i] = 5
    if i < len1:
        day[i] = 26
        time[i] = time1[i]
        lon[i] = data1.variables["PosLong"][i]
        lat[i] = data1.variables["PosLat"][i]
        pres[i] = data1.variables["pstatic"][i]
        O3[i] = data1.variables["Ozone"][i]
        NO[i] = data1.variables["NO"][i]
        NO2[i] = data1.variables["NO2"][i]
        NOy[i] = data1.variables["NOy"][i]
    elif i < len1+len2:
        day[i] = 27
        time[i] = time2[i-len1]
        lon[i] = data2.variables["PosLong"][i-len1]
        lat[i] = data2.variables["PosLat"][i-len1]
        pres[i] = data2.variables["pstatic"][i-len1]
        O3[i] = data2.variables["Ozone"][i-len1]
        NO[i] = data2.variables["NO"][i-len1]
        NO2[i] = data2.variables["NO2"][i-len1]
        NOy[i] = data2.variables["NOy"][i-len1]
    elif i < len1+len2+len3:
        day[i] = 27
        time[i] = time3[i-len1-len2]
        lon[i] = data3.variables["PosLong"][i-len1-len2]
        lat[i] = data3.variables["PosLat"][i-len1-len2]
        pres[i] = data3.variables["pstatic"][i-len1-len2]
        O3[i] = data3.variables["Ozone"][i-len1-len2]
        NO[i] = data3.variables["NO"][i-len1-len2]
        NO2[i] = data3.variables["NO2"][i-len1-len2]
        NOy[i] = data3.variables["NOy"][i-len1-len2]
    else:
        time[i] = time4[i-len1-len2-len3]
        lon[i] = data4.variables["PosLong"][i-len1-len2-len3]
        lat[i] = data4.variables["PosLat"][i-len1-len2-len3]
        pres[i] = data4.variables["pstatic"][i-len1-len2-len3]
        O3[i] = data4.variables["Ozone"][i-len1-len2-len3]
        NO[i] = data4.variables["NO"][i-len1-len2-len3]
        NO2[i] = data4.variables["NO2"][i-len1-len2-len3]
        NOy[i] = data4.variables["NOy"][i-len1-len2-len3]
        day[i] = 28

hours = floor(time/(3600.))
hours[hours>=24] -=24
mins = floor(time/60.)%60
utc = time
utc = utc/(3600.)
doy[utc>=24] += 1
date = hours*100 + mins
doy = day + doy
m=zeros((len_data,8))
#print O3

print m.shape
m[:,0] = year
m[:,1] = month
m[:,2] = doy
m[:,3] = date
m[:,4] = lon
m[:,5] = lat
m[:,6] = pres
m[:,7] = np.nan_to_num(NO)
l = m.tolist()

outkeys = ["Year","Month","DOY","UTC","LON","LAT","PRES","NO"]

out = open('/users/jk/12/xzhang/gcadj_std_O_OSI/runs/v8-02-01/geos5_osino2_0905_2/caribic_no_200905.csv','w')
out_csv = csv.writer(out,delimiter="\t")

for row in l:
    out_csv.writerow(row)
out.close()

print "done"
