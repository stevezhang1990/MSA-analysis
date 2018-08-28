import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

data = Dataset('/users/jk/16/xzhang/TCCON/tccon_gc_msa1_feng.nc')
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_time,len_lat,len_lon = shape(data.variables["ts_co_ap_tccon"][:,:,:])
#len_time1 = len_time-6*24
data_tccon = np.zeros((4,len_time,len_lat,len_lon))
data_tccon[0,:,:,:] = data.variables["ts_co_ap_tccon"][0:len_time,:,:]
data_tccon[1,:,:,:] = data.variables["ts_co_mop_tccon"][0:len_time,:,:]
data_tccon[2,:,:,:] = data.variables["ts_co_msa05_tccon"][0:len_time,:,:]
#data_tccon[3,:,:,:] = data.variables["ts_co_msa06_tccon"][0:len_time1,:,:]
data_tccon[3,:,:,:] = data.variables["ts_co_obs_tccon"][0:len_time,:,:]
data_tccon[1:3,:,39:,:]= 1.2*data_tccon[1:3,:,39:,:]
data_tccon[3,:,24:39,:] = 1.1*data_tccon[3,:,24:39,:]

len_time_plot = len_time
len5 = int(len_time/len_time_plot)
data_tccon_reshape = ma.masked_invalid(np.reshape(data_tccon,(4,len5,len_time_plot,len_lat,len_lon)))
tccon_reshape = np.ma.masked_where(abs(data_tccon_reshape) == 0, data_tccon_reshape)
co_tccon = np.mean(tccon_reshape,axis=1)
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]

mon_counter=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
model_hour=24
#len_time = model_hour*np.sum(mon_counter[0:4])
time = arange((len_time_plot))
fig = plt.figure(1,figsize=(16,9))

plot_counter = 0
for i in range(len_lat):
   for j in range(len_lon):
      if np.sum(co_tccon[3,:,i,j])>0 and plot_counter<9 and lat_gc[i]>50:
         plot_frame = int(plot_counter+331)
         ax=subplot(plot_frame)
         title(r'TCCON XCO comparison, '+str(lat_gc[i])+', '+str(lon_gc[j]),fontsize=11)
         plot(time,co_tccon[0,:,i,j],'k+',time,co_tccon[1,:,i,j],'.',time,co_tccon[2,:,i,j],'.',\
              time,co_tccon[3,:,i,j],'ro',linewidth=2.0)
         xticks([0,24*3,24*7,24*11,24*13,24*16,24*20,24*24],['Nov-2','6','10','14','Jul-2','6','10','14'])
         if plot_counter == 0:
            legend(('A priori','MOPITT','All instrument','TCCON'),loc='upper left',borderpad=0.3,prop={'size':11},frameon=False)
         plot_counter = plot_counter+1
print "check"
show()
