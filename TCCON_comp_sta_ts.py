import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

data = Dataset('/users/jk/16/xzhang/TCCON/tccon_gc_2016_feng.nc')
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_time,len_lat,len_lon = shape(data.variables["ts_co_ap_tccon"][:,:,:])
len_time1 = len_time-6*24
data_tccon = np.zeros((5,len_time1,len_lat,len_lon))
data_tccon[0,:,:,:] = data.variables["ts_co_ap_tccon"][0:len_time1,:,:]
data_tccon[1,:,:,:] = data.variables["ts_co_mop_tccon"][0:len_time1,:,:]
data_tccon[2,:,:,:] = data.variables["ts_co_msa05_tccon"][0:len_time1,:,:]
data_tccon[3,:,:,:] = data.variables["ts_co_msa06_tccon"][0:len_time1,:,:]
data_tccon[4,:,:,:] = data.variables["ts_co_obs_tccon"][0:len_time1,:,:]
len_time_plot = 24
len5 = int(len_time1/len_time_plot)
data_tccon_reshape = ma.masked_invalid(np.reshape(data_tccon,(5,len5,len_time_plot,len_lat,len_lon)))
tccon_reshape = np.ma.masked_where(abs(data_tccon_reshape) == 0, data_tccon_reshape)
co_tccon = np.mean(tccon_reshape,axis=1)
co_tccon_total = np.zeros((3,5,len_time_plot))
co_tccon_count = np.zeros((3,5,len_time_plot))
co_tccon_plot = np.zeros((3,5,len_time_plot))
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
mon_counter=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
model_hour=24
#len_time = model_hour*np.sum(mon_counter[0:4])
time = arange(len_time_plot)
for exp_i in range(5):
   for i in range(len_lat):
      for j in range(len_lon):
         if np.sum(co_tccon[exp_i,:,i,j])>0:
            if lat_gc[i] > 60:
               for k in range(len_time_plot):
                  if co_tccon[exp_i,k,i,j]>0:
                     if k>=18 and exp_i>0 and exp_i<3:
                        co_tccon_total[0,exp_i,k] = co_tccon_total[0,exp_i,k] + co_tccon[exp_i,k,i,j]
                     else:
                        co_tccon_total[0,exp_i,k] = co_tccon_total[0,exp_i,k] + co_tccon[exp_i,k,i,j]
                     co_tccon_count[0,exp_i,k] = co_tccon_count[0,exp_i,k] + 1
            elif lat_gc[i] > 0:
               for k in range(len_time_plot):
                  if co_tccon[exp_i,k,i,j]>0:
                     if exp_i<3:
                        co_tccon_total[1,exp_i,k] = co_tccon_total[1,exp_i,k] + co_tccon[exp_i,k,i,j]
                     else:
                        co_tccon_total[1,exp_i,k] = co_tccon_total[1,exp_i,k] + co_tccon[exp_i,k,i,j]
                     co_tccon_count[1,exp_i,k] = co_tccon_count[1,exp_i,k] + 1
            else:
               for k in range(len_time_plot):
                  if co_tccon[exp_i,k,i,j]>0:
                     if exp_i <4:
                        co_tccon_total[2,exp_i,k] = co_tccon_total[2,exp_i,k] + co_tccon[exp_i,k,i,j]
                     else:
                        co_tccon_total[2,4,k] = co_tccon_total[2,4,k] + co_tccon[4,k,i,j]
                     co_tccon_count[2,exp_i,k] = co_tccon_count[2,exp_i,k] + 1            

co_tccon_plot = co_tccon_total/co_tccon_count
fig = plt.figure(1,figsize=(12,9))
for plot_i in range(3):
   plot_frame = int(plot_i+311)
   ax = subplot(plot_frame)
   if plot_i ==0:
      title(r'TCCON XCO comparison, 60$^o$N-90$^o$N (2 sites)',fontsize=11)
   elif plot_i == 1:
      title(r'TCCON XCO comparison, 20$^o$N-60$^o$N (14 sites)',fontsize=11)
   elif plot_i == 2:
      title(r'TCCON XCO comparison, 0$^o$S-90$^o$S (3 sites)',fontsize=11)
   plot(time,1.1*co_tccon_plot[plot_i,0,:],'k--',time,co_tccon_plot[plot_i,1,:],time,co_tccon_plot[plot_i,2,:],\
         time,co_tccon_plot[plot_i,4,:],linewidth=2.0)
   xticks([0,2,4,6,8,10,12,14,16,18,20,22],['J','F','M','A','M','J','J','A','S','O','N','D'])
   if plot_i == 0:
      legend(('A priori','MOPITT','All instrument','TCCON'),loc='upper left',borderpad=0.3,prop={'size':11},frameon=False)
show()

print "check"
