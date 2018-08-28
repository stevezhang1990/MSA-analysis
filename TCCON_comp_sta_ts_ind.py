import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *

data = Dataset('/users/jk/16/xzhang/TCCON/tccon_gc_2016_feng.nc')
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
len_time,len_lat,len_lon = shape(data.variables["ts_co_ap_tccon"][:,:,:])
len_time1 = len_time#-6*24
data_tccon = np.zeros((5,len_time1,len_lat,len_lon))
data_tccon[0,:,:,:] = data.variables["ts_co_ap_tccon"][0:len_time1,:,:]
data_tccon[1,:,:,:] = data.variables["ts_co_mop_tccon"][0:len_time1,:,:]
data_tccon[2,:,:,:] = data.variables["ts_co_msa05_tccon"][0:len_time1,:,:]
data_tccon[3,:,:,:] = data.variables["ts_co_msa06_tccon"][0:len_time1,:,:]
data_tccon[4,:,:,:] = data.variables["ts_co_obs_tccon"][0:len_time1,:,:]
len_time_plot = len_time
data_tccon_filt = ma.masked_invalid(data_tccon)
tccon_filt = np.ma.masked_where(abs(data_tccon_filt) < 20, data_tccon_filt)
co_tccon = np.ma.masked_where(abs(tccon_filt) > 150, tccon_filt)
lon_gc = infile00.variables["LON"][:]
lat_gc = infile00.variables["LAT"][:]
mon_counter=[31,29,31,30,31,30,31,31,30,31,30,31]
#avg_crt = 6*24
#len_time_day = int(len_time/avg_crt)
len_time_day = 12
co_tccon_day = np.zeros((5,len_time_day,len_lat,len_lon))
model_hour=24
#mon_tick = np.zeros(12)
#for mon_i in range(len(mon_counter)):
   #mon_tick[mon_i] = 24*np.sum(mon_counter[:mon_i])/avg_crt
for day_i in range(len_time_day):
   co_tccon_day[:,day_i,:,:] = np.mean(co_tccon[:,24*np.sum(mon_counter[:day_i]):24*np.sum(mon_counter[:day_i+1]),:,:],axis=1)
#len_time = model_hour*np.sum(mon_counter[0:4])

time = arange(len_time_day)
plot_counter = 0
fig = plt.figure(1,figsize=(16,9))
site = ['Lauder','Wollongong','Reunion','Darwin','Ascension Island','Izana',\
        'JPL','Saga','Lamont','Anmyeondo','Tsukuba','Rikubetsu','Park Falls',\
        'Orleans','Zugspitze','Paris','Karlsruhe','East Trout Lake','Bremen',\
        'Bialystok','Sodanskyla','Ny-Alesund']
site_index = np.array([21,20,18,19,17,16,4,15,5,14,13,12,6,7,8,9,10,0,11,1,2,3])
site1 = ['Boreal forest (4 sites)','North America (3 sites)','Europe (5 sites)',\
         'East Asia (4 sites)','SH Ocean (3 sites)','Oceania (3 sites)']
co_lim = np.array([70,100,100,100,125,110,150,180,140,160,180,140,\
                   140,140,200,145,130,105,130,120,130,120])
co_tccon_gr = np.zeros((6,5,len_time_day))
co_tccon_orb = np.zeros((22,5,len_time_day))
sta_i = 0
for i in range(len_lat):
   for j in range(len_lon):
      if np.sum(co_tccon[4,:,i,j])>0:
         co_tccon_orb[site_index[sta_i],:,:] = co_tccon_day[:,:,i,j]
         sta_i = sta_i + 1

co_tccon_gr[0,:,:] = np.mean(co_tccon_orb[0:4,:,:],axis=0)
co_tccon_gr[1,:,:] = np.mean(co_tccon_orb[4:7,:,:],axis=0)
co_tccon_gr[2,:,:] = np.mean(co_tccon_orb[7:12,:,:],axis=0)
co_tccon_gr[3,:,:] = np.mean(co_tccon_orb[12:16,:,:],axis=0)
co_tccon_gr[4,:,:] = np.mean(co_tccon_orb[16:19,:,:],axis=0)
co_tccon_gr[5,:,:] = np.mean(co_tccon_orb[19:22,:,:],axis=0)
props = dict(boxstyle='square', facecolor='white', alpha=0.0)
fig = plt.figure(1,figsize=(16,9))
for plot_counter in range(6):
   ax=subplot(321+plot_counter)
   title(r'TCCON '+site1[plot_counter],fontsize=11)
   plot(time,co_tccon_gr[plot_counter,0,:],'k+',time,co_tccon_gr[plot_counter,1,:],'.',time,co_tccon_gr[plot_counter,2,:],'.',\
        time,co_tccon_gr[plot_counter,4,:],'ro',linewidth=2.0)
   xticks(time,['J','F','M','A','M','J','J','A','S','O','N','D'])
   textstr = '\n'.join((
      'CTRL='+str(np.around(np.mean(np.ma.masked_where(co_tccon_gr[plot_counter,0,:]==0,co_tccon_gr[plot_counter,0,:])),decimals=1)),
      'MOPITT='+str(np.around(np.mean(np.ma.masked_where(co_tccon_gr[plot_counter,1,:]==0,co_tccon_gr[plot_counter,1,:])),decimals=1)),
      'MSA='+str(np.around(np.mean(np.ma.masked_where(co_tccon_gr[plot_counter,2,:]==0,co_tccon_gr[plot_counter,2,:])),decimals=1)),
      'TCCON='+str(np.around(np.mean(np.ma.masked_where(co_tccon_gr[plot_counter,4,:]==0,co_tccon_gr[plot_counter,4,:])),decimals=1))))
   ax.text(0.98, 0.95, textstr, transform=ax.transAxes, fontsize=11,horizontalalignment='right',
        verticalalignment='top', bbox=props)
   if plot_counter==0:
      legend(('CTRL','MOPITT','MSA','TCCON'),loc='lower center',borderpad=0.3,prop={'size':11},frameon=False, ncol=2)
   xlim([0,12])
   ylim([np.min(co_tccon_gr[plot_counter,:,:])-3,np.max(co_tccon_gr[plot_counter,:,:])+3])
##for i in range(len_lat):
##   for j in range(len_lon):
##      if plot_counter <9:
##         plot_i = 0
##         plot_frame = int(plot_counter+331)
##      elif plot_counter <18:
##         plot_i = 1
##         plot_frame = int(plot_counter+322)
##      else:
##         plot_i = 2
##         plot_frame = int(plot_counter+313)
## 
##      if np.sum(co_tccon[4,:,i,j])>0:
##         fig = plt.figure(plot_i,figsize=(16,9))
##         ax=subplot(plot_frame)
##         title(r'TCCON '+site[plot_counter]+' '+str(lat_gc[i])+', '+str(lon_gc[j]),fontsize=11)
##         plot(time,co_tccon[0,:,i,j],'k+',time,co_tccon[1,:,i,j],'.',time,co_tccon[2,:,i,j],'.',\
##              time,co_tccon[4,:,i,j],'ro',linewidth=2.0)
##         xticks(mon_tick,['J','F','M','A','M','J','J','A','S','O','N','D'])
##         xlim([0,len_time_plot])
##         ylim([40,co_lim[plot_counter]])
##         if plot_counter == 4 or plot_counter == 17:# or plot_counter == 18:
##            legend(('A priori','MOPITT','All instrument','TCCON'),loc='upper left',borderpad=0.3,prop={'size':11},frameon=False)
##         plot_counter = plot_counter+1
print "check"
show()

