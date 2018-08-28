import os
import numpy as np
from get_intmap import *
from netCDF4 import Dataset
from scipy import *
from pylab import *
# SHADOZ observations
data = Dataset('/users/jk/16/xzhang/SHADOZ/shadoz_obs_2016.nc')
# sample file reading modeled lat lon and pressure info
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
ts_pres = Dataset('/users/jk/16/xzhang/TS_annual/geos5_ap_2016/ts_pres_2016.nc')
len_station,len_time,len_lvl = shape(data.variables["shadoz_o3_obs"][:,:,:])
len_time1 = len_time#-6*24
# read in shadoz observation projection on modeled grid for obs, CTRL and assimilation run.
data_shadoz = np.zeros((10,len_station,len_time1,len_lvl))
data_shadoz[0,:,:,:] = data.variables["shadoz_o3_ap"][:,0:len_time1,:]
data_shadoz[1,:,:,:] = data.variables["shadoz_o3_iasio3"][:,0:len_time1,:]
data_shadoz[2,:,:,:] = data.variables["shadoz_o3_msa07"][:,0:len_time1,:]
data_shadoz[3,:,:,:] = data.variables["shadoz_o3_msa08"][:,0:len_time1,:]
data_shadoz[4,:,:,:] = data.variables["shadoz_o3_obs"][:,0:len_time1,:]
lat_gc = data.variables["lat_obs"][:]
lon_gc = data.variables["lon_obs"][:]
ts_p_center = np.mean(np.ma.masked_invalid(ts_pres.variables["ts_p_center"][:,:,:,:]),axis=0)

#for i in range(5):
   #data_shadoz[5+i,:,:,:] = 100*(data_shadoz[i,:,:,:]-data_shadoz[4,:,:,:])/data_shadoz[4,:,:,:]
# masking and averaging to compute annual mean
diff_shadoz_filt = ma.masked_invalid(data_shadoz[:5,:,:,:])
diff_shadoz2 = np.ma.masked_where(abs(diff_shadoz_filt)==0,diff_shadoz_filt)
diff_shadoz = np.ma.masked_where(abs(diff_shadoz2)>300,diff_shadoz2)
#o3_shadoz = np.zeros((5,len_station,len_lvl))
#bias correction on modeled and observed O3
o3_shadoz= np.mean(diff_shadoz,axis=2)
o3_shadoz[3,5,:] = o3_shadoz[3,5,:]-15
k0 = 0.5
k1=8
k2=23
for i in range(len_lvl):
   for j in range(len_station):
      if i>=k2 and j!=5 and i<=k2+int(k1/k0):
         o3_shadoz[3,j,i] = o3_shadoz[3,j,i]-k0*(i-k2)
      elif i>k2+int(k1/k0) and j!=5:
         o3_shadoz[3,j,i] = o3_shadoz[3,j,i]-k0*int(k1/k0)
#for i in range(5):
   #o3_shadoz[i,:,:] = 100*(o3_shadoz_raw[i,:,:]-o3_shadoz_raw[4,:,:])/o3_shadoz_raw[4,:,:]
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
mon_counter=[31,29,31,30,31,30,31,31,30,31,30,31]
mon_tick = np.zeros(12)
model_hour=24
#pres_mod = np.zeros(47)
for mon_i in range(len(mon_counter)):
   mon_tick[mon_i] = 24*np.sum(mon_counter[:mon_i])
#len_time = model_hour*np.sum(mon_counter[0:4])
plot_counter = 0
fig = plt.figure(1,figsize=(16,9))
site = ['Irene, S.Africa','La Reunion Is, France','Suva, Fiji','Pago Pago, Samoa',\
        'Natal, Brazil','Ascension Is','San Cristobal, Ecuador','Nairobi, Kenya',\
        'Kuala Lumpur, Malaysia','Paramaribo, Surinam','Costa Rica',\
        'Hilo, HI','Hanoi, Vietnam']
diff_obs = np.zeros(len_lvl)
# plot vertical profile of SHADOZ ozone
props = dict(boxstyle='square', facecolor='white', alpha=0.0)
for plot_counter in range(len_station):
      if plot_counter <9:
         plot_i = 0
         plot_frame = int(-plot_counter+339)
      elif plot_counter <18:
         plot_i = 1
         plot_frame = int(-plot_counter+348)
      else:
         plot_i = 2
         plot_frame = int(-plot_counter+352)
      fig = plt.figure(plot_i,figsize=(16,12))
      if np.sum(o3_shadoz[1,plot_counter,:])!=0:
         fig = plt.figure(plot_i,figsize=(16,12))
         ax=subplot(plot_frame)
         if lat_gc[plot_counter] > 0:
            lat_ind = 'N'
         else:
            lat_ind = 'S'
         if lon_gc[plot_counter] < 0:
            lon_ind = 'W'
         else:
            lon_ind = 'E'
         title(site[plot_counter]+' '+str(int(abs(lat_gc[plot_counter])))+lat_ind+\
               ', '+str(int(abs(lon_gc[plot_counter])))+lon_ind,fontsize=11)
         i = np.where(lat==lat_gc[plot_counter])
         j = np.where(lon==lon_gc[plot_counter])
         pres_mod = np.reshape(ts_p_center[:,i,j],(47))[1:]
         plot(o3_shadoz[0,plot_counter,1:],pres_mod,'k--',o3_shadoz[1,plot_counter,1:],pres_mod,'b-',o3_shadoz[3,plot_counter,1:],pres_mod,'g-',\
              o3_shadoz[4,plot_counter,1:],pres_mod,'ro',linewidth=2.0)
         textstr = '\n'.join((
            'LT, '+'MT&UT',
            'CTRL='+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[0,plot_counter,1:16]==0,o3_shadoz[0,plot_counter,1:16])),decimals=1))\
            +', '+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[0,plot_counter,16:38]==0,o3_shadoz[0,plot_counter,16:38])),decimals=1)),
            'IASI='+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[1,plot_counter,1:16]==0,o3_shadoz[1,plot_counter,1:16])),decimals=1))\
            +', '+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[1,plot_counter,16:38]==0,o3_shadoz[1,plot_counter,16:38])),decimals=1)),
            'MSA='+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[3,plot_counter,1:16]==0,o3_shadoz[3,plot_counter,1:16])),decimals=1))\
            +', '+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[3,plot_counter,16:38]==0,o3_shadoz[3,plot_counter,16:38])),decimals=1)),
            'OBS='+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[4,plot_counter,1:16]==0,o3_shadoz[4,plot_counter,1:16])),decimals=1))\
            +', '+str(np.around(np.mean(np.ma.masked_where(o3_shadoz[4,plot_counter,16:38]==0,o3_shadoz[4,plot_counter,16:38])),decimals=1))))
         ax.text(0.98, 0.05, textstr, transform=ax.transAxes, fontsize=11,horizontalalignment='right',
                 verticalalignment='bottom', bbox=props)
         #xticks(mon_tick,['J','F','M','A','M','J','J','A','S','O','N','D'])
         xlim(0,150)
         ylim(pres_mod[0],100)
         if plot_counter == 4:# or plot_counter == 17:# or plot_counter == 18:
            legend(('A priori','IASI O3','MSA','SHADOZ'),loc='upper left',borderpad=0.3,prop={'size':11},frameon=False)
         plot_counter = plot_counter+1
print "check"
show()

