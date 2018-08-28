#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
import os

#regulate regional structures
def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
fig = plt.figure(1,figsize=(12,7.5))
molCO = 6.022*10**23 / ( 28.0 )
#CO is reported to be either 46 or 40.03
s2d = 60*60*24
#month indicator from Jan-Dec (v1), March-Feb(v2)
#d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
d2y = np.array([15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16,15,16,14,14])
#mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
mon_counter = np.array([31,30,31,30,31,31,30,31,30,31,31,28])
#sample file to read in lat, lon, pressure information
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon = infile00.variables["LON"][:]
len2 = len(lon)
lat = infile00.variables["LAT"][:]
ps = infile00.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile00.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1 = len(lat)
len3 = 4
len4 = 24
CO_ab = np.zeros((len3,4,24,len1,len2))
emiss_mop = np.zeros((3,len3,len4,len1,len2))
#lists of strings for different file entries
sysdir = ['/users/jk/16/xzhang/','/users/jk/15/xzhang/']
codedir = ['gcadj_std_M_V35/','gcadj_std_M_V35/','gcadj_std_T_3d','gcadj_std_I_3d']
runbuffer = '/runs/v8-02-01/'
rundir = ['geos5_ap_2016','geos5_mop_16','geos5_all_16','geos5_all_16']
rundir2 = ['geos5_ap_2017','geos5_mop_17','geos5_all_17','geos5_all_17']
sf_ind = ['10','20','20']
#retrieve emission file from different diretories
for j in range(4):
    for i in range(12):
        if i<10:
            rund = rundir[j]
            yeard = '.2016'
            month = str(i+3).zfill(2)
        else:
            rund = rundir2[j]
            yeard = '.2017'
            month = str(i-9).zfill(2)
        mondir1 = [month,month+'_07',month+'_08']
        mondir2 = [str(int(month)+20).zfill(2),month+'_27',month+'_28']
        if j!=2:
            sysd = sysdir[0]
        else:
            sysd = sysdir[1]
        path = sysd+codedir[j]+runbuffer+rund
        if j==0:
            path_ap = path+'/ctm.00'+yeard+month+'01.nc'
            path_ap2 = path+'/ctm.00'+yeard+month+'16.nc'
        else:
            path_ap = path+mondir1[j-1]+'/ctm.01'+yeard+month+'01.nc'
            path_ap2 = path+mondir2[j-1]+'/ctm.01'+yeard+month+'16.nc'
        spec_ap = Dataset(path_ap)
        spec_ap2 = Dataset(path_ap2)
        if j>0:
            path_sf = path+mondir1[j-1]+'/gctm.sf.'+sf_ind[j-1]+yeard+month+'01.nc'
            path_sf2 = path+mondir2[j-1]+'/gctm.sf.'+sf_ind[j-1]+yeard+month+'16.nc'
            if os.path.isfile(path_sf) == False:
                sf_ind2 = '15'
                path_sf =  path+mondir1[j-1]+'/gctm.sf.'+sf_ind2+yeard+month+'01.nc'
            if os.path.isfile(path_sf2) == False:
                sf_ind2 = '15'
                path_sf2 =  path+mondir2[j-1]+'/gctm.sf.'+sf_ind2+yeard+month+'16.nc'
            spec_in = Dataset(path_sf)
            spec_in2 = Dataset(path_sf2)
#NOX-BIOB__NOx or ANTHSRCE__NOx or NOX-an-S__NOx or NOX-LI-S__NOx
#The array includes anthropogenic and biomass burning component of the CO emissions
        CO_ab[j,0,2*i,:,:] = spec_ap.variables["CO--SRCE__COanth"][0,:,:]
        CO_ab[j,0,2*i+1,:,:] = spec_ap2.variables["CO--SRCE__COanth"][0,:,:]
        CO_ab[j,1,2*i,:,:] = spec_ap.variables["CO--SRCE__CObb"][0,:,:]
        CO_ab[j,1,2*i+1,:,:] = spec_ap2.variables["CO--SRCE__CObb"][0,:,:]
        if j==1:
            CO_ab[j,2,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_an"][:,:]
            CO_ab[j,2,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_an"][:,:]
            CO_ab[j,3,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_bb"][:,:]
            CO_ab[j,3,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_bb"][:,:]
        elif j>1:
            CO_ab[j,2,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_an"][0,:,:]
            CO_ab[j,2,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_an"][0,:,:]
            CO_ab[j,3,2*i,:,:] = spec_in.variables["IJ-EMS-S__CO_bb"][0,:,:]
            CO_ab[j,3,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__CO_bb"][0,:,:]
    if j==0:
        emiss_mop[0,j,:,:,:] = CO_ab[j,0,:,:,:]
        emiss_mop[1,j,:,:,:] = CO_ab[j,1,:,:,:]
    elif j==1:
        emiss_mop[0,j,:,:,:] = CO_ab[j,0,:,:,:]*(0.9*CO_ab[j,2,:,:,:])
        emiss_mop[1,j,:,:,:] = CO_ab[j,1,:,:,:]*(0.9*CO_ab[j,3,:,:,:])
    else:
        emiss_mop[0,j,:,:,:] = CO_ab[j,0,:,:,:]*CO_ab[j,2,:,:,:]
        emiss_mop[1,j,:,:,:] = CO_ab[j,1,:,:,:]*CO_ab[j,3,:,:,:]
emiss_mop[2,:,:,:,:] = emiss_mop[0,:,:,:,:] + emiss_mop[1,:,:,:,:]
#ems_mop = np.zeros((len3,len4/2,len1,len2))
#for mon_i in range(len4/2):
    #ems_mop[:,mon_i,:,:] = np.mean(emiss_mop[:,2*mon_i:2*(mon_i+1),:,:],axis=1)
sum_mop = np.zeros((3,9,len3,len4))
sum_mop_mon = np.zeros((3,9,len3,len4/2))
fig = plt.figure(1,figsize=(20,14.5))
x_series = arange(12)
title1 = ['Southern American','Southern African','Indian and Southeast Asian','East Asian',\
          'Boreal Euroasian','North American','Australian','Northern African','Global']
#lat: SE Asia only: -10,10; SE Asia + India: -10,20; North America: 10,60; US 48 states: 30,50
#lat: Boreal Euroasia: 55,80; Europe: 35,70
#lon: Boreal Euroasia: -30-180; Europe:-10-60
lon_left=np.array([-80,10,95,95,-30,-130,100,-10,-180])
lon_right = np.array([-40,50,160,150,180,-60,160,50,180])
lat_low = np.array([-50,-30,-10,20,55,10,-50,0,-90])
lat_high = np.array([10,0,20,50,80,80,-10,30,90])
title2 = ['anth','bb','total']
### Compute time series of regional CO emissions
for m in range(3):
    for k in range(9):
        lon_mask = (lon >= lon_left[k]) & (lon <= lon_right[k])
        lat_mask = (lat >= lat_low[k]) & (lat<= lat_high[k])            
        for i in range(len3):
            for j in range(len4):
                sum_mop[m,k,i,j] = (emiss_mop*surface)[m,i,j,lat_mask,:][:,lon_mask].sum()/molCO*s2d*d2y[j]/10**12
                sum_mop[m,k,i,j] = (emiss_mop*surface)[m,i,j,lat_mask,:][:,lon_mask].sum()/molCO*s2d*d2y[j]/10**12
        for k2 in range(len4/2):
            sum_mop_mon[m,k,:,k2] = np.sum(sum_mop[m,k,:,2*k2:2*(k2+1)],axis=1)
        sum_reg = np.sum(sum_mop[m,k,:,:],axis=1)
        fig = plt.figure(m,figsize=(18,12))
        ax = plt.subplot(331+k)
        plot(x_series,sum_mop_mon[m,k,0,:],x_series,sum_mop_mon[m,k,1,:],x_series,sum_mop_mon[m,k,2,:],linewidth=2.0)
        title('CO '+title1[k]+' '+title2[m]+' emissions')
        #xticks([0,1,2,3,4,5,6,7,8,9,10,11],['J','F','M','A','M','J','J','A','S','O','N','D'])
        xticks([0,1,2,3,4,5,6,7,8,9,10,11],['M','A','M','J','J','A','S','O','N','D','J','F'])
        if k==8:
            legend(('CTRL:'+str(int(sum_reg[0])),'MOP:'+str(int(sum_reg[1])),'MSA:'+str(int(sum_reg[3]))),\
                    loc='upper right',borderpad=0.3,prop={'size':11},frameon=False)
        else:
            legend((str(int(sum_reg[0])),str(int(sum_reg[1])),str(int(sum_reg[3]))),loc='upper right',\
                    borderpad=0.3,prop={'size':11},frameon=False)
show()
