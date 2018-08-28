#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate

def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
#rc('text', usetex=True)
#rc('font',size=12)
fig = plt.figure(1,figsize=(12,7.5))
molCO = 6.022*10**23 / ( 68.02 )
#CO is reported to be either 46 or 40.03
s2d = 60*60*24
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
mon_counter = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon = infile00.variables["LON"][:]
len2 = len(lon)
lat = infile00.variables["LAT"][:]
len1 = len(lat)
len3 = 4
len4 = 24
CO_ab = np.zeros((len3,2,24,len1,len2))
emiss_mop = np.zeros((len3,len4,len1,len2))

for j in range(4):
    for i in range(12):
        month = str(i+1).zfill(2)
        month2 = str(i+21).zfill(2)
        index1 = str(8).zfill(2)
        index2 = str(28).zfill(2)
        runbuffer = '/runs/v8-02-01/'
        if j==0:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_M_V35/'
            rundir = 'geos5_ap_20'
            mondir1 = '16'
            mondir2 = '16'
        elif j==1:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_O_V35/'
            rundir = 'geos5_omh_16'
            mondir1 = month
            mondir2 = month2
            sf_ind = '10'
        elif j==2:
            sysdir = '/users/jk/15/xzhang/'
            codedir = 'gcadj_std_T_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_07'
            mondir2 = month+'_27'
            sf_ind = '20'
        elif j==3:
            sysdir = '/users/jk/16/xzhang/'
            codedir = 'gcadj_std_I_3d'
            rundir = 'geos5_all_16'
            mondir1 = month+'_08'
            mondir2 = month+'_28'
            sf_ind = '20'
        path = sysdir+codedir+runbuffer+rundir
        path_ap = path+mondir1+'/ctm.01.2016'+month+'01.nc'
        path_ap2 = path+mondir2+'/ctm.01.2016'+month+'16.nc'
        spec_ap = Dataset(path_ap)
        spec_ap2 = Dataset(path_ap2)
        if j>0:
            path_sf = path+mondir1+'/gctm.sf.'+sf_ind+'.2016'+month+'01.nc'
            path_sf2 = path+mondir2+'/gctm.sf.'+sf_ind+'.2016'+month+'16.nc'
            if os.path.isfile(path_sf) == False:
                sf_ind = '15'
                path_sf =  path+mondir1+'/gctm.sf.'+sf_ind+'.2016'+month+'01.nc'
            if os.path.isfile(path_sf2) == False:
                sf_ind = '15'
                path_sf2 =  path+mondir1+'/gctm.sf.'+sf_ind+'.2016'+month+'16.nc'
            spec_in = Dataset(path_sf)
            spec_in2 = Dataset(path_sf2)
#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
        ps = spec_ap.variables["DAO-FLDS__PS-PBL"][:,:]
        surface = spec_ap.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
        CO_ab[j,0,2*i,:,:] = spec_ap.variables["BIOGSRCE__ISOP"][0,:,:]
        CO_ab[j,0,2*i+1,:,:] = spec_ap2.variables["BIOGSRCE__ISOP"][0,:,:]
        if j==1:
            CO_ab[j,1,2*i,:,:] = spec_in.variables["IJ-EMS-S__ISOP_an"][:,:]
            CO_ab[j,1,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__ISOP_an"][:,:]
        elif j>1:
            CO_ab[j,1,2*i,:,:] = spec_in.variables["IJ-EMS-S__ISOP_an"][0,:,:]
            CO_ab[j,1,2*i+1,:,:] = spec_in2.variables["IJ-EMS-S__ISOP_an"][0,:,:]
        if j==0:
            CO_ems[j,:,:,:] = CO_ab[j,0,:,:,:]
        else:
            CO_ems[j,:,:,:] = CO_ab[j,0,:,:,:]*CO_sf[j,1,:,:,:]
#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx

sum_mop = np.zeros((9,len3,len4))
fig = plt.figure(1,figsize=(16,7.5))
x_series = arange(9)
title1 = ['Southern American','Southern African','Southeast Asian','East Asia',\
          'European','North American','Australian','Northern African','Global']
lon_left=np.array([-80,10,95,95,-10,-130,100,-10,-180])
lon_right = np.array([-40,50,160,150,60,-60,160,50,180])
lat_low = np.array([-50,-30,-10,20,35,10,-50,0,-90])
lat_high = np.array([10,0,10,50,70,60,-10,30,90])
### Compute Southern America Biomass Burning Emissions ### /Or Northern Africa depends on the regions
for k in range(9):
    lon_mask = (lon>=lon_left[k]) & (lon<=lon_right[k])
    lat_mask = (lat>=lat_low[k]) & (lat<=lat_high[k])    
    for i in range(len3):
        for j in range(len4):
            sum_mop[k,i,j] = (emiss_mop*surface)[i,j,lat_mask,:][:,lon_mask].sum()/molISOP*s2d*d2y[j]/10**12
    ax = plt.subplot(331+k)
    plot(x_series,sum_mop[k,0,:],x_series,sum_mop[k,1,:],x_series,sum_mop[k,2,:],linewidth=2.0)
    title('ISOP '+title1[k]+' emissions')
    xticks([0,2,4,6,8,10,12,14,16,18,20,22],['J','F','M','A','M','J','J','A','S','O','N','D'])
    if k==8:### World ###
        legend(('CTRL','MOPITT only','All instrument'),loc='upper right')         
show()
