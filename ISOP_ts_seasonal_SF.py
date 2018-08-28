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
molISOP = 6.022*10**23 / ( 68.12 )
#ISOP is reported to be 68.12
s2d = 60*60*24
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
infile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/ctm.01.20160101.nc')
infile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/ctm.01.20160116.nc')
infile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/ctm.01.20160201.nc')
infile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/ctm.01.20160216.nc')
infile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/ctm.01.20160301.nc')
infile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/ctm.01.20160316.nc')
infile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/ctm.01.20160401.nc')
infile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/ctm.01.20160416.nc')
omhfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601/gctm.sf.10.20160101.nc')
omhfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1601_2/gctm.sf.10.20160116.nc')
omhfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602/gctm.sf.10.20160201.nc')
omhfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1602_2/gctm.sf.10.20160216.nc')
omhfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603/gctm.sf.10.20160301.nc')
omhfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1603_2/gctm.sf.10.20160316.nc')
omhfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604/gctm.sf.10.20160401.nc')
omhfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omh_1604_2/gctm.sf.10.20160416.nc')
allfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_06/gctm.sf.20.20160101.nc')
allfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_26/gctm.sf.15.20160116.nc')
allfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_06/gctm.sf.20.20160201.nc')
allfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_26/gctm.sf.20.20160216.nc')
allfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_06/gctm.sf.15.20160301.nc')
allfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_26/gctm.sf.20.20160316.nc')
allfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_06/gctm.sf.15.20160401.nc')
allfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_26/gctm.sf.15.20160416.nc')


#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2 = shape(omhfile01.variables["IJ-EMS-S__ISOP_an"][:,:])
len3 = 3
len4 = 8
ISOP_an = np.zeros((len4,len1,len2))
ISOP_bb = np.zeros((len4,len1,len2))

SF_omh_an = np.zeros((len4,len1,len2))
SF_omh_bb = np.zeros((len4,len1,len2))

SF_all_an = np.zeros((len4,len1,len2))
SF_all_bb = np.zeros((len4,len1,len2))

ISOP_an[0,:,:] = infile01.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[1,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[2,:,:] = infile02.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[3,:,:] = infile21.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[4,:,:] = infile03.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[5,:,:] = infile23.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[6,:,:] = infile04.variables["BIOGSRCE__ISOP"][0,:,:]
ISOP_an[7,:,:] = infile24.variables["BIOGSRCE__ISOP"][0,:,:]

SF_omh_an[0,:,:] = omhfile01.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[1,:,:] = omhfile21.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[2,:,:] = omhfile02.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[3,:,:] = omhfile21.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[4,:,:] = omhfile03.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[5,:,:] = omhfile23.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[6,:,:] = omhfile04.variables["IJ-EMS-S__ISOP_an"][:,:]
SF_omh_an[7,:,:] = omhfile24.variables["IJ-EMS-S__ISOP_an"][:,:]


SF_all_an[0,:,:] = allfile01.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[1,:,:] = allfile21.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[2,:,:] = allfile02.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[3,:,:] = allfile21.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[4,:,:] = allfile03.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[5,:,:] = allfile23.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[6,:,:] = allfile04.variables["IJ-EMS-S__ISOP_an"][0,:,:]
SF_all_an[7,:,:] = allfile24.variables["IJ-EMS-S__ISOP_an"][0,:,:]

ISOP_ap = ISOP_an
ISOP_omh = SF_omh_an
ISOP_all = SF_all_an

emiss_moa = np.zeros((len3,len4,len1,len2))
emiss_mob = np.zeros((len3,len4,len1,len2))

emiss_moa[0,:,:,:] = ISOP_an
emiss_moa[1,:,:,:] = SF_omh_an*ISOP_an
emiss_moa[2,:,:,:] = SF_all_an*ISOP_an

emiss_omh = emiss_moa 

sum_omh = np.zeros((9,len3,len4))

### Compute Southern America Biomass Burning Emissions ### /Or Northern Africa depends on the regions

lon_mask_SAmerica = (lon >= -80) & (lon <= -40)
lat_mask_SAmerica = (lat >= -50) & (lat<= 10)

for i in range(len3):
    for j in range(len4):
        sum_omh[0,i,j] = (emiss_omh*surface)[i,j,lat_mask_SAmerica,:][:,lon_mask_SAmerica].sum()/molISOP*s2d*d2y[j]/10**12

### Compute Southern Africa Biomass Burning Emissions ###

lon_mask_SAfrica = (lon >= 10) & (lon <= 50)
lat_mask_SAfrica = (lat >= -30) & (lat<= 0)
for i in range(len3):
    for j in range(len4):
        sum_omh[1,i,j] = (emiss_omh*surface)[i,j,lat_mask_SAfrica,:][:,lon_mask_SAfrica].sum()/molISOP*s2d*d2y[j]/10**12


### Compute South East Asian Anthropogenic Emissions ###

lon_mask_SE_Asia = (lon >= 95) & (lon <= 160)
lat_mask_SE_Asia = (lat >= -10) & (lat<= 10)

for i in range(len3):
    for j in range(len4):
        sum_omh[2,i,j] = (emiss_omh*surface)[i,j,lat_mask_SE_Asia,:][:,lon_mask_SE_Asia].sum()/molISOP*s2d*d2y[j]/10**12


### Compute East Asian Anthropogenic Emissions ###

lon_mask_Asia = (lon >= 95) & (lon <= 150)
lat_mask_Asia = (lat >= 20) & (lat<= 50)
for i in range(len3):
    for j in range(len4):
        sum_omh[3,i,j] = (emiss_omh*surface)[i,j,lat_mask_Asia,:][:,lon_mask_Asia].sum()/molISOP*s2d*d2y[j]/10**12

### Compute European Anthropogenic Emissions ###

lon_mask_Europe = (lon >= -10) & (lon <= 60)
lat_mask_Europe = (lat >= 35) & (lat<= 70)
for i in range(len3):
    for j in range(len4):
        sum_omh[4,i,j] = (emiss_omh*surface)[i,j,lat_mask_Europe,:][:,lon_mask_Europe].sum()/molISOP*s2d*d2y[j]/10**12

### Compute North American Anthropogenic Emissions ###

lon_mask_NA = (lon >= -130) & (lon <= -60)
lat_mask_NA = (lat >= 10) & (lat<= 70)
for i in range(len3):
    for j in range(len4):
        sum_omh[5,i,j] = (emiss_omh*surface)[i,j,lat_mask_NA,:][:,lon_mask_NA].sum()/molISOP*s2d*d2y[j]/10**12


### Compute Australian Anthropogenic Emissions ###

lon_mask_Australia = (lon >= 100) & (lon <= 160)
lat_mask_Australia = (lat >= -50) & (lat<= -10)
for i in range(len3):
    for j in range(len4):
        sum_omh[6,i,j] = (emiss_omh*surface)[i,j,lat_mask_Australia,:][:,lon_mask_Australia].sum()/molISOP*s2d*d2y[j]/10**12

### Northern Africa emissions ###
lon_mask_NAfrica = (lon >= -10) & (lon <= 50)
lat_mask_NAfrica= (lat >= 0) & (lat<= 30)
for i in range(len3):
    for j in range(len4):
        sum_omh[7,i,j] = (emiss_omh*surface)[i,j,lat_mask_NAfrica,:][:,lon_mask_NAfrica].sum()/molISOP*s2d*d2y[j]/10**12
    
### World ###
lon_mask_world = (lon >= -180) & (lon <= 180)
lat_mask_world = (lat >= -90) & (lat<= 90)
for i in range(len3):
    for j in range(len4):
        sum_omh[8,i,j] = (emiss_omh*surface)[i,j,lat_mask_world,:][:,lon_mask_world].sum()/molISOP*s2d*d2y[j]/10**12


fig = plt.figure(1,figsize=(16,7.5))
ax = plt.subplot(331)
x_series = arange(8)
plot(x_series,sum_omh[8,0,:],x_series,sum_omh[8,1,:],x_series,sum_omh[8,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP global emissions [Tg/month]")

ax = plt.subplot(332)
plot(x_series,sum_omh[0,0,:],x_series,sum_omh[0,1,:],x_series,sum_omh[0,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP Southern American emissions")
ax = plt.subplot(333)
plot(x_series,sum_omh[1,0,:],x_series,sum_omh[1,1,:],x_series,sum_omh[1,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP Southern African emissions")
ax = plt.subplot(334)
plot(x_series,sum_omh[2,0,:],x_series,sum_omh[2,1,:],x_series,sum_omh[2,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP South East Asian emissions")
ax = plt.subplot(335)
plot(x_series,sum_omh[3,0,:],x_series,sum_omh[3,1,:],x_series,sum_omh[3,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP East Asian emissions")
ax = plt.subplot(336)
plot(x_series,sum_omh[4,0,:],x_series,sum_omh[4,1,:],x_series,sum_omh[4,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP European emissions")
ax = plt.subplot(337)
plot(x_series,sum_omh[5,0,:],x_series,sum_omh[5,1,:],x_series,sum_omh[5,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP North American emissions")
ax = plt.subplot(338)
plot(x_series,sum_omh[6,0,:],x_series,sum_omh[6,1,:],x_series,sum_omh[6,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("ISOP Australian emissions")
ax = plt.subplot(339)
plot(x_series,sum_omh[7,0,:],x_series,sum_omh[7,1,:],x_series,sum_omh[7,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
legend(('CTRL','OMI only','All instrument'),loc='upper right')
title("ISOP Northern African emissions")


show()
