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
molNOx = 6.022*10**23 / ( 14.0 )
#NOx is reported to be 14
s2d = 60*60*24
d2y = np.array([15,16,14,15,15,16,15,15,15,16,15,15,15,16,15,16,15,15,15,16,15,15,15,16])
infile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601/ctm.01.20160101.nc')
infile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601_2/ctm.01.20160116.nc')
infile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602/ctm.01.20160201.nc')
infile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602_2/ctm.01.20160216.nc')
infile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603/ctm.01.20160301.nc')
infile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603_2/ctm.01.20160316.nc')
infile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604/ctm.01.20160401.nc')
infile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604_2/ctm.01.20160416.nc')
omifile01 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601/gctm.sf.10.20160101.nc')
omifile21 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1601_2/gctm.sf.10.20160116.nc')
omifile02 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602/gctm.sf.10.20160201.nc')
omifile22 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1602_2/gctm.sf.10.20160216.nc')
omifile03 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603/gctm.sf.10.20160301.nc')
omifile23 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1603_2/gctm.sf.10.20160316.nc')
omifile04 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604/gctm.sf.10.20160401.nc')
omifile24 = Dataset('/users/jk/16/xzhang/gcadj_std_O_3d/runs/v8-02-01/geos5_omino2_1604_2/gctm.sf.10.20160416.nc')
allfile01 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1601_05/gctm.sf.20.20160101.nc')
allfile21 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1601_25/gctm.sf.15.20160116.nc')
allfile02 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1602_05/gctm.sf.15.20160201.nc')
allfile22 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1602_25/gctm.sf.20.20160216.nc')
allfile03 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1603_05/gctm.sf.15.20160301.nc')
allfile23 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1603_25/gctm.sf.20.20160316.nc')
allfile04 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1604_05/gctm.sf.20.20160401.nc')
allfile24 = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1604_25/gctm.sf.20.20160416.nc')

#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2 = shape(omifile01.variables["IJ-EMS-S__NOX_an"][:,:])
len3 = 3
len4 = 8
NOX_an = np.zeros((len4,len1,len2))
NOX_bb = np.zeros((len4,len1,len2))

SF_omi_an = np.zeros((len4,len1,len2))
SF_omi_bb = np.zeros((len4,len1,len2))

SF_all_an = np.zeros((len4,len1,len2))
SF_all_bb = np.zeros((len4,len1,len2))

NOX_an[0,:,:] = infile01.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[1,:,:] = infile21.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[2,:,:] = infile02.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[3,:,:] = infile21.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[4,:,:] = infile03.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[5,:,:] = infile23.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[6,:,:] = infile04.variables["ANTHSRCE__NOx"][0,:,:]
NOX_an[7,:,:] = infile24.variables["ANTHSRCE__NOx"][0,:,:]

NOX_bb[0,:,:] = infile01.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[1,:,:] = infile21.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[2,:,:] = infile02.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[3,:,:] = infile22.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[4,:,:] = infile03.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[5,:,:] = infile23.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[6,:,:] = infile04.variables["NOX-BIOB__NOx"][0,:,:]
NOX_bb[7,:,:] = infile24.variables["NOX-BIOB__NOx"][0,:,:]

SF_omi_an[0,:,:] = omifile01.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[1,:,:] = omifile21.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[2,:,:] = omifile02.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[3,:,:] = omifile21.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[4,:,:] = omifile03.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[5,:,:] = omifile23.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[6,:,:] = omifile04.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_an[7,:,:] = omifile24.variables["IJ-EMS-S__NOX_an"][:,:]
SF_omi_bb[0,:,:] = omifile01.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[1,:,:] = omifile21.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[2,:,:] = omifile02.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[3,:,:] = omifile22.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[4,:,:] = omifile03.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[5,:,:] = omifile23.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[6,:,:] = omifile04.variables["IJ-EMS-S__NOX_bb"][:,:]
SF_omi_bb[7,:,:] = omifile24.variables["IJ-EMS-S__NOX_bb"][:,:]

SF_all_an[0,:,:] = allfile01.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[1,:,:] = allfile21.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[2,:,:] = allfile02.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[3,:,:] = allfile21.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[4,:,:] = allfile03.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[5,:,:] = allfile23.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[6,:,:] = allfile04.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_an[7,:,:] = allfile24.variables["IJ-EMS-S__NOX_an"][0,:,:]
SF_all_bb[0,:,:] = allfile01.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[1,:,:] = allfile21.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[2,:,:] = allfile02.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[3,:,:] = allfile22.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[4,:,:] = allfile03.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[5,:,:] = allfile23.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[6,:,:] = allfile04.variables["IJ-EMS-S__NOX_bb"][0,:,:]
SF_all_bb[7,:,:] = allfile24.variables["IJ-EMS-S__NOX_bb"][0,:,:]

NOx_ap = NOX_an +NOX_bb
NOx_omi = SF_omi_an+SF_omi_bb-1
NOx_all = SF_all_an+SF_all_bb -1

emiss_moa = np.zeros((len3,len4,len1,len2))
emiss_mob = np.zeros((len3,len4,len1,len2))

emiss_moa[0,:,:,:] = NOX_an
emiss_moa[1,:,:,:] = SF_omi_an * NOX_an
emiss_moa[2,:,:,:] = SF_all_an * NOX_an

emiss_mob[0,:,:,:] = NOX_bb
emiss_mob[1,:,:,:] = SF_omi_bb * NOX_bb
emiss_mob[2,:,:,:] = SF_all_bb * NOX_bb

emiss_omi = emiss_moa + emiss_mob

sum_omi = np.zeros((9,len3,len4))

### Compute Southern America Biomass Burning Emissions ### /Or Northern Africa depends on the regions

lon_mask_SAmerica = (lon >= -80) & (lon <= -40)
lat_mask_SAmerica = (lat >= -50) & (lat<= 10)

for i in range(len3):
    for j in range(len4):
        sum_omi[0,i,j] = (emiss_omi*surface)[i,j,lat_mask_SAmerica,:][:,lon_mask_SAmerica].sum()/molNOx*s2d*d2y[j]/10**12

### Compute Southern Africa Biomass Burning Emissions ###

lon_mask_SAfrica = (lon >= 10) & (lon <= 50)
lat_mask_SAfrica = (lat >= -30) & (lat<= 0)
for i in range(len3):
    for j in range(len4):
        sum_omi[1,i,j] = (emiss_omi*surface)[i,j,lat_mask_SAfrica,:][:,lon_mask_SAfrica].sum()/molNOx*s2d*d2y[j]/10**12


### Compute South East Asian Anthropogenic Emissions ###

lon_mask_SE_Asia = (lon >= 95) & (lon <= 160)
lat_mask_SE_Asia = (lat >= -10) & (lat<= 10)

for i in range(len3):
    for j in range(len4):
        sum_omi[2,i,j] = (emiss_omi*surface)[i,j,lat_mask_SE_Asia,:][:,lon_mask_SE_Asia].sum()/molNOx*s2d*d2y[j]/10**12


### Compute East Asian Anthropogenic Emissions ###

lon_mask_Asia = (lon >= 95) & (lon <= 150)
lat_mask_Asia = (lat >= 20) & (lat<= 50)
for i in range(len3):
    for j in range(len4):
        sum_omi[3,i,j] = (emiss_omi*surface)[i,j,lat_mask_Asia,:][:,lon_mask_Asia].sum()/molNOx*s2d*d2y[j]/10**12

### Compute European Anthropogenic Emissions ###

lon_mask_Europe = (lon >= -10) & (lon <= 60)
lat_mask_Europe = (lat >= 35) & (lat<= 70)
for i in range(len3):
    for j in range(len4):
        sum_omi[4,i,j] = (emiss_omi*surface)[i,j,lat_mask_Europe,:][:,lon_mask_Europe].sum()/molNOx*s2d*d2y[j]/10**12

### Compute North American Anthropogenic Emissions ###

lon_mask_NA = (lon >= -130) & (lon <= -60)
lat_mask_NA = (lat >= 10) & (lat<= 70)
for i in range(len3):
    for j in range(len4):
        sum_omi[5,i,j] = (emiss_omi*surface)[i,j,lat_mask_NA,:][:,lon_mask_NA].sum()/molNOx*s2d*d2y[j]/10**12


### Compute Australian Anthropogenic Emissions ###

lon_mask_Australia = (lon >= 100) & (lon <= 160)
lat_mask_Australia = (lat >= -50) & (lat<= -10)
for i in range(len3):
    for j in range(len4):
        sum_omi[6,i,j] = (emiss_omi*surface)[i,j,lat_mask_Australia,:][:,lon_mask_Australia].sum()/molNOx*s2d*d2y[j]/10**12

### Northern Africa emissions ###
lon_mask_NAfrica = (lon >= -10) & (lon <= 50)
lat_mask_NAfrica= (lat >= 0) & (lat<= 30)
for i in range(len3):
    for j in range(len4):
        sum_omi[7,i,j] = (emiss_omi*surface)[i,j,lat_mask_NAfrica,:][:,lon_mask_NAfrica].sum()/molNOx*s2d*d2y[j]/10**12
    
### World ###
lon_mask_world = (lon >= -180) & (lon <= 180)
lat_mask_world = (lat >= -90) & (lat<= 90)
for i in range(len3):
    for j in range(len4):
        sum_omi[8,i,j] = (emiss_omi*surface)[i,j,lat_mask_world,:][:,lon_mask_world].sum()/molNOx*s2d*d2y[j]/10**12


fig = plt.figure(1,figsize=(16,7.5))
ax = plt.subplot(331)
x_series = arange(8)
plot(x_series,sum_omi[8,0,:],x_series,sum_omi[8,1,:],x_series,sum_omi[8,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx global emissions [Tg/month]")

ax = plt.subplot(332)
plot(x_series,sum_omi[0,0,:],x_series,sum_omi[0,1,:],x_series,sum_omi[0,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx Southern American emissions")
ax = plt.subplot(333)
plot(x_series,sum_omi[1,0,:],x_series,sum_omi[1,1,:],x_series,sum_omi[1,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx Southern African emissions")
ax = plt.subplot(334)
plot(x_series,sum_omi[2,0,:],x_series,sum_omi[2,1,:],x_series,sum_omi[2,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx South East Asian emissions")
ax = plt.subplot(335)
plot(x_series,sum_omi[3,0,:],x_series,sum_omi[3,1,:],x_series,sum_omi[3,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx East Asian emissions")
ax = plt.subplot(336)
plot(x_series,sum_omi[4,0,:],x_series,sum_omi[4,1,:],x_series,sum_omi[4,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx European emissions")
ax = plt.subplot(337)
plot(x_series,sum_omi[5,0,:],x_series,sum_omi[5,1,:],x_series,sum_omi[5,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx North American emissions")
ax = plt.subplot(338)
plot(x_series,sum_omi[6,0,:],x_series,sum_omi[6,1,:],x_series,sum_omi[6,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("NOx Australian emissions")
ax = plt.subplot(339)
plot(x_series,sum_omi[7,0,:],x_series,sum_omi[7,1,:],x_series,sum_omi[7,2,:],linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
legend(('CTRL','OMI only','All instrument'),loc='upper right')
title("NOx Northern African emissions")


show()
