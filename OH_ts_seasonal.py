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
mopfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.10.20160101.nc')
mopfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601_2/ctm.10.20160116.nc')
mopfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602/ctm.10.20160201.nc')
mopfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1602_2/ctm.10.20160216.nc')
mopfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603/ctm.10.20160301.nc')
mopfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1603_2/ctm.10.20160316.nc')
mopfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604/ctm.10.20160401.nc')
mopfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1604_2/ctm.10.20160416.nc')
allfile01 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_06/ctm.20.20160101.nc')
allfile21 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1601_26/ctm.15.20160116.nc')
allfile02 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_06/ctm.20.20160201.nc')
allfile22 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1602_26/ctm.20.20160216.nc')
allfile03 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_06/ctm.15.20160301.nc')
allfile23 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1603_26/ctm.20.20160316.nc')
allfile04 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_06/ctm.15.20160401.nc')
allfile24 = Dataset('/users/jk/16/xzhang/gcadj_std_I_3d/runs/v8-02-01/geos5_all_1604_26/ctm.15.20160416.nc')
infile_am = Dataset('/users/jk/15/xzhang/gcadj_std_T_3d/runs/v8-02-01/geos5_all_1608_AM/ctm.00.20160801.nc')
spi = np.array([11.265,11.265,10.91,10.91,10.99,10.99,10.7,10.7])
transcom = np.array([9.15,9.15,9.02,9.02,9.11,9.11,8.95,8.95])

#XOX-BIOB__XOx or anTHSRCE__XOx or XOX-an-S__XOx or xOX-LI-S__xOx
lon = infile01.variables["LON"][:]
lat = infile01.variables["LAT"][:]
ps = infile01.variables["DAO-FLDS__PS-PBL"][:,:]
surface = infile01.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len5 = shape(mopfile01.variables["CHEM-L_S__OH"][:,:,:])
len3 = 3
len4 = 8
OH = np.zeros((len3,len4,len1,len2,len5))
AM = infile_am.variables["BXHGHT-S__AD"][0:38,:,:]
OH[0,0,:,:,:] = infile01.variables["CHEM-L_S__OH"][:,:,:]
OH[0,1,:,:,:] = infile21.variables["CHEM-L_S__OH"][:,:,:]
OH[0,2,:,:,:] = infile02.variables["CHEM-L_S__OH"][:,:,:]
OH[0,3,:,:,:] = infile21.variables["CHEM-L_S__OH"][:,:,:]
OH[0,4,:,:,:] = infile03.variables["CHEM-L_S__OH"][:,:,:]
OH[0,5,:,:,:] = infile23.variables["CHEM-L_S__OH"][:,:,:]
OH[0,6,:,:,:] = infile04.variables["CHEM-L_S__OH"][:,:,:]
OH[0,7,:,:,:] = infile24.variables["CHEM-L_S__OH"][:,:,:]

OH[1,0,:,:,:] = mopfile01.variables["CHEM-L_S__OH"][:,:,:]
OH[1,1,:,:,:] = mopfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[1,2,:,:,:] = mopfile02.variables["CHEM-L_S__OH"][:,:,:]
OH[1,3,:,:,:] = mopfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[1,4,:,:,:] = mopfile03.variables["CHEM-L_S__OH"][:,:,:]
OH[1,5,:,:,:] = mopfile23.variables["CHEM-L_S__OH"][:,:,:]
OH[1,6,:,:,:] = mopfile04.variables["CHEM-L_S__OH"][:,:,:]
OH[1,7,:,:,:] = mopfile24.variables["CHEM-L_S__OH"][:,:,:]


OH[2,0,:,:,:] = allfile01.variables["CHEM-L_S__OH"][:,:,:]
OH[2,1,:,:,:] = allfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[2,2,:,:,:] = allfile02.variables["CHEM-L_S__OH"][:,:,:]
OH[2,3,:,:,:] = allfile21.variables["CHEM-L_S__OH"][:,:,:]
OH[2,4,:,:,:] = allfile03.variables["CHEM-L_S__OH"][:,:,:]
OH[2,5,:,:,:] = allfile23.variables["CHEM-L_S__OH"][:,:,:]
OH[2,6,:,:,:] = allfile04.variables["CHEM-L_S__OH"][:,:,:]
OH[2,7,:,:,:] = allfile24.variables["CHEM-L_S__OH"][:,:,:]

OH_avg = np.zeros((3,8))

for i in range(len3):
    for j in range(len4):
        OH_avg[i,j] = np.sum(OH[i,j,:,:,:]*AM)/(np.sum(AM)*1e5)
        
fig = plt.figure(1,figsize=(16,7.5))
x_series = arange(8)
plot(x_series,OH_avg[0,:],x_series,OH_avg[1,:],x_series,OH_avg[2,:],x_series,spi,linewidth=2.0)
xticks([0,2,4,6],['JAN','FEB','MAR','APR'])
title("Global mean OH concentrations (1e5 molec/cm3) Jan-April 2016")
legend(('CTRL','MOPITT','All instrument','Spivakovsky et al 2000'),loc='best')
ylim(10,14)

show()
