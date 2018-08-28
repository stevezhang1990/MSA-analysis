#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate,trop_mwmean
import os
import glob
from netCDF4 import Dataset

def draw_screen_mask( lat_mask, lon_mask, m):

    lats = [lat_mask[0],lat_mask[-1],lat_mask[-1],lat_mask[0]]
    lons = [lon_mask[0],lon_mask[0],lon_mask[-1],lon_mask[-1]]
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor='blue', alpha=0.4 )
    gca().add_patch(poly)
# get lat,lon,ps, surface area from sample file
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
ps = infile00.variables["DAO-FLDS__PS-PBL"][0,:,:]
surface = infile00.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len3 = shape(infile00.variables["CHEM-L_S__OH"][:,:,:])
# list of strings for different file entries
ems_mer = np.zeros((26,3,5,len2,len3))
spec_mer = np.zeros((26,2,2,len1,len2,len3))
sysdir = ['/users/jk/15/xzhang/','/users/jk/16/xzhang/']
codedir = ['gcadj_std_T_3d/','gcadj_std_I_3d/','gcadj_std_O_V35/']
codedir2 = ['gcadj_std_M_V35','gcadj_std_M_V35','gcadj_std_M_SOB','gcadj_std_O_3d','gcadj_std_O_V35']
rundir2 = ['geos5_ap_2016/','geos5_mop_1608/','geos5_iasio3_1608_WC','geos5_omino2_1608/','geos5_omh_1608']
rundir2w = ['geos5_all_1608_07/','geos5_all_1608_08/','geos5_all_1608_03/']
rundir1w = ['geos5_all_1608_7d_02/','geos5_all_1608_7d_00/','geos5_all_1608_7d_51/']
rundir1w2 = ['geos5_all_1608_7d_22/','geos5_all_1608_7d_20/','geos5_all_1608_7d_61/']
rundir3d = ['geos5_all_1608_02/','geos5_all_1608_09/','geos5_all_1608_08/']
rundir3d2 = ['01-04/','04-08/','08-11/','11-15/']
runbuffer = '/runs/v8-02-01/'
filedir = ['ctm.01','ctm.20','gctm.sf.20']
filedir2 = ['ctm.00','ctm.10','gctm.sf.10']
filedir3 = ['ctm.00','ctm.15','gctm.sf.15']
mol_mass = np.array([28,14,68.12])
spec_mol = 6.022*10**23/mol_mass
s2d = 60*60*24*31
# there are 26 sub experiments in total
path = ["" for x in range(26)]
path[0] = sysdir[1]+codedir2[0]+runbuffer+rundir2[0]
if len(glob.glob(os.path.join(path[0],filedir[0]+'.201608*.nc')))==0:
    filename=glob.glob(os.path.join(path[0],filedir2[0]+'.201608*.nc'))[0]
else:
    filename=glob.glob(os.path.join(path[0],filedir[0]+'.201608*.nc'))[0]
infile=Dataset(filename)
# read in modeled surface CO, NOx and isoprene emissions, and O3 and OH abundances for CTRL run
ems_mer[0,0,0,:,:] = infile.variables["CO--SRCE__COanth"][0,:,:]
ems_mer[0,0,1,:,:] = infile.variables["CO--SRCE__CObb"][0,:,:]
ems_mer[0,0,2,:,:] = infile.variables["ANTHSRCE__NOx"][0,:,:]
ems_mer[0,0,3,:,:] = infile.variables["NOX-BIOB__NOx"][0,:,:]
ems_mer[0,0,4,:,:] = infile.variables["BIOGSRCE__ISOP"][0,:,:]
spec_mer[0,0,0,:,:,:] = infile.variables["IJ-AVG-S__Ox"][:len1,:,:]
spec_mer[0,0,1,:,:,:] = infile.variables["CHEM-L_S__OH"][:,:,:]
# different combinations of subentries
for i in range(4):
    path_sub = sysdir[1]+codedir2[i+1]+runbuffer
    path[1+i] = path_sub+rundir2[i+1]
for i in range(3):
    if i==0:
        path_sub = sysdir[0]+codedir[i]+runbuffer
    else:
        path_sub = sysdir[1]+codedir[i]+runbuffer
    path[5+i] = path_sub+rundir2w[i]
    path[8+i] = path_sub+rundir1w[i]
    path[11+i] = path_sub+rundir1w2[i]
    for j2 in range(4):
        path[14+4*i+j2] = path_sub+rundir3d[i]+rundir3d2[j2]
# same as before, but for all the assimilation runs.
for j in range(25):
    for k in range(3):
        if len(glob.glob(os.path.join(path[j+1],filedir[k]+'.201608*.nc')))==0:
            if len(glob.glob(os.path.join(path[j+1],filedir2[k]+'.201608*.nc')))==0:
                filename=glob.glob(os.path.join(path[j+1],filedir3[k]+'.201608*.nc'))[0]
            else:
                filename=glob.glob(os.path.join(path[j+1],filedir2[k]+'.201608*.nc'))[0]
        else:
            filename=glob.glob(os.path.join(path[j+1],filedir[k]+'.201608*.nc'))[0]
        infile=Dataset(filename)
        #print filename
        if k<2:
            ems_mer[j+1,k,0,:,:] = infile.variables["CO--SRCE__COanth"][0,:,:]
            ems_mer[j+1,k,1,:,:] = infile.variables["CO--SRCE__CObb"][0,:,:]
            ems_mer[j+1,k,2,:,:] = infile.variables["ANTHSRCE__NOx"][0,:,:]
            ems_mer[j+1,k,3,:,:] = infile.variables["NOX-BIOB__NOx"][0,:,:]
            ems_mer[j+1,k,4,:,:] = infile.variables["BIOGSRCE__ISOP"][0,:,:]
            spec_mer[j+1,k,0,:,:,:] = infile.variables["IJ-AVG-S__Ox"][:len1,:,:]
            spec_mer[j+1,k,1,:,:,:] = infile.variables["CHEM-L_S__OH"][:,:,:]
        elif j<4 and j!=1:
            ems_mer[j+1,k,0,:,:] = infile.variables["IJ-EMS-S__CO_an"][:,:]
            ems_mer[j+1,k,1,:,:] = infile.variables["IJ-EMS-S__CO_bb"][:,:]
            ems_mer[j+1,k,2,:,:] = infile.variables["IJ-EMS-S__NOX_an"][:,:]
            ems_mer[j+1,k,3,:,:] = infile.variables["IJ-EMS-S__NOX_bb"][:,:]
            ems_mer[j+1,k,4,:,:] = infile.variables["IJ-EMS-S__ISOP_an"][:,:]
        elif j!=1:
            ems_mer[j+1,k,0,:,:] = infile.variables["IJ-EMS-S__CO_an"][0,:,:]
            ems_mer[j+1,k,1,:,:] = infile.variables["IJ-EMS-S__CO_bb"][0,:,:]
            ems_mer[j+1,k,2,:,:] = infile.variables["IJ-EMS-S__NOX_an"][0,:,:]
            ems_mer[j+1,k,3,:,:] = infile.variables["IJ-EMS-S__NOX_bb"][0,:,:]
            ems_mer[j+1,k,4,:,:] = infile.variables["IJ-EMS-S__ISOP_an"][0,:,:]

# masked all the invalid and bad data when calculating the mean
ems_exp = np.zeros((14,3,5,len2,len3))
spec_exp = np.zeros((14,2,2,len1,len2,len3))
ems_exp[0:8,:,:,:,:] = np.ma.masked_where(ems_mer[0:8,:,:,:,:]<=0,ems_mer[0:8,:,:,:,:])
spec_exp[0:8,:,:,:,:,:] = np.ma.masked_where(spec_mer[0:8,:,:,:,:,:]<=0,spec_mer[0:8,:,:,:,:,:])

for i in range(3):
    ems_exp[8+i,:,:,:,:] = np.mean(np.ma.masked_where(ems_mer[[8+i,11+i],:,:,:,:]<=0,\
                                                        ems_mer[[8+i,11+i],:,:,:,:]),axis=0)
    ems_exp[11+i,:,:,:,:] = np.mean(np.ma.masked_where(ems_mer[[14+i,17+i,20+i,23+i],:,:,:,:]<=0,\
                                                        ems_mer[[14+i,17+i,20+i,23+i],:,:,:,:]),axis=0)
    spec_exp[8+i,:,:,:,:,:] = np.mean(np.ma.masked_where(spec_mer[[8+i,11+i],:,:,:,:,:]<=0,\
                                                           spec_mer[[8+i,11+i],:,:,:,:,:]),axis=0)
    spec_exp[11+i,:,:,:,:,:] = np.mean(np.ma.masked_where(spec_mer[[14+i,17+i,20+i,23+i],:,:,:,:,:]<=0,\
                                                           spec_mer[[14+i,17+i,20+i,23+i],:,:,:,:,:]),axis=0)
# compute total emissions including anth and bb emissions with their corresponding SF.
ems_exp_tot = np.zeros((14,3,len2,len3))
ems_exp_tot[0,0,:,:] = ems_exp[0,0,0,:,:] + ems_exp[0,0,1,:,:]
ems_exp_tot[0,1,:,:] = ems_exp[0,0,2,:,:] + ems_exp[0,0,3,:,:]
ems_exp_tot[0,2,:,:] = ems_exp[0,0,4,:,:]
ems_exp_tot[1:,0,:,:] = ems_exp[1:,0,0,:,:]*ems_exp[1:,2,0,:,:]+ ems_exp[1:,0,1,:,:]*ems_exp[1:,2,1,:,:]
ems_exp_tot[1:,1,:,:] = ems_exp[1:,0,2,:,:]*ems_exp[1:,2,2,:,:]+ ems_exp[1:,0,3,:,:]*ems_exp[1:,2,3,:,:]
ems_exp_tot[1:,2,:,:] = ems_exp[1:,0,4,:,:]*ems_exp[1:,2,4,:,:]
fig = plt.figure(1,figsize=(8,10))
x_series = arange(14)
title1 = ['Southern American','Southern African','Indian and Southeast Asian','East Asian',\
          'Boreal Euroasian','North American','Australian','Northern African','Global']
#lat: SE Asia only: -10,10; SE Asia + India: -10,20; North America: 10,60; US 48 states: 30,50
#lat: Boreal Euroasia: 55,80; Europe: 35,70
#lon: Boreal Euroasia: -30-180; Europe:-10-60
lon_left=np.array([-80,10,95,95,-30,-130,100,-10,-180])
lon_right = np.array([-40,50,160,150,180,-60,160,50,180])
lat_low = np.array([-50,-30,-10,20,55,10,-50,0,-90])
lat_high = np.array([10,0,20,50,80,80,-10,30,90])
#title2 = ['anth','bb','total']
title2 = ['CO','NOx','ISOP',r'$O_3$','OH']
sum_ems = np.zeros((14,5))
# compute tropospheric mean followed by Spivakovsky et al 2000 for O3 and OH
sum_ems[0,3] = trop_mwmean(len1,lat,lon,spec_exp[0,0,0,:,:,:],ps,spec='O3')
sum_ems[0,4] = trop_mwmean(len1,lat,lon,spec_exp[0,0,1,:,:,:],ps,spec='OH')
for i in range(13):
    sum_ems[i+1,3] = trop_mwmean(len1,lat,lon,spec_exp[i+1,1,0,:,:,:],ps,spec='O3')
    sum_ems[i+1,4] = trop_mwmean(len1,lat,lon,spec_exp[i+1,1,1,:,:,:],ps,spec='OH')
sum_ems[10,3] = sum_ems[1,3]
# compute regional and global emissions
for m in range(5):
    ax=plt.subplot(511+m)
    k = 8
    lon_mask = (lon >= lon_left[k]) & (lon <= lon_right[k])
    lat_mask = (lat >= lat_low[k]) & (lat<= lat_high[k])
    if m<3:
        for i in range(14):
            sum_ems[i,m] = (ems_exp_tot*surface)[i,m,lat_mask,:][:,lon_mask].sum()/spec_mol[m]*s2d/10**12
    sum_ems[2,:] = sum_ems[0,:]
    plot(x_series,sum_ems[:,m],'ro')
    title('Total '+title1[k]+' '+title2[m]+' emissions/conc',fontsize=11)
    if m==4:
        xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13],['AP','MOP','IAS','OMI','OMH','T2W','I2W','O2W','T1W','I1W','O1W','T3D','I3D','O3D'])
# compute bar plot on illustration of the sensitivity runs using data assimilation window. 
sum_ems_bar = np.zeros((5,5))
sum_ems_bar[0,:] = sum_ems[0,:]
sum_ems_bar[1,0] = sum_ems[1,0]
sum_ems_bar[1,1] = sum_ems[3,1]
sum_ems_bar[1,2] = sum_ems[4,2]
sum_ems_bar[1,3] = sum_ems[2,3]
sum_ems_bar[1,4] = sum_ems[1,4]
sum_ems_bar[2,:] = sum_ems[6,:]
sum_ems_bar[3,:] = sum_ems[8,:]
sum_ems_bar[4,:] = sum_ems[11,:]
x2 = np.arange(5)
fig =plt.figure(2,figsize=(8,10))
tick_info = ['MOP','OMI','OMI','IASI','MOP']
for j in range(5):
    ax=plt.subplot(511+j)
    plot(x2,sum_ems_bar[:,j],'bo')
    title('Total '+title2[j]+' emissions/conc',fontsize=11)
    xticks([0,1,2,3,4],['AP',tick_info[j],'2 Week','1 week','3 Day'])
show()
