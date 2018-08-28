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
#rc('text', usetex=True)
#rc('font',size=12)
infile00 = Dataset('/users/jk/16/xzhang/gcadj_std_M_V35/runs/v8-02-01/geos5_mop_1601/ctm.01.20160101.nc')
lon = infile00.variables["LON"][:]
lat = infile00.variables["LAT"][:]
ps = infile00.variables["DAO-FLDS__PS-PBL"][0,:,:]
surface = infile00.variables["DXYP__DXYP"][:,:] * (100**2) # cm^2
len1,len2,len3 = shape(infile00.variables["CHEM-L_S__OH"][:,:,:])
ems_mer = np.zeros((23,3,5,len2,len3))
spec_mer = np.zeros((23,2,2,len1,len2,len3))
sysdir = ['/users/jk/15/xzhang/','/users/jk/16/xzhang/']
codedir = ['gcadj_std_T_3d/','gcadj_std_I_3d/']
codedir2 = ['gcadj_std_M_V35','gcadj_std_M_V35','gcadj_std_M_SOB','gcadj_std_O_3d','gcadj_std_O_V35']
rundir2 = ['geos5_ap_2016/','geos5_mop_1608/','geos5_iasio3_1608_WC','geos5_omino2_1608/','geos5_omh_1608']
#T_3d: 00-15
rundir2w = ['geos5_all_1608_07/','geos5_all_1608_08/']
rundir2w2 = ['geos5_all_1608_05/','geos5_all_1608_06/']
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
path = ["" for x in range(23)]
path[0] = sysdir[1]+codedir2[0]+runbuffer+rundir2[0]
if len(glob.glob(os.path.join(path[0],filedir[0]+'.201608*.nc')))==0:
    filename=glob.glob(os.path.join(path[0],filedir2[0]+'.201608*.nc'))[0]
else:
    filename=glob.glob(os.path.join(path[0],filedir[0]+'.201608*.nc'))[0]
infile=Dataset(filename)
ems_mer[0,0,0,:,:] = infile.variables["CO--SRCE__COanth"][0,:,:]
ems_mer[0,0,1,:,:] = infile.variables["CO--SRCE__CObb"][0,:,:]
ems_mer[0,0,2,:,:] = infile.variables["ANTHSRCE__NOx"][0,:,:]
ems_mer[0,0,3,:,:] = infile.variables["NOX-BIOB__NOx"][0,:,:]
ems_mer[0,0,4,:,:] = infile.variables["BIOGSRCE__ISOP"][0,:,:]
spec_mer[0,0,0,:,:,:] = infile.variables["IJ-AVG-S__Ox"][:len1,:,:]
spec_mer[0,0,1,:,:,:] = infile.variables["CHEM-L_S__OH"][:,:,:]

for i in range(4):
    path_sub = sysdir[1]+codedir2[i+1]+runbuffer
    path[1+i] = path_sub+rundir2[i+1]
for i in range(2):
    if i==0:
        path_sub = sysdir[0]+codedir[i]+runbuffer
        for j in range(16):
            run_ind = str(j).zfill(2)
            path[5+j] = path_sub+'geos5_all_1608_'+run_ind
    else:
        path_sub = sysdir[1]+codedir[i]+runbuffer
        path[20+i] = path_sub+rundir2w[i]
        path[21+i] = path_sub+rundir2w2[i]

for j in range(22):
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
ems_exp = np.ma.masked_where(ems_mer<=0,ems_mer)
spec_exp = np.ma.masked_where(spec_mer<=0,spec_mer)

ems_exp_tot = np.zeros((23,3,len2,len3))
ems_exp_tot[0,0,:,:] = ems_exp[0,0,0,:,:] + ems_exp[0,0,1,:,:]
ems_exp_tot[0,1,:,:] = ems_exp[0,0,2,:,:] + ems_exp[0,0,3,:,:]
ems_exp_tot[0,2,:,:] = ems_exp[0,0,4,:,:]
ems_exp_tot[1:,0,:,:] = ems_exp[1:,0,0,:,:]*ems_exp[1:,2,0,:,:]+ ems_exp[1:,0,1,:,:]*ems_exp[1:,2,1,:,:]
ems_exp_tot[1:,1,:,:] = ems_exp[1:,0,2,:,:]*ems_exp[1:,2,2,:,:]+ ems_exp[1:,0,3,:,:]*ems_exp[1:,2,3,:,:]
ems_exp_tot[1:,2,:,:] = ems_exp[1:,0,4,:,:]*ems_exp[1:,2,4,:,:]
fig = plt.figure(1,figsize=(12,10))
x_series = arange(23)
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
sum_ems = np.zeros((23,5))
sum_ems[0,3] = trop_mwmean(len1,lat,lon,spec_exp[0,0,0,:,:,:],ps,spec='O3')
sum_ems[0,4] = trop_mwmean(len1,lat,lon,spec_exp[0,0,1,:,:,:],ps,spec='OH')
for i in range(22):
    sum_ems[i+1,3] = trop_mwmean(len1,lat,lon,spec_exp[i+1,1,0,:,:,:],ps,spec='O3')
    sum_ems[i+1,4] = trop_mwmean(len1,lat,lon,spec_exp[i+1,1,1,:,:,:],ps,spec='OH')
sum_ems[9,3] = sum_ems[1,3]
for m in range(5):
    ax=plt.subplot(511+m)
    k = 8
    lon_mask = (lon >= lon_left[k]) & (lon <= lon_right[k])
    lat_mask = (lat >= lat_low[k]) & (lat<= lat_high[k])
    if m<3:
        for i in range(23):
            sum_ems[i,m] = (ems_exp_tot*surface)[i,m,lat_mask,:][:,lon_mask].sum()/spec_mol[m]*s2d/10**12
    sum_ems[2,:3] = sum_ems[0,:3]
    plot(x_series,sum_ems[:,m],'ro')
    title('Total '+title1[k]+' '+title2[m]+' emissions/conc',fontsize=11)
    if m==4:
        xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],\
                   ['AP','MOP','IAS','OMI','OMH','T00','T01','T02','T03','T04','T05','T06','T07','T08',
                    'T09','T10','T11','T12','T13','T14','T15','I08','I06'])
sum_ems_bar = np.zeros(6)
sum_ems_bar[0] = sum_ems[0,0]
sum_ems_bar[1] = sum_ems[1,0]
sum_ems_bar[2] = sum_ems[16,0]
sum_ems_bar[3] = sum_ems[21,0]
sum_ems_bar[4] = sum_ems[12,0]
sum_ems_bar[5] = sum_ems[19,0]
x2 = np.arange(6)
fig = plt.figure(2,figsize=(8,6))
plt.bar(x2, sum_ems_bar, width=0.35)
xticks([0,1,2,3,4,5],['AP','MOP',r'$\gamma$=1',r'$\gamma$=4',r'$\gamma$=16',r'$\gamma$=64'])
show()
