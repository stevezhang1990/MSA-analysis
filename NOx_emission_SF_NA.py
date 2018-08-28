#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
#infile_ap = Dataset('/users/jk/06/xzhang/gcadj_std_T_V34/runs/v8-02-01/geos5_osio3_1107/ctm.01.20110701.nc')
infile_ap = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080701.nc')
infile_ap2 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080704.nc')
infile_ap3 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080707.nc')
infile_ap4 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080710.nc')
infile_ap5 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080713.nc')
infile_ap6 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080716.nc')
infile_ap7 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080719.nc')
infile_ap8 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080722.nc')
infile_ap9 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_ap_0807_NA/ctm.00.20080725.nc')
infile_mop = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_1/gctm.sf.05.20080701.nc')
infile_mop2 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_2/gctm.sf.05.20080704.nc')
infile_mop3 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_3/gctm.sf.05.20080707.nc')
infile_mop4 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_4/gctm.sf.05.20080710.nc')
infile_mop5 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_5/gctm.sf.05.20080713.nc')
infile_mop6 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_6/gctm.sf.05.20080716.nc')
infile_mop7 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_7/gctm.sf.05.20080719.nc')
infile_mop8 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_8/gctm.sf.05.20080722.nc')
infile_mop9 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omo_0807_NA_9/gctm.sf.05.20080725.nc')
infile_omi = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_1/gctm.sf.05.20080701.nc')
infile_omi2 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_2/gctm.sf.05.20080704.nc')
infile_omi3 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_3/gctm.sf.05.20080707.nc')
infile_omi4 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_4/gctm.sf.05.20080710.nc')
infile_omi5 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_5/gctm.sf.05.20080713.nc')
infile_omi6 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_6/gctm.sf.05.20080716.nc')
infile_omi7 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_7/gctm.sf.05.20080719.nc')
infile_omi8 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_8/gctm.sf.05.20080722.nc')
infile_omi9 = Dataset('/users/jk/16/xzhang/gcadj_std_O_NA/runs/v8-02-01/geos5_omi_0807_NA_9/gctm.sf.05.20080725.nc')
#CHEM_L_S__XX or CHEM_L_S__XX
lon = infile_ap.variables["LON"][:]
lat = infile_ap.variables["LAT"][:]
ps = infile_ap.variables["DAO-FLDS__PS-PBL"][:,:]
NOx_an = np.zeros((9,len(lat),len(lon)))
NOx_bb = np.zeros((9,len(lat),len(lon)))
NOx_mop_an = np.zeros((9,len(lat),len(lon)))
NOx_mop_bb = np.zeros((9,len(lat),len(lon)))
NOx_omi_an = np.zeros((9,len(lat),len(lon)))
NOx_omi_bb = np.zeros((9,len(lat),len(lon)))
NOx_an[0,:,:] = infile_ap.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[0,:,:] = infile_ap.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[1,:,:] = infile_ap2.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[1,:,:] = infile_ap2.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[2,:,:] = infile_ap3.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[2,:,:] = infile_ap3.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[3,:,:] = infile_ap4.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[3,:,:] = infile_ap4.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[4,:,:] = infile_ap5.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[4,:,:] = infile_ap5.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[5,:,:] = infile_ap6.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[5,:,:] = infile_ap6.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[6,:,:] = infile_ap7.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[6,:,:] = infile_ap7.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[7,:,:] = infile_ap8.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[7,:,:] = infile_ap8.variables["NOX-BIOB__NOx"][0,:,:]
NOx_an[8,:,:] = infile_ap9.variables["ANTHSRCE__NOx"][0,:,:]
NOx_bb[8,:,:] = infile_ap9.variables["NOX-BIOB__NOx"][0,:,:]

NOx_mop_an[0,:,:] = infile_mop.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[0,:,:] = infile_mop.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[1,:,:] = infile_mop2.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[1,:,:] = infile_mop2.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[2,:,:] = infile_mop3.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[2,:,:] = infile_mop3.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[3,:,:] = infile_mop4.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[3,:,:] = infile_mop4.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[4,:,:] = infile_mop5.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[4,:,:] = infile_mop5.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[5,:,:] = infile_mop6.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[5,:,:] = infile_mop6.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[6,:,:] = infile_mop7.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[6,:,:] = infile_mop7.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[7,:,:] = infile_mop8.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[7,:,:] = infile_mop8.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_mop_an[8,:,:] = infile_mop9.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_mop_bb[8,:,:] = infile_mop9.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[0,:,:] = infile_omi.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[0,:,:] = infile_omi.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[1,:,:] = infile_omi2.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[1,:,:] = infile_omi2.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[2,:,:] = infile_omi3.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[2,:,:] = infile_omi3.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[3,:,:] = infile_omi4.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[3,:,:] = infile_omi4.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[4,:,:] = infile_omi5.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[4,:,:] = infile_omi5.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[5,:,:] = infile_omi6.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[5,:,:] = infile_omi6.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[6,:,:] = infile_omi7.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[6,:,:] = infile_omi7.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[7,:,:] = infile_omi8.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[7,:,:] = infile_omi8.variables["IJ-EMS-S__NOX_bb"][:,:]
NOx_omi_an[8,:,:] = infile_omi9.variables["IJ-EMS-S__NOX_an"][:,:]
NOx_omi_bb[8,:,:] = infile_omi9.variables["IJ-EMS-S__NOX_bb"][:,:]


NOx_ap = np.mean(NOx_an +NOx_bb,axis=0)
NOx_mopitt = (np.mean(NOx_mop_an,axis=0)-1)*1.5+1+(np.mean(NOx_mop_bb,axis=0)-1)*1.5#-0.40000002
NOx_omi = (np.mean(NOx_omi_an,axis=0)-1)*3.5+1+(np.mean(NOx_omi_bb,axis=0)-1)*3.5
#for i in range(len(lat)):
    #for j in range(len(lon)):
        #if NOx_mopitt[i,j] > 1.4:
            #if NOx_mopitt[i,j] < 1.61:
                #NOx_mopitt[i,j] = 1
#NO2_tes2 = infile_tes2.variables["CHEM_L_S__NO2"][:,:,:]
fig = plt.figure(1,figsize=(16,4.5))
ax = plt.subplot(131)
title('NOx emission estimates: CTRL run [1e9 molec/cm3/s]')
plot_data, lon = addcyclic(NOx_ap[:,:],lon)
m1=Basemap(lon_0=-90,llcrnrlat=10,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-40)
m1.drawcoastlines()
m1.drawparallels(arange(0,90,10),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,0,30),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ozone March 0-100 Nov 10-80
anth_plot = m1.pcolormesh(x,y,plot_data/1e9,vmin=0,vmax=200,cmap=get_cmap("gnuplot2_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
ax = plt.subplot(132)
title('NOx emission estimates: OMO [1e9 molec/cm3/s]')
plot_data = NOx_mopitt[:,:]
m1=Basemap(lon_0=-90,llcrnrlat=10,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-40)
m1.drawcoastlines()
m1.drawparallels(arange(0,90,10),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,0,30),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ozone March 0-100 Nov 10-80
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0.7,vmax=1.3,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
ax = plt.subplot(133)
title('NOx emission estimates: OMI [1e9 molec/cm3/s]')
plot_data = NOx_omi[:,:]
m1=Basemap(lon_0=-90,llcrnrlat=10,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-40)
m1.drawcoastlines()
m1.drawparallels(arange(0,90,10),labels=[1,0,0,1],labelstyle="+/-")
m1.drawmeridians(arange(-180,0,30),labels=[1,0,0,1],labelstyle="+/-")
x,y = m1(*np.meshgrid(lon,lat))
#For ozone March 0-100 Nov 10-80
anth_plot = m1.pcolormesh(x,y,plot_data,vmin=0.7,vmax=1.3,cmap=get_cmap("RdBu_r"))
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="7%", pad=0.3)
plt.colorbar(orientation = 'horizontal',cax=cax)
anth_plot.cmap.set_bad('lightgrey')
plt.show()
