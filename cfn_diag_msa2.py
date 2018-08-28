#!/usr/bin/python
from scipy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.colors import LogNorm
from geos_interpolate import geos_interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import numpy as np

path = '/users/jk/16/xzhang/gcadj_std_O_V35/runs/v8-02-01/geos5_omino2_omh_1608/diagadj/'
cfn_filename = ['mop','iasio3','osi','omino2','omich2o','mlshno3']
#IT_spec = ['01','05']
for i in range(6):
    for j in range(5):
        IT_spec = str(j+1).zfill(2)
        filename = path+'cfn_'+cfn_filename[i]+'.'+IT_spec+'.m'
        if os.path.isfile(filename):
            cfn_spec = np.sum(loadtxt(path+'cfn_'+cfn_filename[i]+'.'+IT_spec+'.m',dtype=float))
            print cfn_filename[i], IT_spec, cfn_spec
        
