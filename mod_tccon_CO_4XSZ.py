import os
import numpy as np

AP = [0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
    1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
    4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
    7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
    1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
    1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
    2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
    2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
    1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
    7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01,
    1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00,
    6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02]

BP = [1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
    9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
    8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
    7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
    6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
    4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
    2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
    6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00]

def get_intmap(model_vertical_levels, p_center_on_site, gc_psurf, tccon_vertical_levels, tccon_p_center, obs_ps):
   HINTERPZ= np.ndarray(shape=(model_vertical_levels, tccon_vertical_levels), dtype=np.float64, order='C')
   HINTERPZ[:,:]= 0.0

   DELTA_SURFP   = 0.5 * ( obs_ps -gc_psurf )
   for LGC in range(model_vertical_levels):
      p_center_on_site[LGC] = ( BP[LGC] + BP[LGC+1])* DELTA_SURFP + p_center_on_site[LGC]
   for LTM in range(tccon_vertical_levels):
      for LGC in range(model_vertical_levels-1):
         LOW = p_center_on_site[LGC+1]
         HI = p_center_on_site[LGC]
         if (tccon_p_center[LTM]<= HI) and (tccon_p_center[LTM]>LOW):
            DIFF = HI - LOW
            HINTERPZ[LGC+1,LTM] = (HI-tccon_p_center[LTM])/DIFF
            HINTERPZ[LGC,LTM] = (tccon_p_center[LTM]-LOW)/DIFF
   for LTM in range(tccon_vertical_levels):
      if tccon_p_center[LTM] > p_center_on_site[0]:
         HINTERPZ[0,LTM] = 1.0
         HINTERPZ[1:model_vertical_levels,LTM] = 0.0
   return(HINTERPZ)

model_vertical_levels=47
model_lon=72
model_lat=46
model_hour=24
tccon_vertical_levels=71
obs_max=200
month_counter=(31, 28, 31,30,31,30,31,31, 30,31,30,31)
ints=('00', '01', '02', '03', '04', '05', '06', '07', '08', '09','10',
      '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
      '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31')
model_co=np.ndarray(shape=(model_hour, model_vertical_levels, model_lat, model_lon), dtype=np.float64, order='C')
model_p_center=np.ndarray(shape=(model_hour, model_vertical_levels, model_lat,model_lon), dtype=np.float64, order='C')
model_p_edge=np.ndarray(shape=(model_hour, model_vertical_levels, model_lat,model_lon), dtype=np.float64, order='C') #the top edge is not included
pwt=np.ndarray(shape=(obs_max,tccon_vertical_levels), dtype=np.float64, order='C')
tccon_ak_ts=np.ndarray(shape=(obs_max,tccon_vertical_levels), dtype=np.float64, order='C')
tccon_co_ap=np.ndarray(shape=(obs_max,tccon_vertical_levels), dtype=np.float64, order='C')
tccon_pressure=np.ndarray(shape=(obs_max,tccon_vertical_levels), dtype=np.float64, order='C')
lev=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
tccon_hour=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
tccon_minute=np.ndarray(shape=(obs_max), dtype=np.int, order='C')
XCO=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
tccon_lon=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
tccon_lat=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
station_id=np.ndarray(shape=(obs_max), dtype='|S2', order='C')
tccon_psurf=np.ndarray(shape=(obs_max), dtype=np.float64, order='C')
GC_CO_NATIVE = np.ndarray(shape=(model_vertical_levels), dtype=np.float64, order='C')
CO_ap = np.ndarray(shape=(tccon_vertical_levels), dtype=np.float64, order='C')
GC_CO = np.ndarray(shape=(tccon_vertical_levels), dtype=np.float64, order='C')
CO_pert = np.ndarray(shape=(tccon_vertical_levels), dtype=np.float64, order='C')
p_center_on_site= np.ndarray(shape=(model_vertical_levels), dtype=np.float64, order='C')
tccon_p_center= np.ndarray(shape=(tccon_vertical_levels), dtype=np.float64, order='C')
MAP = np.ndarray(shape=(model_vertical_levels, tccon_vertical_levels), dtype=np.float64, order='C')
outfilen='tccon_mod_obs.txt'
fo=open(outfilen, 'a')
#SPECIFY YEAR, MONTH, AND DAY
for iy in [2010]:
   for im in range(12):
      for id in range(month_counter[im]):
         #define iy, im correctly, or you need check the availability of this file using isfile 
         mod_file =''  #file name of 3d fld of CO and edge pressure
         model_co= #CO MR in range of model_hour, model_vertical_levels, model_lat, model_lon
         model_p_edge= #Edge presssure from surf to bottom of 47 layer
         # READ TCCON FILE
         FILE_NAME='/users/jk/15/tccondataforxyz/tcco_xco_'+str(iy)+ints[im+1]+ints[id+1]+'.dat'
         file_by_day=0
         if os.path.isfile(FILE_NAME):
            with open(FILE_NAME) as f:
               while True:
                  line=f.readline()
                  if ("" == line):

                     f.close()
                     break;
                  list1=line.split()
                  tccon_hour[file_by_day]=int(list1[2])
                  tccon_minute[file_by_day]=int(list1[3])
                  XCO[file_by_day]=float(list1[4])*1.0e+9
                  tccon_lon[file_by_day]=float(list1[5])
                  tccon_lat[file_by_day]=float(list1[6])
                  tccon_psurf[file_by_day]=float(list1[11])
                  station_id[file_by_day]=list1[12]
                  line=f.readline()
                  for obs_lvl_counter0 in range(71):
                     pwt[file_by_day,obs_lvl_counter0]=float(line.split()[obs_lvl_counter0])
                  line=f.readline()
                  for obs_lvl_counter0 in range(71):
                     tccon_ak_ts[file_by_day,obs_lvl_counter0]=float(line.split()[obs_lvl_counter0])
                  line=f.readline()
                  for obs_lvl_counter0 in range(71):
                     tccon_co_ap[file_by_day,obs_lvl_counter0]=float(line.split()[obs_lvl_counter0])*1.0e+9
                  line=f.readline()
                  for obs_lvl_counter0 in range(71):
                     tccon_pressure[file_by_day,obs_lvl_counter0]=float(line.split()[obs_lvl_counter0])
                  file_by_day=file_by_day+1
            # LOCATE TCCON LOCATIONS IN GC MODEL
            for file_counter in range(file_by_day):
               lat_gc_projection=tccon_lat[file_counter]+88
               if lat_gc_projection<=0:
                  model_lat_index=0
               else:
                  model_lat_index=int(lat_gc_projection/4)+1
               lon_gc_projection=tccon_lon[file_counter]+182.5
               if lon_gc_projection>=360:
                  model_lon_index=0
               else:
                  model_lon_index=int(lon_gc_projection/5)
            # LOCATE TCCON TIME IN GC MODEL   
               hour=int(tccon_hour[file_counter])
            # LOCATE TCCON OBS: CO MIXING RATIO AND PRESSURE IN GC MODEL   
               for gc_lvl_counter in range(model_vertical_levels):
                  GC_CO_NATIVE[gc_lvl_counter]=model_co[hour,gc_lvl_counter,model_lat_index,model_lon_index]
                  if gc_lvl_counter<(model_vertical_levels-1):
                     p_center_on_site[gc_lvl_counter]=0.5*model_p_edge[hour,gc_lvl_counter,model_lat_index,model_lon_index] \
                     +0.5*model_p_edge[hour,gc_lvl_counter+1,model_lat_index,model_lon_index]
                  else:
                     p_center_on_site[gc_lvl_counter]=0.5*model_p_edge[hour,gc_lvl_counter, model_lat_index, model_lon_index]+0.5*0.01


               for obs_lvl_counter0 in range(tccon_vertical_levels):
		            tccon_p_center[obs_lvl_counter0] = tccon_pressure[file_counter, obs_lvl_counter0]
             # LOCATE GC SURFACE PRESSURE  
               gc_psurf =model_p_edge[hour,0, model_lat_index, model_lon_index]
             # MAPPING FROM GC LVLS TO TCCON LVLS, NEED #GC_LVLS,PCENTER ON TCCON GRID, GC_PSURF, #TCCON LVLS, TCCON PRESSURE, TCCON PSURF
               MAP[:, :]=get_intmap(model_vertical_levels,  p_center_on_site, gc_psurf, tccon_vertical_levels, tccon_p_center, tccon_psurf[file_counter])
             # COMPUTE from GC_CO_NATIVE TO GC_CO  
               for obs_lvl_counter in range(tccon_vertical_levels):
                  GC_CO[obs_lvl_counter] = 0.0
                  for L in range(model_vertical_levels):
                     GC_CO[obs_lvl_counter] = GC_CO[obs_lvl_counter] + MAP[L,obs_lvl_counter]*GC_CO_NATIVE[L]
               
               for L in range (tccon_vertical_levels):
                  CO_ap[L]=tccon_co_ap[file_counter, L]*1.0672
                  CO_pert[L] = GC_CO[L] - CO_ap[L]
               
               CO_HAT=0.0
               for obs_lvl_counter in range(tccon_vertical_levels):
                  CO_HAT=CO_HAT+pwt[file_counter,obs_lvl_counter]*CO_ap[obs_lvl_counter]\
                  +pwt[file_counter,obs_lvl_counter]*tccon_ak_ts[file_counter,obs_lvl_counter]*CO_pert[obs_lvl_counter]
                  
               fo.write('%d %d  %d  %d  %d %f  %f %f %f %s\n'%( iy, im+1, id+1, tccon_hour[file_counter], tccon_minute[file_counter], \
                                                                tccon_lat[file_counter],tccon_lon[file_counter], XCO[file_counter], CO_HAT, station_id[file_counter]))
