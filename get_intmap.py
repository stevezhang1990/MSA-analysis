import os
import numpy as np
# AP, BP value for GEOS5, GEOSFP reduced vertical level grids
AP = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
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
 6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

BP = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
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
 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

# compute latitudinal mean of the modeled concentrations
def get_atom_chart(data_in):
    len_exp1,len_spec1,len_alt1,len_lon1 = shape(data_in)
    lvl_crt = int(len_alt1/2)
    data_in_chart = np.zeros((len_exp1,len_spec1,13))
    data_in_chart[:,:,0] = np.mean(np.mean(data_in[:,:,:,:],axis=2),axis=2)
    data_in_chart[:,:,1] = np.mean(np.mean(data_in[:,:,0:lvl_crt,:],axis=2),axis=2)
    data_in_chart[:,:,2] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,:],axis=2),axis=2)
    data_in_chart[:,:,3] = np.mean(np.mean(data_in[:,:,0:lvl_crt,150:180],axis=2),axis=2)
    data_in_chart[:,:,4] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,150:180],axis=2),axis=2)
    data_in_chart[:,:,5] = np.mean(np.mean(data_in[:,:,0:lvl_crt,120:150],axis=2),axis=2)
    data_in_chart[:,:,6] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,120:150],axis=2),axis=2)
    data_in_chart[:,:,7] = np.mean(np.mean(data_in[:,:,0:lvl_crt,60:120],axis=2),axis=2)
    data_in_chart[:,:,8] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,60:120],axis=2),axis=2)
    data_in_chart[:,:,9] = np.mean(np.mean(data_in[:,:,0:lvl_crt,30:60],axis=2),axis=2)
    data_in_chart[:,:,10] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,30:60],axis=2),axis=2)
    data_in_chart[:,:,11] = np.mean(np.mean(data_in[:,:,0:lvl_crt,0:30],axis=2),axis=2)
    data_in_chart[:,:,12] = np.mean(np.mean(data_in[:,:,lvl_crt:len_alt,0:30],axis=2),axis=2)
    return data_in_chart      

# convert model vertical grid to obs vertical grid
def get_intmap(mod_l, mod_pc, mod_ps, obs_l, obs_pc, obs_ps, mode='MR'):
   if mode == 'MR':
      HINTERPZ= np.ndarray(shape=(mod_l, obs_l), dtype=np.float64, order='C')
      HINTERPZ[:,:]= 0.0
      DELTA_SURFP   = 0.5 * ( obs_ps -mod_ps )
      for LGC in range(mod_l):
         mod_pc[LGC] = ( BP[LGC] + BP[LGC+1])* DELTA_SURFP + mod_pc[LGC]
      for LTM in range(obs_l):
         for LGC in range(mod_l-1):
            LOW = mod_pc[LGC+1]
            HI = mod_pc[LGC]
            if (obs_pc[LTM]<= HI) and (obs_pc[LTM]>LOW):
               DIFF = HI - LOW
               HINTERPZ[LGC+1,LTM] = (HI-obs_pc[LTM])/DIFF
               HINTERPZ[LGC,LTM] = (obs_pc[LTM]-LOW)/DIFF
      for LTM in range(obs_l):
         if obs_pc[LTM] > mod_pc[0]:
            HINTERPZ[0,LTM] = 1.0
            HINTERPZ[1:mod_l,LTM] = 0.0
   if mode == 'PC':
      HINTERPZ= np.ndarray(shape=(mod_l+1, obs_l+1), dtype=np.float64, order='C')
      HINTERPZ[:,:]= 0.0
      for LGC in range(mod_l):
         for LTM in range(obs_l):
            LOW = mod_pc[LGC+1]
            HI = mod_pc[LGC]
            if HI <= obs_pc[LTM]:
               if LOW >= obs_pc[LTM+1]:
                  HINTERPZ[LGC,LTM] = 1.0
               elif HI >= obs_pc[LTM+1]:
                  DIFF = HI - LOW
                  HINTERPZ[LGC,LTM] = (HI-obs_pc[LTM+1])/DIFF
                  HINTERPZ[LGC,LTM+1] = (obs_pc[LTM+1]-LOW)/DIFF
            elif HI > obs_ps:
               HINTERPZ[LGC,0] = 1.0
   return(HINTERPZ)

# calculate pressure levels based on sigma levels
def geos5_pressure(gc_lat,gc_lon,ps):
    tot_level = 48
    p = np.zeros((tot_level,len(gc_lat),len(gc_lon)))
    for i in range(len(gc_lat)):
        for j in range(len(gc_lon)):
           for k in range(tot_level):
               p[k,i,j] = AP[k] + ps[i,j]*BP[k]
    return(p)

# get lat lon index of the model based on observed lat and lon.
def get_gc_ij(obs_lat,obs_lon):
   lat_gc_projection = obs_lat+88
   if lat_gc_projection <=0:
      mod_lat_ind = 0
   else:
      mod_lat_ind = int(lat_gc_projection/4)+1
   lon_gc_projection = obs_lon+182.5
   if lon_gc_projection>=360:
      mod_lon_ind = 0
   else:
      mod_lon_ind = int(lon_gc_projection/5)
   return(mod_lat_ind,mod_lon_ind)
# get lat lon index of the nested grid model on NA
def get_gc_ij_NA(obs_lat,obs_lon):
   lat_gc_projection = obs_lat-9.75
   if lat_gc_projection <=0:
      mod_lat_ind = 0
   else:
      mod_lat_ind = int(lat_gc_projection/0.5)
   lon_gc_projection = obs_lon+140+0.3333
   if lon_gc_projection<=0:
      mod_lon_ind = 0
   else:
      mod_lon_ind = int(lon_gc_projection/0.6667)
   return(mod_lat_ind,mod_lon_ind)
# get vert lvl index of the model based on observed altitude in km      
def get_gc_lvl(obs_alt):
   GEOS5_ALT = np.array([0.058,0.189,0.320,0.454,0.589,0.726,0.864,1.004,1.146,1.29,
                      1.436,1.584,1.759,1.988,2.249,2.517,2.792,3.074,3.439,3.896,
                      4.375,4.879,5.413,5.98,6.585,7.237,7.943,8.846,9.936,11.021,
                      12.086,13.145,14.17,15.198,16.222,17.243,18.727,20.836])
   result = int(0)
   for i in range(len(GEOS5_ALT)-1):
      if (obs_alt > GEOS5_ALT[i]) and (obs_alt < GEOS5_ALT[i+1]):
         if -GEOS5_ALT[i] + obs_alt <= GEOS5_ALT[i+1] - obs_alt:
            result = i
         else:
            result = int(i+1)
      elif obs_alt > GEOS5_ALT[-1]:
         result = int(len(GEOS5_ALT))
      elif obs_alt < GEOS5_ALT[0]:
         result = int(0)
   return(result)
# get vert lvl index for the atom observations
def get_atom_obs(obs1,obs2):
   obs = 0
   if obs1 > 0 and obs1 < 900:
      if obs2 <= 0:
         obs = obs1
      elif obs2 > 900:
         obs = obs1
      else:
         obs = 0.5*(obs1 + obs2)
   elif obs2 > 0 and obs2 < 900:
      obs = obs2
   return(obs)
