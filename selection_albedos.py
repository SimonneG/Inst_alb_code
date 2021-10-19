#!/usr/bin/env python3

import os
import sys
import argparse
import warnings
import string
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
import numpy as np
from mpl_toolkits.basemap import Basemap
import scipy as sp
from scipy.interpolate import griddata
#from pyhdf.SD
import matplotlib as mpl 
import matplotlib.pyplot as plt
import sys
import netCDF4
from netCDF4 import Dataset
from pyhdf.SD import SD,SDC
from functools import reduce
#import julian
#import datetime

parser = argparse.ArgumentParser(description='Plot POLDER/PARASOL data')
parser.add_argument('-f','--filename', help='Polder filename', required=True)
parser.add_argument('-c','--ceresdate',help='Date for CERES file',required=True)
args = vars(parser.parse_args()) 
filename=args['filename']
ceres_date=args['ceresdate']
#day=float(filename[29:31])
#month=float(filename[32:34])
month=int(filename[74:76])
day=int(filename[77:79])
if month==4:
    day=day+31
elif month==5:
    day=day+61
elif month==7:
    day=day+30
elif month==8:
    day=day+61
elif month==10:
    day=day+30
elif month==11:
    day=day+61
elif month==1:
    day=day+31
elif month==2:
    day=day+62
#filename='/LARGE/PARASOL/RB2-HDF.v1.01/2005/2005_07_01/POLDER3_L2B-RGB-014008M_2005-07-01T12-37-59_V1-01.h5'
folder_cer = '/LARGE/CERES/CERES_SSF1deg-1H_Aqua-MODIS_Ed4.1/'
file_cer = 'CERES_SSF1deg-1H_Aqua-MODIS_Ed4.1_Subset_'+ceres_date+'.nc'
dataset=Dataset(folder_cer+file_cer)
#print(dataset)
#sys.exit()

# ------- POLDER --------------
prod = 'albedo_SW'
file_in=h5py.File(filename)

# Getting POLDER data
# albedo : offset of 0., scale factor of 0.0001
alb_link = file_in["Data_Fields"][prod]
alb_sf = alb_link.attrs.get("scale_factor")
alb_off = alb_link.attrs.get("add_offset")
# SZA : offset of 0.2, scale factor of 0.004 
SZA_link=file_in["Data_Fields"]['mus']
SZA_sf = SZA_link.attrs.get("scale_factor")
SZA_off = SZA_link.attrs.get("add_offset")
# latitude and longitude
lat_pol_link = file_in["Geolocation_Fields"]["Latitude"]
lat_pol = lat_pol_link[420:660,:]
lon_pol_link = file_in["Geolocation_Fields"]["Longitude"]
lon_pol = lon_pol_link[420:660,:]
# time : minutes and hour
hour_link=file_in["Data_Fields"]['utc_hour']
minutes_link=file_in["Data_Fields"]['utc_minute']

alb_pol = alb_link[420:660,:]*alb_sf+alb_off
SZA_pol = SZA_link[420:660,:]*SZA_sf+SZA_off

land_sea_flag = file_in["Geolocation_Fields"]["land_sea_flag"][420:660,:]
# uncomment below to have sea or land
#alb_pol[np.where(land_sea_flag >= 2)] = np.nan
# 1 : ocean 2 : ocean+coasts, 3 : coasts5 : land + coasts, 6 : land 

# minutes and hours to call right ceres files
hour_TEMP = hour_link[420:660,:] # 
min_TEMP = minutes_link[420:660,:] #

min_h,max_h=hour_link.attrs.get("actual_range")
min_m, max_m=minutes_link.attrs.get("actual_range")
h_cer2 = int(np.round(max_h,0))
h_cer = int(np.round(min_h,0))
#print(min_h,min_m,max_h,max_m)

for ind_h in range(h_cer,h_cer2+1):
    alb_cer=dataset.variables['toa_alb_clr_1h'][(day-1)*24+ind_h,:]
    alb_cer_A = alb_cer[:,0:180]
    alb_cer_B = alb_cer[:,180:360]
    alb_cer = np.append(alb_cer_B,alb_cer_A,1)

    time=dataset.variables['time'][ind_h]

    sza_cer=dataset.variables['solar_zen_angle_1h'][(day-1)*24+ind_h,:]
    sza_cer_A = sza_cer[:,0:180]
    sza_cer_B = sza_cer[:,180:360]
    sza_cer = np.append(sza_cer_B,sza_cer_A,1)

    tfv_cer=dataset.variables['fov_time_1h'][(day-1)*24+ind_h,:]
    tfv_cer_A = tfv_cer[:,0:180]
    tfv_cer_B = tfv_cer[:,180:360]
    tfv_cer = np.append(tfv_cer_B,tfv_cer_A,1)

    aux_oc_cer=dataset.variables['aux_ocean_1h'][(day-1)*24+ind_h,:]
    aux_A = aux_oc_cer[:,0:180]
    aux_B = aux_oc_cer[:,180:360]
    aux_oc_cer = np.append(aux_B,aux_A,1)
    #res_oc = np.where(aux_oc_cer>=50) # 100 is ocean !
    #res_oc_lat = res_oc[0]
    #res_oc_lon = res_oc[1]
    #alb_cer[res_oc_lat,res_oc_lon] = -999.0
    #alb_cer[np.where(aux_oc_cer<=10)] = -999.0

    lon_cer=dataset.variables['lon'][:]
    lon_cer_A = lon_cer[0:180]
    lon_cer_B = lon_cer[180:360]
    lon_cer = np.append(lon_cer_B,lon_cer_A)

    lat_cer=dataset.variables['lat'][:]

    result_pol = np.where(alb_pol<0.5)
    result_alb = np.where(np.logical_and(np.greater_equal(alb_cer,0),np.less_equal(alb_cer,1))) # lat : de -20 à +20, lon de 0 à 180 puis -180 à 0

# polder  0----------180----------360 *6  ->   -180----------0------------180
# ceres   0----------180----------360     ->     0---------180/-180----------0 
# il faut une correspondance :
    # l_lat_cer = 0  --> k_lat_pol = 6*(40-0) = 6*(40-k_lat_cer)
    # k_lon_cer = 0  --> k_lon_pol = 6*(0+180) = 6*(k_lon_cer + 180)

    result_alb_lat=np.unique(result_alb[0])
    result_alb_lon=np.unique(result_alb[1])

    range_lat = np.ndarray.tolist(result_alb_lat)
    range_lon = np.ndarray.tolist(result_alb_lon) 

    alb_pol[alb_pol>0.99]=np.nan
    alb_pol = np.ma.masked_array(alb_pol,np.isnan(alb_pol))


    vec_cer = []
    vec_pol = []
    for k in range_lat: #latitude
        for l in range_lon: #longitude
            if alb_cer[k,l]>0:
                if alb_pol[6*(39-k),l*6]>0:
    #            vec_cer = np.append(vec_cer,alb_cer[k,l])
     #           vec_pol = np.append(vec_pol,np.mean(alb_pol[6*(39-k):6*(39-k)+6,l*6:l*6+6]))
           #    vec_pol = np.append(vec_pol,alb_pol[6*(39-k),l*6])
#              print(k,l)

                    print(alb_cer[k,l],np.nanmean(alb_pol[6*(39-k):6*(39-k)+6,l*6:l*6+6]),lat_cer[k],lon_cer[l],np.nanmean(lat_pol[6*(39-k):6*(39-k)+6,l*6:l*6+6]),np.nanmean(lon_pol[6*(39-k):6*(39-k)+6,l*6:l*6+6]),tfv_cer[k,l],np.nanmean(hour_TEMP[6*(39-k):6*(39-k)+6,l*6:l*6+6]),np.nanmean(min_TEMP[6*(39-k):6*(39-k)+6,l*6:l*6+6]),sza_cer[k,l],np.nanmean(SZA_pol[6*(39-k):6*(39-k)+6,l*6:l*6+6]))

sys.exit()


#test : il faudrait faire np.mean(alb_pol[6*(39-k):6*(39-k)+6),l*6:l*6+6])

#for i in range (0,6):
             #   for j in range (0,6):

                  #vec_alb_cer = np.append(vec_alb_cer,alb_cer[k,l-180])
                   #vec_alb_pol = np.append(vec_alb_pol,alb_pol[k*6+i,l*6+j])
sys.exit()
# ______________END OF PROGRAM : beginning of plot _________________


lonMin=-180
lonMax=180
latMin=-90
latMax=90
res=1/6

lonGrid = np.arange(lonMin,lonMax,res)
latGrid = np.arange(latMin,latMax,res)

lonGrid,latGrid = np.meshgrid(lonGrid,latGrid)

# ---------------- grid values of product --------------------------

#prodGrid = griddata((lon.ravel(),lat.ravel()),prod_Data.ravel(),(lonGrid,latGrid), method = 'linear')

#------------------- plot with Basemap -------------------------------------
#sys.exit()

m = Basemap(projection='sinu',lon_0=0,resolution='l')
m.drawcoastlines(color='k',linewidth=0.3)
#plot_Prod = np.ma.masked_invalid(prod_Data)
#plot_Prod = np.ma.masked_invalid(prodGrid)
cmap = plt.cm.jet
cmap.set_over('w')

fig = plt.gcf()

# Set figure width and height (inches)
fig.set_size_inches(16.1,10.)

x, y= m(lonGrid,latGrid)
im=m.pcolormesh(x,y,prod_Data,cmap=cmap)


plt.rc('text', usetex=True)
plt.rc('font',family='serif')

cbar=m.colorbar(im, 'bottom')   #---------------------COLORBAR-----
cbar.set_label(''+units+'',fontsize='10') #-----------------------


plt.clim(0,0.8)
#plt.clim(np.min(range_Data), np.max(range_Data))

plt.savefig('ffilename', dpi=100)
#plt.show()

