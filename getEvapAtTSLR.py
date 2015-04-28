# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 12:22:09 2015

@author: ajh235
"""

from osgeo import gdal
import numpy as np
import iris

evap_file = '/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/3223_gt20060815.nc'
evap = iris.load_cube(evap_file)

tslr_file = '/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/timeSinceLastRain_gdal.nc'
tslr = iris.load_cube(tslr_file)

tsl_data = np.zeros(evap.shape)

tslr = gdal.Open(tslr_file)
for i in range(tslr.RasterCount):
    print i
    band = tslr.GetRasterBand(i+1)
    tsl_data[i,:,:] = band.ReadAsArray().astype(np.float)
    
    
tslr_data = tslr.GetRasterBand(1)


