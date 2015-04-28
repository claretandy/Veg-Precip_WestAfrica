# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 11:35:00 2015

@author: ajh235
"""

import iris
inrp = iris.load_cube("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/4203.pp")
llgrid = iris.load_cube("/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils/km4/qrparm.landfrac_4km_ll.nc")
outll = inrp.regrid(llgrid, iris.analysis.Linear(extrapolation_mode='mask'))