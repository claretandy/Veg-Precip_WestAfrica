import iris
import datetime
from iris.time import PartialDateTime

dt = datetime.datetime(2006, 8, 16)
dt_constraint = iris.Constraint(time=lambda cell: cell >= dt)

precip_cubelist = iris.load('/Volumes/MYBOOK/DataLocal/ModelData/WAFR/4203.pp')
precip = precip_cubelist[0]

olr_cubelist = iris.load('/Volumes/MYBOOK/DataLocal/ModelData/WAFR/2205.pp')
olr = olr_cubelist[0]

with iris.FUTURE.context(cell_datetime_objects=True):
    precip_ss = precip.extract(dt_constraint)
    olr_ss    = olr.extract(dt_constraint)

print("Saving precip ...")
iris.save(precip_ss, '/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/4203.nc', zlib=True)
print("Saving OLR ...")
iris.save(olr_ss, '/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/5min/2205.nc', zlib=True)

