import xarray as xr
import numpy as np
from os import listdir
from os.path import isfile, join

#NOTE: it does not matter which file you use here, they are all on the same grid

data_dir = './CMEMS_data/'
files = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
ds = xr.open_dataset(data_dir + files[0])

R = 6370e3 #radius of Earth

cell_size = float(ds.latitude[1] - ds.latitude[0])

# Get cell edges (coordinates are for cell centers)
lat_adjusted = [float(i) for i in (ds.latitude - (cell_size/2))]
lat_adjusted.append(lat_adjusted[-1]+cell_size)
lat_adjusted = np.array(lat_adjusted)

lon_adjusted = [float(i) for i in (ds.longitude - (cell_size/2))]
lon_adjusted.append(lon_adjusted[-1]+cell_size)
lon_adjusted = np.array(lon_adjusted)

# Lat & lon in radians
phi = lat_adjusted/180.*np.pi
lam = lon_adjusted/180.*np.pi

cell_areas = R*R*np.diff(np.sin(phi))[:, np.newaxis]*np.diff(lam) / 1000000
np.savetxt(data_dir + 'CMEMS_area_file_lat_-10_35_lon_-165_-115.csv',cell_areas,delimiter=",")
