import numpy as np
import datetime
import scipy
import fabmos.transport.tmm
import fabmos
import netCDF4

tm_config_dir = "."  # directory with a TM configuration from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/

domain = fabmos.transport.tmm.create_domain(tm_config_dir)
lon_tgt, lat_tgt = domain.lon[1::2, 1::2], domain.lat[1::2, 1::2]

#NN interpolation:

with netCDF4.Dataset("BiogeochemData/dust.orca.nc") as nc, netCDF4.Dataset("BiogeochemData/dust.MITgcm_2.8.nc", "w") as ncout:
    ncout.createDimension("x", domain.nx)
    ncout.createDimension("y", domain.ny)
    ncout.createDimension("time_counter", len(nc.dimensions["time_counter"]))
    ncout.createVariable("lon", lon_tgt.dtype, ("y", "x"))[...] = lon_tgt
    ncout.createVariable("lat", lat_tgt.dtype, ("y", "x"))[...] = lat_tgt
    newnctime = ncout.createVariable("time_counter", float, ("time_counter",))
    newnctime.units = "days since 2000-01-01 00:00:00"
    newnctime.calendar = "360_day"
    newnctime[...] = netCDF4.date2num([datetime.datetime(2000, m+1, 16) for m in range(12)], newnctime.units)

    lon = np.asarray(nc.variables["nav_lon"]).ravel()
    lat = np.asarray(nc.variables["nav_lat"]).ravel()
    lonex = np.concatenate((lon - 360, lon, lon + 360))
    latex = np.tile(lat, 3)
    points = np.stack((lonex, latex), axis=-1)
    tri = scipy.spatial.Delaunay(points)
    for name, ncvar in nc.variables.items():
        if name not in ("nav_lon", "nav_lat", "time_counter"):
            values = np.reshape(ncvar, ncvar.shape[:-2] + (-1,))
            values = np.tile(values.T, [3, 1])
            ip = scipy.interpolate.LinearNDInterpolator(tri, values)
            values_ip = ip(lon_tgt, lat_tgt)
            values_ip = np.moveaxis(values_ip, -1, 0)
            newncvar = ncout.createVariable(name, ncvar.dtype, ("time_counter", "y", "x"))
            newncvar[...] = values_ip
            newncvar.coordinates = "lon lat"
            for k in ncvar.ncattrs():
                if k not in ["_FillValue"]:
                    setattr(newncvar, k, getattr(ncvar, k))
