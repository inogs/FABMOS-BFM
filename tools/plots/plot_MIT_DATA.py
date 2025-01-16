import sys
import numpy  as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as NC
import shapely.vectorized
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Script to manage restart in a MITgcm+BFM simulation
    It sets automatically start and end times, and restarts
    at the end of the simulation.
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--indicator', '-id',
                                type = str,
                                required = True,
                                default="richnessNORM",
                                help = ''' Choice are: richness, richnessNORM,  shannon, shannonNORM, PO4, TCHL'''

                                )

    parser.add_argument(   '--input_file', '-i',
                                type = str,
                                required = True,
                                help = ''' input file name'''

                                )
    parser.add_argument(   '--temporal_frame', '-tf',
                                type = int,
                                required = False,
                                default=0,
                                help = ''' temporal frame '''
                                )

    parser.add_argument(   '--depth_frame', '-df',
                                type = int,
                                required = False,
                                default=0,
                                help = ''' depth frame '''
                                )

    return parser.parse_args()

def richness(indata,threshold):
    res=(indata > threshold).sum()
    if res>threshold:
       return res
    else:
       return np.nan

def richnessNORM(indata,threshold):
    indata = indata[indata>threshold]
    TOT = indata.sum()
    if TOT < threshold:
       return np.nan
    indata_norm = indata/TOT
    res=(indata_norm > threshold).sum()
    if res>threshold:
       return res
    else:
       return np.nan

def shannon(indata,threshold):
    SUM = indata.sum()
    if SUM < threshold:
       return np.nan

# compute TOT, excluding extinct organisms
    TOT = 0
    for dat in indata:
       if dat >threshold:
           TOT+= dat

# compute shannon index
    S = 0
    for dat in indata: 
       if dat >threshold:
           S+= - dat/TOT*np.log(dat/TOT)

    return S

def shannonNORM(indata,threshold):
    SUM = indata.sum()
    if SUM < threshold:
        return np.nan

# compute TOT, excluding extinct organisms
    TOT = 0
    R   = 0
    for dat in indata:
       if dat >threshold:
           TOT+= dat
           R  += 1.
# compute shannon index
    S = 0
    for dat in indata: 
       if dat >threshold:
           S+= - dat/TOT*np.log(dat/TOT)

    return S/np.log(R)

args = argument()

data =  pd.read_csv("AMT_29_lonlat.csv", sep=",",skiprows = 0, engine='python')

lat=data['lat'].values
lon=-data['lon'].values

#Example of MIT-DARWIN
#netcdf biomass_depth_integrated_2003 {
#dimensions:
#	longitude = 144 ;
#	latitude = 90 ;
#	PlanktonNum = 51 ;
#variables:
#	double longitude(longitude) ;
#       double latitude(latitude) ;
#       double PhyNum(PlanktonNum) ;
#       double CarbonBiomass(latitude, longitude, PlanktonNum) ;
#       double DIN(latitude, longitude) ;
#       double PO4(latitude, longitude) ;
#       double FET(latitude, longitude) ;
#       double SIL(latitude, longitude) ;
year='clim'#args.input_file.split("_")[3][0:4]
t=args.temporal_frame
k=args.depth_frame

#with NC.Dataset(args.input_file,"r") as ncMM:
ncMM = NC.Dataset(args.input_file,"r")
jpi     = ncMM.dimensions['x'].size # Longitude
jpj     = ncMM.dimensions['y'].size # Latitude
lonM    = ncMM.variables['lon'][:,:]
latM    = ncMM.variables['lat'][:,:]
lonMD   = ncMM.variables['lon'][0,:]
latMD   = ncMM.variables['lat'][:,0]
PO4     = ncMM.variables['N1_p'][t,k,:,:]
TCHL    = np.zeros(np.shape(PO4))
CarbonBiomass = np.zeros(np.shape(PO4))
numPFT = 0
for idv,var in enumerate(ncMM.variables.keys()):
   if 'Chl' in var:
        TCHL+=ncMM.variables[var][t,k,:,:]
   if 'P' in var:
       if '_c' in var:
           CarbonBiomass+=ncMM.variables[var][t,k,:,:]
   if 'P' in var:
       if '_c' in var:
        numPFT+=1 
   if 'Z' in var:
       if '_c' in var:
        numPFT+=1
        
TCHL[TCHL<0]=np.nan

PFT_Biomasses = np.zeros((numPFT,jpj,jpi))
PFT_names = []
id_pft = 0
for idv,var in enumerate(ncMM.variables.keys()):
   if 'P' in var:
        if '_c' in var:
            PFT_Biomasses[id_pft,:,:]=ncMM.variables[var][t,k,:,:]
            PFT_names.append(var)
            id_pft+=1
   if 'Z' in var:
        if '_c' in var:
            PFT_Biomasses[id_pft,:,:]=ncMM.variables[var][t,k,:,:]
            PFT_names.append(var)
            id_pft+=1
  

PHYTO_R=np.zeros((jpj,jpi)) 
#richness, richnessNORM,  shannon
indicator=args.indicator
if indicator == 'richness':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=richness(PFT_Biomasses[:,jj,ji],0.001)
elif indicator == 'richnessNORM':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=richnessNORM(PFT_Biomasses[:,jj,ji],0.001)
elif indicator == 'shannon':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=shannon(PFT_Biomasses[:,jj,ji],0.001)
elif indicator == 'shannonNORM':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=shannonNORM(PFT_Biomasses[:,jj,ji],0.001)
elif indicator == 'PO4':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=PO4[jj,ji]
elif indicator == 'TCHL':
    for jj in range(jpj):
        for ji in range(jpi):
            PHYTO_R[jj,ji]=TCHL[jj,ji]
else:
    print("ERROR: indicator name not found")
    sys.exit()


###################
## PLOTTING
###################

plt.figure(figsize=(12,6))
# Atlantic Ocean
# create the projections
orthoATL = ccrs.Orthographic(central_longitude=-30, central_latitude=-5)
geo = ccrs.Geodetic()

# create the geoaxes for an orthographic projection
ax = plt.subplot(131,projection=orthoATL)

# transform lat/lons points to othographic points
points = orthoATL.transform_points(geo, lon, lat)

# plot native orthographic points
ax.plot(points[:, 0], points[:, 1], 'ro')

xx, yy = np.meshgrid(lonMD, latMD)
#xx=lon
#yy=lat
#mask = shapely.vectorized.contains(gdf.dissolve().geometry.item(), xx, yy)

ec=ax.pcolormesh(xx, yy,PHYTO_R,transform=ccrs.PlateCarree())

# plot north pole for reference (with a projection transform)
ax.plot([0], [90], 'b^', transform=geo)

# add coastlines for reference
ax.coastlines(resolution='50m')
ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')
ax.set_global()
cbar = plt.colorbar(ec)
cbar.set_label(indicator)

# AMT
# create the projections
# create the geoaxes for an orthographic projection
ax = plt.subplot(132)

itp = RegularGridInterpolator((latMD, lonMD), PHYTO_R)

# amt points
amt_crd=np.array(list(zip(lat,360+lon)))

f_interped = itp(amt_crd)
ax.plot(f_interped, lat)


#ax.set_ylim(50.,-50.)

# Pacific Ocean
# create the projections
orthoPAC = ccrs.Orthographic(central_longitude=150, central_latitude=-5)
geo = ccrs.Geodetic()

# create the geoaxes for an orthographic projection
ax = plt.subplot(133,projection=orthoPAC)

xx, yy = np.meshgrid(lonMD, latMD)
#xx=lon
#yy=lat
#mask = shapely.vectorized.contains(gdf.dissolve().geometry.item(), xx, yy)

ec=ax.pcolormesh(xx, yy,PHYTO_R,transform=ccrs.PlateCarree())

# add coastlines for reference
ax.coastlines(resolution='50m')
ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')
ax.set_global()
cbar = plt.colorbar(ec)
cbar.set_label(indicator)

mytitle=indicator + ' ' + str(year)
plt.suptitle(mytitle,fontsize=14)
filename=indicator + '_t_'+str(t)+ '_d_' + str(k)+ '_i_' + indicator +'.png'
plt.savefig(filename)
plt.show()
