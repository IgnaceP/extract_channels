""" Script to retrieve geospatial data on channel extent """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon as shPol
from osgeo import gdal
import cv2

NDVI = 'Sunderbands_NDVI.tif'

def importRaster(file):
    Raster = gdal.Open(file)
    Band=Raster.GetRasterBand(1)
    Array=Band.ReadAsArray()
    return(Raster, Band, Array)

NDVI_geo, NDVI_bands, NDVI = importRaster(NDVI)

# treshold value to distinguish water based on NDVI
binary_treshold = -0.5
simplefy_parameter = 1
Tresh = np.zeros(np.shape(NDVI))

#######################
### Pre-processing  ###
#######################
# apply filter to cancel noise without affecting edges
blur = cv2.bilateralFilter(NDVI,3,5,5)

# make a binary image
(binary_treshold, binary) = cv2.threshold(src =blur,
                                          thresh = binary_treshold,
                                          maxval = 255,
                                          type = cv2.THRESH_BINARY)


# convert to proper data type
binary = np.uint8(binary)

# contour
contours = cv2.findContours(image = binary, mode = cv2.RETR_LIST, method = cv2.CHAIN_APPROX_SIMPLE)
contours = contours[0]
contours_len = []
print("Found %d objects." % len(contours))

# list how long each contour is
for (i, c) in enumerate(contours):
    contours_len.append(len(c))
    print("\tSize of contour %d: %d" % (i, len(c)))


# create empty lists to store the polygons
patches = []
sh_Pols = []

# minimum length of contour to generate polygon
tresh = 20
# generate polygons
binary_treshold = 0
print('Generating Polygons...')
for C in contours:
    binary_treshold += 1

    if len(C)> tresh:
        # adapt the dimensions
        C = C[:,0,:]

        # create a polygon
        pol = Polygon(C, closed = True, facecolor = 'red', edgecolor = 'red', alpha = 0.05) # plt polygon
        patches.append(pol)

        shpol = shPol(C) # shapely polygon
        shpol = shpol.simplify(simplefy_parameter)
        # single polygon
        if type(shpol.buffer(0)) is shPol:
            if shpol.buffer(0).length > tresh:
                sh_Pols.append(shpol.buffer(0))
        # multipolygon
        else:
            for s in list(shpol.buffer(0)):
                if s.length > tresh:
                    sh_Pols.append(s)

# create a collection of polygons
p = PatchCollection(patches, match_original = True, alpha = 0.5, edgecolor = 'Darkgreen', facecolor = 'Maroon')

# create a figure
fig, ax = plt.subplots(figsize = (15,15))
ax.imshow(NDVI, cmap = 'Greys')
ls = dir(ax)
ax.add_collection(p)
fig.savefig('test.png')

"""
# get the major polygon
sh_Pols_len = []
for s in sh_Pols:
    sh_Pols_len.append(s.length)

C_maj = sh_Pols[np.argmax(sh_Pols_len)]

# get the interiors
interiors = []
binary_treshold = 0


print('Getting the interior boundaries...')
def checkWithin_init(c,C_maj):
    if c.within(C_maj):
        if c.area < C_maj.area:
            return c

def checkWithin(c, interiors_init):
    int_init_counter = 0
    for int_init in interiors_init:
        if c.within(int_init):
            int_init_counter += 1
    if int_init_counter <= 1:
        return c

# multiprocessing to
pool = mp.Pool(mp.cpu_count())
interiors_init = pool.starmap(checkWithin_init, [(c, C_maj) for c in sh_Pols])
interiors_init = [i for i in interiors_init if i]
interiors = pool.starmap(checkWithin, [(c, interiors_init) for c in interiors_init])
interiors = [i for i in interiors if i]
pool.close()

print(f'Got them! There are {len(interiors)} interior polygons')

# translate to latlon
interiors_latlon = []

for i in interiors:
    pol_latlon = fun.shapelyInd2shapelyLatLon(i, LatLon)
    interiors_latlon.append(pol_latlon)

pol_maj_latlon = fun.shapelyInd2shapelyLatLon(C_maj, LatLon)

path = 'Google_Earth_Files/ChannelsMachala.kml'
fun.shapely2KML(pol_maj_latlon,pol_int = interiors_latlon, path = path, name = 'test')
"""
