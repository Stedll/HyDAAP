import numpy as np
from osgeo import gdal
import scipy
import pandas as pd
import cv2
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

def readraster(filename):
    raster = gdal.Open(filename)
    return raster

def pol2cart(coords):
    y = coords[:,0] * np.cos(coords[:,1])
    x = -coords[:,0] * np.sin(coords[:,1])
    return(np.vstack((x,y)).T)

def interpolation(zz, cor_y, cor_x):
    cor_y = cor_y - 1
    cor_x = cor_x - 1
    return zz[np.floor(cor_y).astype(int), np.floor(cor_x).astype(int)]*abs((cor_y-np.floor(cor_y))*(cor_x-np.floor(cor_x))) + zz[np.floor(cor_y).astype(int), np.ceil(cor_x).astype(int)]*abs((cor_y-np.floor(cor_y))*(cor_x-np.ceil(cor_x))) \
            + zz[np.ceil(cor_y).astype(int), np.floor(cor_x).astype(int)]*abs((cor_y-np.ceil(cor_y))*(cor_x-np.floor(cor_x))) + zz[np.ceil(cor_y).astype(int), np.ceil(cor_x).astype(int)]*abs((cor_y-np.ceil(cor_y))*(cor_x-np.ceil(cor_x))) 
    

def createvertexarray(raster, radius, lat, lon):
    transform = raster.GetGeoTransform()

    width = raster.RasterXSize
    height = raster.RasterYSize

    bounds = [transform[0], transform[3], transform[0] + (width * transform[1]), transform[3] + (height * transform[5])]

    zz = raster.ReadAsArray()

    rho = np.arange(500, radius, 10)
    phi = np.arange(0, 2*np.pi, 2*np.pi/(360*6))

    rr, pp = np.meshgrid(rho, phi)
    polar_coords = np.vstack((rr,pp)).reshape([2, -1]).transpose()

    cart_coords = pol2cart(polar_coords)
    cart_coords = cart_coords + [lat, lon]
    
    valid = (cart_coords[:,0] < bounds[2]) & (cart_coords[:,0] > bounds[0]) & (cart_coords[:,1] < bounds[1]) & (cart_coords[:,1] > bounds[3])
    cart_coords = cart_coords[valid]
    polar_coords = polar_coords[valid]

    cc_x, cc_y = np.split(cart_coords, 2, axis=1)
    interpolated = interpolation(zz, (cc_y - transform[3])/transform[5], (cc_x - transform[0])/transform[1])
    vertices_cil = np.hstack((polar_coords, interpolated))

    if (lat < bounds[2]) & (lat > bounds[0]) & (lon < bounds[1]) & (lon > bounds[3]):
        dem_heigth = interpolation(zz, (lon - transform[3])/transform[5], (lat - transform[0])/transform[1])
    else:
        dem_heigth = None

    return vertices_cil, bounds, dem_heigth

def cil2spher(rtz, z_off):
    ptsnew = rtz.copy()
    ptsnew[:,2] = np.arctan2(rtz[:,2] - z_off, rtz[:,0])
    return ptsnew

def crop_circular(vertices, center_x, center_y, offset_z, radius):
    vertices_off = vertices - np.array([center_x, center_y, -offset_z])
    center_index = np.argmin(np.linalg.norm(vertices_off[:,:-1], axis=1))
    vertices_off = vertices_off - np.array([0, 0, vertices_off[center_index, 2]])
    masked = vertices_off[np.linalg.norm(vertices_off[:,:-1], axis=1) <= radius]
    indexes = np.argwhere(np.linalg.norm(vertices_off[:,:-1], axis=1) <= radius)
    return masked, indexes

def extract_skyline_from_dem(path, usr_lat, usr_lon, z_off, out_height = 2000, radius=10000):
    corners = None

    raster = readraster(path)
    vertices, corners, z_off_dem = createvertexarray(raster, radius, usr_lat, usr_lon)

    if usr_lat > corners[2] or usr_lat < corners[0] or usr_lon > corners[1] or usr_lon < corners[3]:
        raise ValueError("Error, exif coordinates not inside specified raster, change dem or specify manually the location")

    sphere = cil2spher(vertices, z_off_dem - z_off)
    numbins = 360*6
    skyline=np.empty(numbins)

    df = pd.DataFrame(sphere)
    df = df.sort_values([1, 2], ascending=[True, False]).drop_duplicates([1]).reset_index(drop=True)
    skyline = df[2].to_numpy()
    bin_idx = (df[1].to_numpy()*(numbins)/(2*np.pi)).astype(int)

    skyline = np.flip(skyline)
    skyline = scipy.signal.medfilt(np.tile(skyline, 3), kernel_size=5)[len(skyline):len(skyline)*2]
    skyline = np.rint(((skyline + (np.pi/2)) * (numbins/(2 * np.pi)))).astype(int)
    image = np.zeros((int(numbins/2), numbins))

    image[-skyline, bin_idx] = 255

    image = cv2.resize(image, (0,0), fx=out_height/np.shape(image)[0], fy=out_height/np.shape(image)[0])
    _, image = cv2.threshold(image,1,255,cv2.THRESH_BINARY)
    return image.astype(np.uint8), [raster.RasterXSize, raster.RasterYSize, raster.GetGeoTransform()[0], raster.GetGeoTransform()[3]-raster.GetGeoTransform()[1]*raster.RasterYSize, raster.GetGeoTransform()[1], raster.GetRasterBand(1).GetNoDataValue()]

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()