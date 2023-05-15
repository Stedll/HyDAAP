import os
os.chdir(os.path.dirname(__file__))
import cv2
import time
import numpy as np
import osgeo.gdal as gdal
import scipy.spatial.transform as transform

import dem_functions as DF
import image_functions as IF
import spherical_functions as SF


def gradient_descent(w_init, vert_off, dem_path, data_path, bandwidth, max_iterations, threshold, learning_rate, momentum):
    w = np.asarray(w_init)
    w_history = w

    old_corr = function_eval(w, vert_off, dem_path, bandwidth, data_path)
    f_history = old_corr
    delta_w = np.zeros(w.shape)

    i = 0
    diff = 1.0e9
    delta = 10

    while  i<max_iterations and diff>threshold:
        print("iteration:", i)
        gradient = np.asarray(function_gradient(w, vert_off, dem_path, bandwidth, data_path, old_corr))
        delta_w = learning_rate * gradient * delta + momentum * delta_w
        print(delta_w)
        w = w + delta_w
        old_corr = function_eval(w, vert_off, dem_path, bandwidth, data_path)
        # store the history of w and f
        w_history = np.vstack((w_history, w))
        f_history = np.vstack((f_history, old_corr))
        
        # update iteration number and diff between successive values
        # of objective function
        i += 1
        diff = np.absolute(f_history[-1]-f_history[-2])
    
    return w_history,f_history

def function_gradient(position, vert_off, dem_path, bandwidth, data_path, corr_old):
    delta = 10 # UTM -> meters
    
    from multiprocessing.pool import ThreadPool
    pool = ThreadPool(processes=4)
    t1 = pool.apply_async(function_eval, (position + [delta, 0], vert_off, dem_path, bandwidth, data_path, '1'))
    t2 = pool.apply_async(function_eval, (position + [-delta, 0], vert_off, dem_path, bandwidth, data_path, '2'))
    t3 = pool.apply_async(function_eval, (position + [0, delta], vert_off, dem_path, bandwidth, data_path, '3'))
    t4 = pool.apply_async(function_eval, (position + [0, -delta], vert_off, dem_path, bandwidth, data_path, '4'))

    corr1 = t1.get()
    corr2 = t2.get()
    corr3 = t3.get()
    corr4 = t4.get()

    return np.asarray([(corr1-corr2)/(2*delta), (corr3-corr4)/(2*delta)])


def function_eval(position, vert_off, dem_path, bandwidth, data_path, modifier=''):
    easting = position[0]
    northing = position[1]
    
    print("evaluating correlation at", easting, northing)
    try:
        dem_skyline_360, dem_parameters = DF.extract_skyline_from_dem(dem_path, easting, northing, vert_off, 1200, 30000)
    except ValueError:
        print("Error, exif coordinates not inside specified raster, change dem or specify manually the location")
        quit()

    signal_path = SF.get_coeffs(dem_skyline_360, os.path.join(data_path, "tmp", "signal"+modifier+".dat"), bandwidth)
    corr = os.popen(" ".join(["../SOFT/bin/arss_code_single ", signal_path, os.path.join(data_path, "tmp", "pattern.dat"), str(bandwidth), str(bandwidth), str(bandwidth)])).read().split(',')[1]
    print("correlation", corr)
    os.remove(signal_path)

    return float(corr)

def main():
    print("code not yet implemented for standalone functionalities")

    start = time.time()
    data_path = os.path.join("..", "data")
    img_path = os.path.join(data_path, "img", "hyp", "labwindow_06_el5_BSQ.hdr")
    acq_path = os.path.join(data_path, "img", "hyp", "labwindow_06_el5.azg")
    dem_path = [os.path.join(data_path, "dem", "w50565_s10/w50565_s10.tif"), os.path.join(data_path, "dem", "w51065_s10/w51065_s10.tif")]

    name_of_run = "test"

    vFOV = np.deg2rad(20)
    
    lat = None
    lon = None
    #northing, easting, _, _ = utm.from_latlon(lat, lon)
    easting = 666320.498
    northing = 5103487.682

    vert_off = 2

    visualization = False
    dem = True
    hyperspectral = True
    optimization = True
    projection = False
    correlation = False

    if dem:
        if len(dem_path)>1:
            os.makedirs(os.path.join(data_path, "dem", "merged"), exist_ok=True)
            g = gdal.Warp(os.path.join(data_path, "dem", "merged", "demfile.tif"), dem_path, format="GTiff") # if you want
            g = None # Close file and flush to disk
            dem_path = os.path.join(data_path, "dem", "merged", "demfile.tif")

    if hyperspectral:
        print("extracting skyline from hyperspectral")
        spectral_skyline_360, img_hyp, size, limit, ang_off, FOVs = IF.spectral_read3(img_path, acq_path, vFOV, 0, 1200)
        print("done")

    if hyperspectral and dem and optimization:
        bandwidth = 300
        print("calculating spherical coefficients")
        pattern_path = SF.get_coeffs(spectral_skyline_360, os.path.join(data_path, "tmp", "pattern.dat"), bandwidth)
        print("done")

        max_iterations = 30
        threshold = 0
        learning_rate = 1
        momentum = 0.1

        print("starting position optimization")
        w_history, f_history = gradient_descent([easting, northing], vert_off, dem_path, data_path, bandwidth, max_iterations, threshold, learning_rate, momentum)
        print(w_history)
        print(f_history)
        print("done")


if __name__ == '__main__':
    main()
