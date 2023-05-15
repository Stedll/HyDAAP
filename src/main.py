import exiftool
import utm
import os
import cv2
import time
import numpy as np
import osgeo.gdal as gdal
import scipy.spatial.transform as transform

os.chdir(os.path.dirname(__file__))

import dem_functions as DF
import image_functions as IF
import spherical_functions as SF
#import position_optimizer as PO
import json_writer as JW
import dem_reprojector as DR
import geotiff_writer as GW


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
    easting = 666229.654
    northing = 5103761.858

    vert_off = 2

    visualization = True
    dem = True
    hyperspectral = True
    optimization = True
    projection = False
    correlation = True
    mapping=True
    writer=False

    # control flow
    b_dem_skyline = False
    b_hyp_skyline = False
    b_correlation = False

    if dem:
        if len(dem_path)>1:
            os.makedirs(os.path.join(data_path, "dem", "merged"), exist_ok=True)
            g = gdal.Warp(os.path.join(data_path, "dem", "merged", "demfile.tif"), dem_path, format="GTiff") # if you want
            g = None # Close file and flush to disk
            dem_path = os.path.join(data_path, "dem", "merged", "demfile.tif")

        print("extracting skyline from dem")
        try:
            dem_skyline_360, dem_parameters = DF.extract_skyline_from_dem(dem_path, easting, northing, vert_off, 1200, 30000)

        except ValueError:
            print("Error, coordinates not inside specified raster, change dem or specify manually the location")
            quit()
        print("done")

        b_dem_skyline = True
        if visualization:
            cv2.imshow("dem skyline", dem_skyline_360)
            cv2.waitKey(1)

    if hyperspectral:
        print("extracting skyline from hyperspectral")
        spectral_skyline_360, img_hyp, size, limit, ang_off, FOVs = IF.spectral_read3(img_path, acq_path, vFOV, 0, 1200)
        print("done")
        
        b_hyp_skyline = True
        if visualization:
            cv2.imshow("spectral skyline", spectral_skyline_360)
            cv2.waitKey(1)

    if b_hyp_skyline and b_dem_skyline and correlation:
        bandwidth = 400
        print("calculating spherical coefficients")
        signal_path = SF.get_coeffs(dem_skyline_360, os.path.join(data_path, "tmp", "signal.dat"), bandwidth)
        pattern_path = SF.get_coeffs(spectral_skyline_360, os.path.join(data_path, "tmp", "pattern.dat"), bandwidth)
        print("done")

        print("spherical cross-correlation")
        ret = os.popen(" ".join(["../SOFT/bin/arss_code_single ", signal_path, pattern_path, str(bandwidth), str(bandwidth), str(400)])).read().split(',')
        print("correlation:", ret[1])
        print("done")

        b_correlation = True

        euler = np.asarray(ret[2:]).astype(float)
        ypr = transform.Rotation.from_euler('ZYZ', euler).as_euler('ZYX', degrees=True)

        if visualization:
            photo_rotated = SF.rotate_coeffs(spectral_skyline_360, euler, True)[:-1,:-1]

            dem_skyline_360_filterd = np.zeros_like(dem_skyline_360, dtype=np.float64)
            dem_skyline_360_filterd[photo_rotated<=30] = dem_skyline_360[photo_rotated<=30]
            dem_skyline_360_filterd[spectral_skyline_360<=30] = dem_skyline_360_filterd[spectral_skyline_360<=30]
            img = cv2.merge((dem_skyline_360_filterd, dem_skyline_360_filterd + spectral_skyline_360, dem_skyline_360_filterd + photo_rotated))

            cv2.imshow("spherical matching", cv2.resize(img, None, fx=0.5, fy=0.5))
            cv2.waitKey(1)

    if projection and b_correlation:
        print("creating matlab input files in:")
        matlab_out_dir = JW.create_input(dem_path, img_path, os.path.join(data_path, "tmp"), name_of_run, ypr[0]+0.5, ypr[1], -ypr[2], FOVs[0], FOVs[1], size[0], size[1], dem_parameters, easting, northing, str(vert_off))
        print(matlab_out_dir)
        print("done")

        print("matlab run")
        cmd = "matlab -nodisplay -r \"run PRACTISE/PRACTISE_spher.m;exit\""
        os.system(cmd)
        print("done")

    if mapping:
        print("reprojection")
        correspondences = DR.reproject(img_hyp, size, limit, ang_off, "../data/out/test/PRACTISE")
        print("done")  

        if visualization and b_hyp_skyline:
            bands=img_hyp.read_bands([50,150,160])    
            bands = cv2.rotate(bands, cv2.ROTATE_90_COUNTERCLOCKWISE)
            bands = bands[:,limit[0]:limit[1]]

            for i in range(bands.shape[2]):
                percentile = np.percentile(bands[:,:,i], 99)
                bands[bands[:,:,i] >= percentile, i] = percentile
                bands[:,:,i] = cv2.normalize(bands[:,:,i],  None, 0, 255, cv2.NORM_MINMAX)

            bands = np.rint(bands).astype(np.uint8)
            bands[:,:,0]=255-bands[:,:,0]
            for center in list(zip(correspondences[3],correspondences[2])):
                bands = cv2.circle(bands, center, 2, (0,0,255), 2)
                
            cv2.imwrite("projection.png", cv2.resize(bands, None, fx=0.5, fy=1))
            cv2.imshow("test", cv2.resize(bands, None, fx=0.5, fy=1))
            cv2.waitKey(0)

    if writer:
        print("writing GeoTIFF file")
        GW.writer2(correspondences, img_hyp, limit, dem_path, dem_parameters)
        print("done")

    end = time.time()
    print("execution time:", end-start)

if __name__ == '__main__':
    main()