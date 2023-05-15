def interpolation(zz, cor_y, cor_x):
    import numpy as np
    cor_y = cor_y - 1
    cor_x = cor_x - 1
    
    height =  zz[np.floor(cor_y).astype(int), np.floor(cor_x).astype(int)]*abs((cor_y-np.floor(cor_y))*(cor_x-np.floor(cor_x))) + zz[np.floor(cor_y).astype(int), np.ceil(cor_x).astype(int)]*abs((cor_y-np.floor(cor_y))*(cor_x-np.ceil(cor_x))) \
            + zz[np.ceil(cor_y).astype(int), np.floor(cor_x).astype(int)]*abs((cor_y-np.ceil(cor_y))*(cor_x-np.floor(cor_x))) + zz[np.ceil(cor_y).astype(int), np.ceil(cor_x).astype(int)]*abs((cor_y-np.ceil(cor_y))*(cor_x-np.ceil(cor_x))) 
    
    return height
    
def get_height(dictionary, name = 'camera'):
    from osgeo import gdal
    import struct
    ##### open DEM file and get projection information
    DemDs = gdal.Open(dictionary['dem_file_tif'])
    gt=DemDs.GetGeoTransform()
    dem_band=DemDs.GetRasterBand(1)
    ##### identify pixel position of point in raster
    mx = float(dictionary[name+'_DEMCRS_E'])
    my = float(dictionary[name+'_DEMCRS_N'])
    px = (mx - gt[0]) / gt[1] #x pixel position in raster
    py = (my - gt[3]) / gt[5] #y pixel position in raster  
    ##### collect height value (if outside of extent, return: 0)


    if(px>=0 and px<DemDs.RasterXSize and py>=0 and py<DemDs.RasterYSize):
        zz = DemDs.ReadAsArray()
        height = interpolation(zz, py, px)
    else:
        print(name + " point is outside of DEM raster extent.\nThe DEM height is artificially set to 0")
        print("Please check your coordinates to confirm if input was correct.")
        height = 0
    # return height value to main procedure
    return height

def calc_target_from_lookdir(direction, distance = 1000):
    # import required libraries
    import math, sys
    ##### calculate coordinate offset
    if(0<=direction<45) or (315<=direction<=360):
        dy=distance*math.cos(math.radians(direction))
        dx=dy*math.tan(math.radians(direction))
    elif(45<=direction<135):
        dx=distance*math.sin(math.radians(direction))
        dy=-dx*math.tan(math.radians(direction-90))
    elif(135<=direction<225):
        dy=distance*math.cos(math.radians(direction))
        dx=dy*math.tan(math.radians(direction-180))
    elif(225<=direction<315):
        dx=distance*math.sin(math.radians(direction))
        dy=-dx*math.tan(math.radians(direction-270))
    else:
        print('Error: direction angle not within possible range (0, 360)')
        sys.exit()      
    # return offset to main procedure
    return dx, dy

def calc_orientation (dictionary):
    # import required libraries
    import math
    ##### calculate orientation angle from coordinate offset
    dx = float(dictionary['target_DEMCRS_E']) -  float(dictionary['camera_DEMCRS_E'])
    dy = float(dictionary['target_DEMCRS_N']) -  float(dictionary['camera_DEMCRS_N'])
    if (dy>0):
        if(dx>0):
            return math.degrees(math.atan(float(dx)/float(dy)))
        else:
            return math.degrees(math.atan(float(dx)/float(dy))) + 360
    if (dy<0):
            return (math.degrees(math.atan(float(dx)/float(dy))) + 180)
    if (dy==0):
        if(dx>0):
            return 90
        else:
            return 270

def calc_vertical_params(dictionary, offset=True, z=True, pitch = True, height_diff=True):
    # import required libraries
    import math
    
    ##### get DEM height at target point
    dem_height = get_height(dictionary, name='target')
    dictionary['var_add']['target_DEM_height'] = dem_height
    
    ##### calculate remaining parameters
    if (pitch==False) and height_diff:
        # derive height difference from pitch angle
        dictionary['var_add']['height_diff_target_cam'] = float(dictionary['var_add']['distance'])*math.tan(math.radians(float(dictionary['var_add']['pitch_angle_deg'])))
        height_diff=False
    if(z):
        if(height_diff):
            # derive absolute target height and height difference from target offset
            dictionary['var_add']['target_DEMCRS_Z'] = dem_height + float(dictionary['target_offset'])
            dictionary['var_add']['height_diff_target_cam'] = float(dictionary['var_add']['target_DEMCRS_Z']) - float(dictionary['var_add']['camera_DEMCRS_Z'])
        elif(offset):
            # derive absolute target height and target offset from height difference
            dictionary['var_add']['target_DEMCRS_Z'] = float(dictionary['var_add']['camera_DEMCRS_Z']) + float(dictionary['var_add']['height_diff_target_cam'])
            dictionary['target_offset'] = dictionary['var_add']['target_DEMCRS_Z'] - dem_height
    else:
        # derive height difference and target offset from absolute target height
        dictionary['var_add']['height_diff_target_cam'] = float(dictionary['var_add']['target_DEMCRS_Z']) - float(dictionary['var_add']['camera_DEMCRS_Z'])
        dictionary['target_offset'] = float(dictionary['var_add']['target_DEMCRS_Z']) - dem_height
    if pitch:
        # derive pitch angle from height difference
        dictionary['var_add']['pitch_angle_deg'] = math.degrees(math.atan(float(dictionary['var_add']['height_diff_target_cam'])/float(dictionary['var_add']['distance'])))
    
    # return new dictionary to main procedure
    return dictionary

def get_distance_tc(dictionary):
    # import required libraries
    import math
    # calculate distance
    distance = math.sqrt((dictionary['target_DEMCRS_E']-dictionary['camera_DEMCRS_E'])**2+
                         (dictionary['target_DEMCRS_N']-dictionary['camera_DEMCRS_N'])**2)
    # return distance value to main procedure
    return distance

def create_input(dem_file_in, img_file_in, out_dir, name_of_run, yaw, pitch, roll, hFOV, vFOV, img_rows, img_cols, dem_parameters, camera_utm_lon, camera_utm_lat, camera_offset):
    # import required libraries and submodules of georef_webcam
    import collections, json, os, subprocess
    from osgeo import gdal
     
    ##### create dictionary with parameters required for PRACTISE
    var_mand = ['path', 'image_file', 'dem_file', 'gcp_file', 'camera_DEMCRS_E', 'camera_DEMCRS_N', 'camera_offset', 
                'target_DEMCRS_E', 'target_DEMCRS_N', 'target_offset', 'roll_angle_deg', 
                'hFOV', 'vFOV', 'img_rows', 'img_cols', 'buffer_around_camera_m', 'var_add', 'dem_ncols', 'dem_nrows', 'dem_xllcorner', 'dem_yllcorner', 'dem_cellsize', 'dem_nodata_value']
    input_dict = collections.OrderedDict((k, None) for k in var_mand)
    input_dict['var_add'] = collections.OrderedDict()
    
    ##### add existing default values
    parameters = ['hFOV', 'vFOV', 'img_rows', 'img_cols', 'roll_angle_deg', 'buffer_around_camera_m', 'camera_offset', 'target_offset', 'dem_ncols', 'dem_nrows', 'dem_xllcorner', 'dem_yllcorner', 'dem_cellsize', 'dem_nodata_value']
    default = [hFOV, vFOV, img_rows, img_cols, roll, 100.0, 0, 0, dem_parameters[0], dem_parameters[1], dem_parameters[2], dem_parameters[3], dem_parameters[4], dem_parameters[5]]

    for i in range(len(parameters)):
        input_dict[parameters[i]] = default[i]

    ############### define directories to output & input ######################
    # path to store PRACTISE output
    input_dict['path'] = out_dir        
    
    dem_file = dem_file_in
    print("if raster not in EPSG:4326 CRS it will be converted to it")
    if gdal.Info(dem_file, format='json')['coordinateSystem']['wkt'].rsplit("CONVERSION")[0].rsplit("EPSG")[-1].rsplit("]")[0].rsplit(",")[-1] != "4326":
        dem_file_new = os.path.join(os.path.dirname(dem_file),'WGS84_'+os.path.basename(dem_file)[0:-3]+'tif')
        subprocess.call(['gdalwarp', '-of', 'GTiff', '-t_srs', 'EPSG:4326',  dem_file, dem_file_new])
        dem_file = dem_file_new
    
    input_dict['dem_file_tif'] = dem_file

    ##### decrease spatial resolution of DEM file, if you want to fasten the projection calculation
    subsample = False
    if(subsample):
        resolution = float(input('Which spatial resolution should the resampled DEM have? \n'))
        dem_file_int2 = dem_file[0:-4]+'_res' + str(int(resolution)) + '.tif'
        #dem_file_new = dem_file[0:-4]+'_res' + str(int(resolution)) + '.asc'
        subprocess.call(['gdalwarp', '-of', 'GTiff', '-tr', str(resolution), str(resolution), '-r', 'bilinear', dem_file, dem_file_int2])
        #subprocess.call(['gdal_translate', '-of', 'AAIGrid', '-a_nodata', '-3.402823466e+38', dem_file_int2, dem_file_new])
        input_dict['dem_file_tif'] = dem_file_int2
        dem_file = dem_file_new

    input_dict['dem_file'] = dem_file
    input_dict['image_file'] = img_file_in
    
    #same epsg as DEM
    input_dict['var_add']['camera_epsg'] = ''
    input_dict['var_add']['camera_Easting'] = float(camera_utm_lon)
    input_dict['var_add']['camera_Northing'] = float(camera_utm_lat)
    input_dict['var_add']['camera_x'] = input_dict['var_add']['camera_Easting']
    input_dict['var_add']['camera_y'] = input_dict['var_add']['camera_Northing']

    input_dict['camera_DEMCRS_E']=input_dict['var_add']['camera_x']
    input_dict['camera_DEMCRS_N']=input_dict['var_add']['camera_y']
    
    ##### get DEM height at point
    dem_height = get_height(input_dict, 'camera')
    input_dict['var_add']['camera_DEM_height'] = dem_height
    print("The DEM height at the camera location is: "+str(int(dem_height))+ 'm')
    input_dict['camera_offset'] = camera_offset
    input_dict['var_add']['camera_DEMCRS_Z'] = dem_height + float(input_dict['camera_offset'])

    ##### collect orientation angle
    input_dict['var_add']['yaw_angle_deg'] = float((yaw+180)%360)
                                                   
    ##### derive target point location from it
    dx,dy = calc_target_from_lookdir(float(input_dict['var_add']['yaw_angle_deg']))
    input_dict['target_DEMCRS_E'] =  float(input_dict['camera_DEMCRS_E']) + dx
    input_dict['target_DEMCRS_N'] =  float(input_dict['camera_DEMCRS_N']) + dy

    input_dict['var_add']['target_epsg'] = input_dict['var_add']['camera_epsg']

    input_dict['var_add']['target_Easting']=input_dict['target_DEMCRS_E']
    input_dict['var_add']['target_Northing']=input_dict['target_DEMCRS_N']
    input_dict['var_add']['target_x']=input_dict['target_DEMCRS_E']
    input_dict['var_add']['target_y']=input_dict['target_DEMCRS_N']

    # aux: calculate distance between camera and target
    input_dict['var_add']['distance'] = get_distance_tc(input_dict)
    
    input_dict['var_add']['pitch_angle_deg'] = pitch
    input_dict = calc_vertical_params(input_dict, pitch = False)

	############### write resulting dict to json file #########################
    print('Collected parameters are stored in json file:')
    print(os.path.join(out_dir,name_of_run+".json"))
    jsonfile = json.dumps(input_dict)
    f = open(os.path.join(out_dir,name_of_run+".json"),"w")
    f.write(jsonfile)
    f.close()

    print("creating matlab input file:")
    return create_matlab(os.path.join(out_dir,name_of_run+".json"), name_of_run)

def create_matlab(json_path, name_of_run):
    # import required libraries and submodules of georef_webcam
    import os, shutil, json
    with open(json_path) as read_file:
        input_dict = json.load(read_file)
    
    PRACTISE_input=[]
    
    PRACTISE_input.append("fin_demW=" + "'../"+input_dict['dem_file']+"';")
    # define working directory of PRACTISE, where it stores input & output 
    PRACTISE_input.append("fin_folder=" + "'../"+os.path.join(input_dict['path'].replace('tmp', 'out'), str(name_of_run), 'PRACTISE')+"/';")
    # subdirectory where input by PRACTISE is copied to (e.g. image & GCP file)
    PRACTISE_input.append("fin_im=" + "'../"+input_dict['image_file']+"';")

    # insert camera and target point position
    PRACTISE_input.append("cam(:,1)=" + "["+str(input_dict['camera_DEMCRS_E'])+","+str(input_dict['camera_DEMCRS_N'])+"];")
    PRACTISE_input.append("cam(:,2)=" + "["+str(input_dict['target_DEMCRS_E'])+","+str(input_dict['target_DEMCRS_N'])+"];")
    
    # insert height offset values for camera and target
    PRACTISE_input.append("cam_off=" + "["+str(input_dict['camera_offset'])+","+str(input_dict['target_offset'])+"];")
    
    # insert buffer range & roll angle
    PRACTISE_input.append("buffer_radius=" +str(input_dict['buffer_around_camera_m'])+";")
    PRACTISE_input.append("cam_rol=" +str(input_dict['roll_angle_deg'])+";")
    
    # insert horizontal and vertical FOV
    PRACTISE_input.append("cam_hFOV=" +str(float(input_dict['hFOV']))+";")
    PRACTISE_input.append("cam_vFOV=" +str(float(input_dict['vFOV']))+";")

    # insert image size
    PRACTISE_input.append("pix_c=" +str(input_dict['img_cols'])+";")
    PRACTISE_input.append("pix_r=" +str(input_dict['img_rows'])+";")

    with open(os.path.join(input_dict['path'], str(name_of_run) +'.m'), 'w') as matlab_file:
        for line in PRACTISE_input:
            matlab_file.write(line+'\n') 

    with open(os.path.join('PRACTISE','Input_PRACTISE.m'), 'w') as matlab_file:
        matlab_file.write("Input_PRACTISE_m='" + os.path.join('..', input_dict['path'], str(name_of_run)) + "';\n")
        matlab_file.write("run(Input_PRACTISE_m)")
            
    # return file name of input file (without .m file extension as PRACTISE does not need it)
    return os.path.join(input_dict['path'].replace('tmp', 'out'), str(name_of_run), 'PRACTISE')

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()