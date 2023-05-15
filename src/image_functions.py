def spectral_read(path, az_vel, VFOV, el_ang, az_ang):
    from spectral.io import envi
    import numpy as np
    import glob
    import time
    img = envi.open(path)

    img_width = img.shape[0]
    img_height = img.shape[1]

    hyp_index = np.empty((4,) + np.shape(img)[:-1])
    # hyp index 0 and 1 are the row and column coordinates of each pixel, images are rotated by 90°, each column is a pixel of the hyperspectral camera, the azimuth pan correspond to the rows

    hyp_index[0,:,:], hyp_index[1,:,:] = np.meshgrid(range(img_height), range(img_width))

    # hyp index 2 and 3 correspond to the direction of arrival of each pixel in star tracker angle coordinates

    hyp_index[2,:,:] = (hyp_index[1,:,:] - (img_width/2) + (az_ang/2)) * az_vel * int(img.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6
    hyp_index[3,:,:] = np.arctan(((hyp_index[0,:,:] - (img_height-1)/2)/((img_height-1)/2)) * np.tan(VFOV/2)) + el_ang

    return img, hyp_index, img_width * az_vel * int(img.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6

def spectral_read2(path, acq_path, vFOV, az_off, resolution, bands=[50,150,160]):
    import spectral
    import numpy as np
    import glob
    import cv2
    import csv
    
    timestamp = []
    azimuth = []
    elevation = []

    with open(acq_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            timestamp.append(float(row[0]))
            azimuth.append(float(row[1]))
            elevation.append(float(row[2]))

    az_vel = np.deg2rad(np.mean(np.diff(azimuth[5:-5])/np.diff(timestamp[5:-5])))
    el_ang = np.deg2rad(np.mean(elevation)-2)

    img = spectral.open_image(path)

    img_width = img.shape[0]
    img_height = img.shape[1]

    bands = img.read_bands(bands)
    bands = cv2.rotate(bands, cv2.ROTATE_90_COUNTERCLOCKWISE)

    for i in range(bands.shape[2]):
        percentile = np.percentile(bands[:,:,i], 99)
        bands[bands[:,:,i] >= percentile, i] = percentile
        bands[:,:,i] = cv2.normalize(bands[:,:,i],  None, 0, 255, cv2.NORM_MINMAX)

    bands = np.rint(bands).astype(np.uint8)
    bands[:,:,0]=255-bands[:,:,0]
    bands = cv2.cvtColor(bands, cv2.COLOR_BGR2GRAY)

    import matplotlib.pyplot as plt

    threshold = np.std(bands, axis=0)<10

    indices = np.concatenate(([0], np.nonzero(threshold[1:] != threshold[:-1])[0]+1))
    labels = threshold[indices]
    arrays_len = [len(part) for part in np.split(threshold, indices[1:])]

    temp=[]
    for idx, index in enumerate(indices):
        if arrays_len[idx]>len(threshold)/100:
            if temp:
                if labels[idx] != threshold[temp[-1]]:
                    temp.append(index)
            else:
                temp.append(index)

    indices = np.asarray(temp)

    arrays = np.split(threshold, indices[1:])
    arrays_len = np.asarray([len(part) for part in arrays])
    labels = np.asarray(threshold[indices])
    
    limits = (indices[labels==False][np.argmax(arrays_len[labels==False])], indices[labels==False][np.argmax(arrays_len[labels==False])] + arrays_len[labels==False][np.argmax(arrays_len[labels==False])])

    img_width = limits[1]-limits[0]

    hyp_index = np.empty((4,) + (img_width, img_height))
    # hyp index 0 and 1 are the row and column coordinates of each pixel, images are rotated by 90°, each column is a pixel of the hyperspectral camera, the azimuth pan correspond to the rows

    hyp_index[0,:,:], hyp_index[1,:,:] = np.meshgrid(range(img_height), range(img_width))

    # hyp index 2 and 3 correspond to the direction of arrival of each pixel in star tracker angle coordinates

    hyp_index[2,:,:] = (hyp_index[1,:,:] - (img_width/2) + (az_off/2)) * az_vel * int(img.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6
    hyp_index[3,:,:] = np.arctan(((hyp_index[0,:,:] - (img_height-1)/2)/((img_height-1)/2)) * np.tan(vFOV/2)) - el_ang

    spectral_skyline_360 = spectral_skyline(bands[:,limits[0]:limits[1]], hyp_index, resolution)

    return img, [int(img_height), int(img_width)], limits, [el_ang, az_off], spectral_skyline_360, hyp_index, img_width * az_vel * int(img.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6

def spectral_read3(path, acq_path, vFOV, az_off, resolution, bands=[50,150,160]):
    import spectral
    import numpy as np
    import glob
    import cv2
    import csv

    ########## calculate real azimuthal speed and elevation from azg file #########
        
    timestamp = []
    azimuth = []
    elevation = []

    with open(acq_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
            timestamp.append(float(row[0]))
            azimuth.append(float(row[1]))
            elevation.append(float(row[2]))

    az_vel = np.deg2rad(np.mean(np.diff(azimuth[5:-5])/np.diff(timestamp[5:-5])))
    el_ang = np.deg2rad(np.mean(elevation)-1.5)

    ########## load hyperspectral image and compute limits  #########

    img_hyp = spectral.open_image(path)

    img_width = img_hyp.shape[0]
    img_height = img_hyp.shape[1]

    
    # load and normalize bands in 99th percentile
    bands = img_hyp.read_bands(bands)
    bands = cv2.rotate(bands, cv2.ROTATE_90_COUNTERCLOCKWISE)

    for i in range(bands.shape[2]):
        percentile = np.percentile(bands[:,:,i], 99)
        bands[bands[:,:,i] >= percentile, i] = percentile
        bands[:,:,i] = cv2.normalize(bands[:,:,i],  None, 0, 255, cv2.NORM_MINMAX)

    bands = np.rint(bands).astype(np.uint8)
    
    # invert blue channel to same behaviour as NIR channels
    bands[:,:,0]=255-bands[:,:,0]

    # convert bands to grayscale
    bands = cv2.cvtColor(bands, cv2.COLOR_BGR2GRAY)

    # threshold on standard deviation for limit computation
    threshold = np.std(bands, axis=0)<10

    # stuff to remove boundaries (walls) at the extremes of the images
    indices = np.concatenate(([0], np.nonzero(threshold[1:] != threshold[:-1])[0]+1))
    labels = threshold[indices]
    arrays_len = [len(part) for part in np.split(threshold, indices[1:])]

    temp=[]
    for idx, index in enumerate(indices):
        if arrays_len[idx]>len(threshold)/100:
            if temp:
                if labels[idx] != threshold[temp[-1]]:
                    temp.append(index)
            else:
                temp.append(index)

    indices = np.asarray(temp)

    arrays = np.split(threshold, indices[1:])
    arrays_len = np.asarray([len(part) for part in arrays])
    labels = np.asarray(threshold[indices])

    # horizontal limit computation    
    hor_lim = (indices[labels==False][np.argmax(arrays_len[labels==False])], indices[labels==False][np.argmax(arrays_len[labels==False])] + arrays_len[labels==False][np.argmax(arrays_len[labels==False])])

    img_width = hor_lim[1]-hor_lim[0]

    hyp_index = np.empty((4,) + (img_width, img_height))

    # hyp index 0 and 1 are the row and column coordinates of each pixel, images are rotated by 90°, each column is a pixel of the hyperspectral camera, the azimuth pan correspond to the rows
    hyp_index[0,:,:], hyp_index[1,:,:] = np.meshgrid(range(img_height), range(img_width))

    # hyp index 2 and 3 correspond to the direction of arrival of each pixel in star tracker angle coordinates
    hyp_index[2,:,:] = (hyp_index[1,:,:] - (img_width/2) + (az_off/2)) * az_vel * int(img_hyp.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6
    hyp_index[3,:,:] = np.arctan(((hyp_index[0,:,:] - (img_height-1)/2)/((img_height-1)/2)) * np.tan(vFOV/2)) - el_ang

    spectral_skyline_360 = spectral_skyline(bands[:,hor_lim[0]:hor_lim[1]], hyp_index, resolution)

    return spectral_skyline_360, img_hyp, [int(img_height), int(img_width)], hor_lim, [el_ang, az_off, az_vel], [img_width * az_vel * int(img_hyp.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6, vFOV]

def image_skyline(path, vertical_res, HFOV, resolution):
    import cv2
    import numpy as np
    src = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    blue_channel = src[:,:,0]
    blue_channel = cv2.resize(blue_channel, (0,0), fx=vertical_res/np.shape(blue_channel)[0], fy=vertical_res/np.shape(blue_channel)[0])
    blue_channel = cv2.normalize(blue_channel, None, 0, 255, cv2.NORM_MINMAX)

    temp = cv2.morphologyEx(blue_channel.copy(), cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5,5)))
    temp = cv2.morphologyEx(temp, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (10,10)))
    temp = cv2.medianBlur(temp, 9)
    edges = cv2.Canny(temp,50,100)
    edges = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_RECT, (3,3)))

    _, labels, stats, _ = cv2.connectedComponentsWithStats(edges, None, None, None, 8)
    candidates = sorted(np.argwhere(stats[:,cv2.CC_STAT_WIDTH]>=np.shape(edges)[0]/10).reshape(-1), key = lambda x: stats[x,cv2.CC_STAT_WIDTH], reverse=True)
    candidates = np.delete(candidates, 0)
    skyline = (labels==candidates[0]).astype(np.uint8)*255
    skyline = cv2.morphologyEx(skyline, cv2.MORPH_DILATE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3)))

    image = prosp_to_equirectangular(skyline, HFOV, resolution)
    return image

def spectral_skyline(bands, hyp_index, size):
    import cv2
    import numpy as np

    temp = cv2.GaussianBlur(bands, (11,11), 1)
    temp = cv2.medianBlur(temp, 11)

    _, thresh = cv2.threshold(bands,80,255,cv2.THRESH_BINARY)

    thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_RECT, (9,9)))

    edges = cv2.Canny(thresh,20,100) 
    

    _, labels, stats, _ = cv2.connectedComponentsWithStats(edges, None, None, None, 8)

    candidates = sorted(np.argwhere(stats[:,cv2.CC_STAT_WIDTH]>=np.shape(edges)[0]/10).reshape(-1), key = lambda x: stats[x,cv2.CC_STAT_WIDTH], reverse=True)
    candidates = np.delete(candidates, 0)
    
    points = []
    for candiate in candidates[:5]:
        points = points + [[hyp_index[2, point[1], point[0]], hyp_index[3, point[1], point[0]]] for point in np.argwhere(labels==candiate)]
    points = np.asarray(points)

    position = np.vstack([np.digitize(points[:,1], np.arange(-np.pi/2, np.pi/2, np.pi/size)), np.digitize(points[:,0], np.arange(-np.pi, np.pi, np.pi/size))])
    image_hyp = np.zeros((size, size*2))

    image_hyp[position[0], position[1]]=255
    
    #img = cv2.merge((bands, bands, bands))
    #for candidate in candidates[:5]:
    #    skyline = (labels==candidate).astype(np.uint8)*255
    #    skyline = cv2.morphologyEx(skyline, cv2.MORPH_DILATE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11,11)))
    #    img[skyline!=0, 2] = 255

    #cv2.imshow("test", bands)
    #cv2.waitKey(0)   
    #quit()

    return image_hyp

def prosp_to_equirectangular(img, hFOV, out_height = 2000):
    import numpy as np
    import cv2
    img_height = img.shape[0]
    img_width = img.shape[1]
    out_width = out_height*2

    vFOV = float(img_height) / img_width * hFOV

    w_len = np.tan(np.radians(hFOV/ 2.0))
    h_len = np.tan(np.radians(vFOV/ 2.0))

    x,y = np.meshgrid(np.linspace(-180, 180, out_width),np.linspace(90, -90, out_height))
    
    x_map = np.cos(np.radians(x)) * np.cos(np.radians(y))
    y_map = np.sin(np.radians(x)) * np.cos(np.radians(y))
    z_map = np.sin(np.radians(y))

    xyz = np.stack((x_map,y_map,z_map), axis=2)
    xyz = xyz.reshape([out_height , out_width, 3])
    inverse_mask = np.where(xyz[:,:,0]>0,np.uint8(1),np.uint8(0))
    xyz[:,:] = xyz[:,:]/np.repeat(xyz[:,:,0][:, :, np.newaxis], 3, axis=2)
    
    
    lon_map = np.where((-w_len<xyz[:,:,1])&(xyz[:,:,1]<w_len)&(-h_len<xyz[:,:,2])&(xyz[:,:,2]<h_len),(xyz[:,:,1]+w_len)/2/w_len*img_width,0)
    lat_map = np.where((-w_len<xyz[:,:,1])&(xyz[:,:,1]<w_len)&(-h_len<xyz[:,:,2])&(xyz[:,:,2]<h_len),(-xyz[:,:,2]+h_len)/2/h_len*img_height,0)
    mask = np.where((-w_len<xyz[:,:,1])&(xyz[:,:,1]<w_len)&(-h_len<xyz[:,:,2])&(xyz[:,:,2]<h_len),np.uint8(1),np.uint8(0))

    persp = cv2.remap(img, lon_map.astype(np.float32), lat_map.astype(np.float32), cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)
    
    mask = mask * inverse_mask
    persp = persp * mask
    _, persp = cv2.threshold(persp,1,255,cv2.THRESH_BINARY)
    persp = cv2.morphologyEx(persp, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11,11)))
    persp = cv2.morphologyEx(persp, cv2.MORPH_DILATE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3)))
    
    return persp   

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()