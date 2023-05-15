def reproject(spectral, size, limits, offsets, matlab_out_dir):
    import scipy as sp
    import os
    import numpy as np

    img_height = size[0]
    img_width = size[1]

    dem_to_spher = sp.io.loadmat(os.path.join(matlab_out_dir, "projection.mat"))['demW'] 
    # 0-1 longitude and latitude
    # 2   altitude
    # 3-4 raster row and column
    # 5-6 spherical directions for each visible DEM point

    #print(limits)
    col_index = np.rint(np.asarray((dem_to_spher[5] / (offsets[2] * int(spectral.metadata['description'].split('\n')[0].split(' ')[-1]) / 1e6)) + img_width/2)).astype(int)
    row_index = np.rint((img_height/2)*np.tan(offsets[0]-dem_to_spher[6])/np.tan(np.deg2rad(20)/2) + img_height/2).astype(int)
    #np.deg2rad(offsets[0])

    #print(np.min(row_index))
    #print(np.max(row_index))
    #print(np.min(col_index))
    #print(np.max(col_index))

    hyp_col_index = col_index[np.logical_and(np.logical_and(col_index>=limits[0], col_index<limits[1]), np.logical_and(row_index>=0, row_index<img_height))]
    hyp_row_index = row_index[np.logical_and(np.logical_and(col_index>=limits[0], col_index<limits[1]), np.logical_and(row_index>=0, row_index<img_height))]

    dem_col_index = dem_to_spher[4, np.logical_and(np.logical_and(col_index>=limits[0], col_index<limits[1]), np.logical_and(row_index>=0, row_index<img_height))].astype(int)
    dem_row_index = dem_to_spher[3, np.logical_and(np.logical_and(col_index>=limits[0], col_index<limits[1]), np.logical_and(row_index>=0, row_index<img_height))].astype(int)

    #print(np.min(hyp_col_index))
    #print(np.max(hyp_col_index))
    #print(np.min(hyp_row_index))
    #print(np.max(hyp_row_index))

    indexes = np.stack([dem_row_index, dem_col_index, hyp_row_index, hyp_col_index])

    return indexes

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()