def writer(mapping, img_hyp, limit, dem_path, dem_parameters): # EXTREMELY slow, 1.5 horus to write 200 bands
    import cv2
    import numpy as np
    import osgeo.gdal as gdal

    raster = gdal.Open(dem_path)

    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create("../data/out/output.tif", dem_parameters[0], dem_parameters[1], img_hyp.shape[-1]+1, gdal.GDT_Float32)
    outdata.SetGeoTransform(raster.GetGeoTransform())
    outdata.SetProjection(raster.GetProjection())

    band = outdata.GetRasterBand(1)
    band.WriteArray(raster.ReadAsArray())

    for i in range(img_hyp.shape[-1]):
        print(i)
        band = img_hyp.read_band(i)
        band = cv2.rotate(band, cv2.ROTATE_90_COUNTERCLOCKWISE)
        layer = np.zeros(np.flip(dem_parameters[:2]))
        layer.fill(np.nan)
        layer[tuple(mapping[:2])] = band[limit[0]:limit[1]][tuple(mapping[2:])]
        band = outdata.GetRasterBand(i+2)
        band.WriteArray(layer)
        band.SetNoDataValue(dem_parameters[-1])

def writer2(mapping, img_hyp, limit, dem_path, dem_parameters, name_of_run): # slow but not EXTREMELY slow, 20 minutes to write 200 bands
    import cv2
    import numpy as np
    import rasterio

    with rasterio.Env(GDAL_CACHEMAX=60000000000) as env:
        raster = rasterio.open(dem_path)
        with rasterio.open("../data/out/output_"+name_of_run+".tif", 'w', 'GTiff', dem_parameters[0], dem_parameters[1], img_hyp.shape[-1]+1, raster.crs, raster.transform, rasterio.float32, dem_parameters[-1]) as outfile:
            outfile.write(raster.read(1), 1)
            for i in range(0, img_hyp.shape[-1]):
                print(i+1, "/", img_hyp.shape[-1])
                band = img_hyp.read_band(i)
                band = cv2.rotate(band, cv2.ROTATE_90_COUNTERCLOCKWISE)
                layer = np.zeros(np.flip(dem_parameters[:2]))
                layer.fill(np.nan)
                layer[tuple(mapping[:2])] = band[limit[0]:limit[1]][tuple(mapping[2:])]
                outfile.write(layer, i+1)

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()