% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cam, cam_rc] = Cam_PRACTISE(cam, header, dem, cam_off)
%   Name:       Cam_PRACTISE
%   Purpose:    Calculate camera position using longitude, latitude and camera
%               offset height
%                   1) altitude in DEM raster and world coordinates
%                   2) in raster (rows, cols)
%
%   Output:     cam (longitude, latitude and altitude of view and target point)
%               cam_rc (cols and rows in DEM raster, altitude)
%   Input:      cam (longitude, latitude and altitude of view and target point, 
%                 --optional without altitude)
%               header (values od DEM ASCII-file header [cf. ESRI ArcGIS]:
%                 ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
%   --optional  dem (DEM raster from ASCII-file [cf. ESRI ArcGIS])
%               cam_off (camera viewpoint and targetpoint offset height)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Camera positions 
if nargin<2 || nargin==3 || nargin>4
    error('Wrong number of input arguments in function Cam_PRACTISE.')
else
    %   viewpoint in rows and columns
    cam_rc=NaN(3,2);
    cam_rc(1,1)=round((cam(1,1)-header(3,1)+(0.5*header(5,1)))/header(5,1)); % longitude -> col
    cam_rc(2,1)=round((-cam(2,1)+header(4,1)+(header(2,1)*header(5,1)) ...
        +(0.5*header(5,1)))/header(5,1)); % latitude -> row
    if nargin==4
        cam(3,1)=dem(cam_rc(2,1),cam_rc(1,1))+cam_off(1); % DEM(latitude,longitude)
    end
    cam_rc(3,1)=cam(3,1)/header(5,1);
    %   target in rows and columns
        cam_rc(1,2)=round((cam(1,2)-header(3,1)+(0.5*header(5,1)))/header(5,1)); % longitude -> col
        cam_rc(2,2)=round((-cam(2,2)+header(4,1)+(header(2,1)*header(5,1)) ...
            +(0.5*header(5,1)))/header(5,1)); % latitude -> row
    if nargin==4
        cam(3,2)=dem(cam_rc(2,2),cam_rc(1,2))+cam_off(2); % DEM(latitude,longitude)
    end
    cam_rc(3,2)=cam(3,2)/header(5,1);
end