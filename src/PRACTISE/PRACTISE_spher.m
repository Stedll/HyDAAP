% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
    % incorporating adapted versions of
    %
    % Aronica, Bates and Horritt, 2002 
    % "Assessing the uncertainty in distributed model predictions using
    % observed binary pattern information within GLUE" (Hydrol Process)
    %
    % Corripio, 2004 
    % "Snow surface albedo estimation using terrestrial photography"
    % (Int J Rem Sens)
    %
    % Haerer, Bernhardt, Schulz and Corripio (2012)
    % "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.1.0)"
    % (GMD)
    %
    % Salvatori et al., 2011
    % "Snow cover monitoring with images from digital camera systems"
    % (ItJRS)
    %
    % Tolson and Shoemaker, 2007 
    % "Dynamically dimensioned search algorithm for computationally
    % efficient watershed model calibration" (WRR)
    % 
    % Wang, Robinson and White, 2000 
    % "Generating Viewsheds without Using Sightliines" (PE&RS)
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 22/12/2015)
%       automatically edited by georef_webcam by Sebastian Buchelt (University of Würzburg, 15/06/2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name:       PRACTISE
%   Purpose:    Main file of PRACTISE
%   Comment:    This file starts and controls the complete software tool,
%               be careful with any changes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
tic
%
disp(' ')
disp('Starting PRACTISE - Photo Rectification And ClassificaTIon SoftwarE')
% load input configuration
Input_PRACTISE_m='Input_PRACTISE';
run(Input_PRACTISE_m)
%%%%%%%%%%%%%%%%%%%%%%%%%% check & administration %%%%%%%%%%%%%%%%%%%%%%%%%
% Camera positions (viewpoint and targetpoint) defined?
if ~exist('cam', 'var')
    error(['Camera view- and targetpoint: longitude and latitude ', ...
     'positions are not defined.'])
elseif size(cam, 2)<2
    error('Camera targetpoint: longitude and latitude position is ', ...
     'not defined.')
end
% check if Unix and change "\" to "/" in "fin_..." and "fout_..." in case 
list=who('-regexp', 'fin_');
for i=1:length(list)
    assignin('base',char(list(i)),Unix_Slash_PRACTISE(evalin('base',char(list(i)))));
end
clear list
% Output folder defined and exists?
if ~exist(fin_folder, 'dir')
    mkdir(fin_folder);
end

disp(['The used input file ''', Input_PRACTISE_m, '.m'' will be ', ...
    'saved and a logfile will be produced in ''', fin_folder, ''''])
diary([fin_folder, 'logfile_PRACTISE.txt']);
copyfile([Input_PRACTISE_m, '.m'], fin_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading and preparing input')
% load
%   DEM
%       read LINK file (if existing)
if strcmp(fin_demW(end-8:end), '.dem.LINK')
    fid=fopen(fin_demW);
    if fid==-1
        error('Linking file to DEM could not be found.')
    end
    fin_demW=fgetl(fid);
    fin_demW=[fin_demW, fgetl(fid)];
    fin_demW=Unix_Slash_PRACTISE(fin_demW);
    fclose(fid);
end
%       read DEM file
disp(['Load DEM ''', strrep(fin_demW, 'asc', 'tif'), ''''])

[demWrcz,R] = readgeoraster(strrep(fin_demW, 'asc', 'tif'),'OutputType','double');

headerv_W(2,1) = R.RasterSize(1);
headerv_W(1,1) = R.RasterSize(2);
headerv_W(3,1) = R.XWorldLimits(1);
headerv_W(4,1) = R.YWorldLimits(1);
if R.CellExtentInWorldX == R.CellExtentInWorldY
    headerv_W(5,1) = R.CellExtentInWorldX;
else
    error("Raster does not have square cells, processing cannot continue")
    return
end

headerv_W(6,1) = georasterinfo(strrep(fin_demW, 'asc', 'tif')).MissingDataIndicator;
demWrcz = standardizeMissing(demWrcz,headerv_W(6,1));

clear ans fid

X0=[cam(1,1), cam(2,1), cam_off(1), cam_rol, cam(1,2), cam(2,2), ...
     cam_off(2), cam_hFOV, cam_vFOV];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% viewshed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate camera parameters (raster position and altitude)
[cam, cam_rc]=Cam_PRACTISE(cam, headerv_W, demWrcz, cam_off);
% Specified camera angles viewshed
toc
%   normalize DEM to viewpoint altitude
demWvs=(demWrcz/headerv_W(5,1))-cam_rc(3,1);
%       create radial buffer zone around camera location
%           j=col i=row cam_rc(1,1)=col cam_rc(2,1)=row
[j i]=meshgrid(1:size(demWvs, 2), 1:size(demWvs, 1)); 
buffer_radius_rc=buffer_radius/headerv_W(5,1); 
buffer_cam=sqrt((j-cam_rc(1,1)).^2+(i-cam_rc(2,1)).^2)<=buffer_radius_rc;
demWvs(buffer_cam==1)=min(min(demWvs))-1; % smaller than minimum existing value
clear i j buffer_radius_rc
demWvs(cam_rc(2,1),cam_rc(1,1))=0;
%   Cam properties & viewshed raster initialisation
disp('Start calculation of specified camera angles (+ 0.1rad) viewshed')
%       field of view (use maximum range, if cam_rol~=0 + 0.1rad)
roldummy=deg2rad(cam_rol);%*pi/180; % deg to rad
% FOV_vert=atan((0.5*cam_hei*abs(cos(roldummy))+0.5*cam_wid*abs(sin(roldummy)))/cam_foc)+0.1;
FOV_vert = abs(cam_hFOV * sin(roldummy)) + abs(cam_vFOV * cos(roldummy));
% FOV_hor=atan((0.5*cam_wid*abs(cos(roldummy))+0.5*cam_hei*abs(sin(roldummy)))/cam_foc)+0.1;
FOV_hor = abs(cam_hFOV * cos(roldummy)) + abs(cam_vFOV * sin(roldummy));
clear roldummy
%       calculate azimuth angles (clockwise from Azi(1,1) to Azi(2,1),
%       with N=0, E=pi/2, S=pi, W=3*pi/2)
%           delta DEM in the viewing direction [DEM pixel] -> Cam target-Cam position
N0_rc(1,1)=cam_rc(1,2)-cam_rc(1,1); % delta longitude (=delta cols)
N0_rc(2,1)=cam_rc(2,1)-cam_rc(2,2); % delta latitude (=-delta rows)        
N0_rc(3,1)=cam_rc(3,2)-cam_rc(3,1); % delta altitude
%           calculate angles in viewing direction
if N0_rc(2,1)>0 && N0_rc(1,1)>0 % delta latitude, delta longitude > 0 = NE
    N0_rc_hor=pi/2-atan(N0_rc(2,1)/N0_rc(1,1)); 
elseif N0_rc(2,1)<0 && N0_rc(1,1)>0 % delta latitude < 0 & delta longitude > 0 = SE
    N0_rc_hor=pi/2-atan(N0_rc(2,1)/N0_rc(1,1)); 
elseif N0_rc(2,1)<0 && N0_rc(1,1)<0 % delta latitude, delta longitude < 0 = SW
    N0_rc_hor=3*pi/2-atan(N0_rc(2,1)/N0_rc(1,1));
elseif N0_rc(2,1)>0 && N0_rc(1,1)<0 % delta latitude > 0 & delta longitude < 0 = NW
    N0_rc_hor=3*pi/2-atan(N0_rc(2,1)/N0_rc(1,1)); 
elseif N0_rc(2,1)==0 && N0_rc(1,1)==0
    error(['Camera position and target at the same xy ', ... 
        'coordinate, observing the sky or the soil in nadir ', ...
        'position does not work yet.'])
elseif N0_rc(2,1)==0
    if N0_rc(1,1)>0
        N0_rc_hor=pi/2;
    elseif N0_rc(1,1)<0
        N0_rc_hor=3*pi/2;
    end
elseif N0_rc(1,1)==0
    if N0_rc(2,1)>0
        N0_rc_hor=0;
    elseif N0_rc(2,1)<0
        N0_rc_hor=pi;
    end
end
%           calculate azimuth boundaries (with N=0, E=pi/2, ...)
Azi(1,1)=N0_rc_hor-FOV_hor;
if Azi(1,1)<0
    Azi(1,1)=2*pi+Azi(1,1);
end
Azi(2,1)=N0_rc_hor+FOV_hor;
if Azi(2,1)>2*pi
    Azi(2,1)=Azi(2,1)-2*pi;
end
%           define viewshed sectors with NE=1, SE=2, SW=3, NW=4
for i=1:2
    if Azi(i,1)>=0 && Azi(i,1)<pi/2 %N-E
        Azi(i,2)=1;
    elseif Azi(i,1)>=pi/2 && Azi(i,1)<pi %E-S
        Azi(i,2)=2;
    elseif Azi(i,1)>=pi && Azi(i,1)<3*pi/2 %S-W
        Azi(i,2)=3;
    elseif Azi(i,1)>=3*pi/2 && Azi(i,1)<2*pi %W-N
        Azi(i,2)=4;
    end
end
if Azi(1,1)>Azi(2,1) % Azimuth angle including N (2*pi to 0)
    azi_sec=[Azi(1,2):4, 1:Azi(2,2)];
else
    azi_sec=[Azi(1,2):Azi(2,2)];
end
%    azi_sec_name=['Northeast'; 'Southeast'; 'Southwest'; 'Northwest'];
%    azi_sec_name=azi_sec_name(azi_sec,:);
%    disp('The viewshed for sector(s) ')
%    disp(azi_sec_name)
disp(['will be calculated and no sight barriers within ', ...
    num2str(buffer_radius),  'm distance to the camera location are assumed.'])
%    clear(azi_sec_name)
%       calculate vertical angles (from Vert(1,1) to Vert(2,1); upwards to downwards)
%           upwards=+, horizontal=0, downwards=-
N0_rc_vert=atan(N0_rc(3,1)/norm(N0_rc(1:2,1))); % norm(N0_rc(1:2,1)) is never 0
Vert(1,1)=N0_rc_vert+FOV_vert;
Vert(2,1)=N0_rc_vert-FOV_vert;
if Vert(1,1)>pi/2
    Vert(1,1)=pi/2;
    disp(['Warning: Sky observation angle adjusted to pi/2 as ', ...
        'higher values are not usable. If wanted, use 2 different ', ...
        'viewshed calculations and intersect them!'])
elseif Vert(2,1)<-pi/2
    Vert(2,1)=-pi/2;
    disp(['Warning: Surface observation angle adjusted to -pi/2 ', ...
        'as lower values are not usable. If wanted, use 2 different ', ...
        'viewshed calculations and intersect them!'])
end
clear FOV_hor FOV_vert  
%       viewshed raster initialization
viewW=-1*ones(headerv_W(2,1), headerv_W(1,1));
viewW(isnan(demWvs))=NaN;
%   initialization of temporary rasters and variables
%       R=maximal normalized altitude of visibile pixel
%       index=number-indexed raster (row 1 & col 1 = 1, row 2 & col 1 = 2, ...) 
R=viewW; 
index=reshape([1:headerv_W(2,1)*headerv_W(1,1)], headerv_W(2,1), []);
r=cam_rc(2,1);
c=cam_rc(1,1);
%   viewshed 
%       2 process steps for specified camera angles computation: 
%           1) algorithm assumes that all pixels in the chosen sectors
%              that lie inbetween the vertical angles are visible
%           2) but azimuth angles are saved in case of vertical visibility and
%              taken out later if they do not fulfill the azimuth conditions
%       calculation of the viewpoint and the 8 neighbouring pixels 
viewW(r-1,c+1)=pi/4; %NE
viewW(r+1,c+1)=3/4*pi; %SE
viewW(r+1,c-1)=pi+pi/4; %SW
viewW(r-1,c-1)=pi+3/4*pi; %NW
viewW(r-1,c)=0; %N
viewW(r+1,c)=pi; %S
viewW(r,c-1)=3*pi/2; %W
viewW(r,c+1)=pi/2; %E
for i=-1:1
    for j=-1:1
        if i~=0 || j~=0
            vertdummy=atan(demWvs(r+i,c+j)/sqrt(i^2+j^2));
            if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                viewW(r+i,c+j)=-1;
            end
        else
            if Vert(2,1)>-pi/2
                viewW(r,c)=-1;
            else
                viewW(r,c)=N0_rc_hor; 
            end
        end
    end
end
R(r-1:r+1,c-1:c+1)=demWvs(r-1:r+1,c-1:c+1);
%       calculation of the N, E, S and W directions, seen from the viewpoint
%           using proxy vectors
for dir=1:4
    if dir==1 % N
%               flipud to S direction (=from viewpoint to the south)
        demdummy=flipud(demWvs(1:r-1,c));
        viewdummy=flipud(viewW(1:r-1,c));
        Rdummy=flipud(R(1:r-1,c));
        indexdummy=flipud(index(1:r-1,c));
        azidummy=0;
    elseif dir==2 % E
%               transpose to S direction
        demdummy=demWvs(r,c+1:headerv_W(1,1))'; 
        viewdummy=viewW(r,c+1:headerv_W(1,1))';
        Rdummy=R(r,c+1:headerv_W(1,1))';
        indexdummy=index(r,c+1:headerv_W(1,1))';
        azidummy=pi/2;
    elseif dir==3 % S
        demdummy=demWvs(r+1:headerv_W(2,1),c);
        viewdummy=viewW(r+1:headerv_W(2,1),c);
        Rdummy=R(r+1:headerv_W(2,1),c);
        indexdummy=index(r+1:headerv_W(2,1),c);
        azidummy=pi;
    elseif dir==4 % W
%               transpose and flipud to S direction
        demdummy=flipud(demWvs(r,1:c-1)');
        viewdummy=flipud(viewW(r,1:c-1)');
        Rdummy=flipud(R(r,1:c-1)');
        indexdummy=flipud(index(r,1:c-1)');
        azidummy=3*pi/2;
    end
%           S direction or rotated/mirrored other directions to S direction
    for l=2:size(demdummy, 1)
        Z=Rdummy(l-1,1)*l/(l-1);
        if isfinite(demdummy(l,1)) 
            if demdummy(l,1)>Z
                vertdummy=atan(demdummy(l)/l);
                if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                    viewdummy(l,1)=-1; % invisible
                else
                    viewdummy(l,1)=azidummy; % vertically visible
                end
                Rdummy(l,1)=demdummy(l,1);
            else
                viewdummy(l,1)=-1; % invisible
                Rdummy(l,1)=Z;
            end
        else
            viewdummy(l,1)=NaN; % like invisible, but NaN
            Rdummy(l,1)=Z;
        end
    end
    viewW(indexdummy)=viewdummy;
    R(indexdummy)=Rdummy;
    clear demdummy viewdummy Rdummy indexdummy
end
clear l dir azidummy 
%       calculation of the chosen sectors
%           using proxy arrays
%               W-NW + NW-N = W-N, S-SW + SW-W = S-W,
%               N-NE + NE-E = N-E, E-SE + SE-S = E-S)
for l=1:length(azi_sec)
    dir=azi_sec(l);
    if dir==1 % NE-E & N-NE
%               leftrright array (NE flipped to NW)
        demdummy=fliplr(demWvs(1:r,c:headerv_W(1,1)));
        viewdummy=fliplr(viewW(1:r,c:headerv_W(1,1)));
        Rdummy=fliplr(R(1:r,c:headerv_W(1,1)));
        indexdummy=fliplr(index(1:r,c:headerv_W(1,1)));
        secdummy=-pi/2; % angle clockwise from N=0
    elseif dir==2 % E-SE & SE-S
%               upsidedown + leftrright array (SE flipped to NW)
        demdummy=flipud(fliplr(demWvs(r:headerv_W(2,1),c:headerv_W(1,1))));
        viewdummy=flipud(fliplr(viewW(r:headerv_W(2,1),c:headerv_W(1,1))));
        Rdummy=flipud(fliplr(R(r:headerv_W(2,1),c:headerv_W(1,1))));
        indexdummy=flipud(fliplr(index(r:headerv_W(2,1),c:headerv_W(1,1))));
        secdummy=pi/2; 
    elseif dir==3 % W-SW & S-SW 
%               upsidedown array (SW flipped to NW)
        demdummy=flipud(demWvs(r:headerv_W(2,1),1:c));
        viewdummy=flipud(viewW(r:headerv_W(2,1),1:c));
        Rdummy=flipud(R(r:headerv_W(2,1),1:c));
        indexdummy=flipud(index(r:headerv_W(2,1),1:c));
        secdummy=-3*pi/2;
    elseif dir==4 % W-NW & NW-N
        demdummy=demWvs(1:r,1:c);
        viewdummy=viewW(1:r,1:c);
        Rdummy=R(1:r,1:c);
        indexdummy=index(1:r,1:c);
        secdummy=3*pi/2; 
    end
    rr=size(demdummy, 1);
    cc=size(demdummy, 2);
    m=rr-2;
    n=cc-1;
%           W-NW and NW-N = NW sector or rotated/mirrored other sectors to NW sector
    while ( n>0 ) % cols
        if ( m>0 ) % rows
            if ( rr-m<cc-n ) % "south of NW-diagonal"
                Z=(rr-m)*(Rdummy(m,n+1)-Rdummy(m+1,n+1))+(cc-n)*((rr-m)*(Rdummy(m+1,n+1)-Rdummy(m,n+1))+Rdummy(m,n+1))/(cc-n-1);
            elseif ( rr-m>cc-n ) % "north of NW-diagonal"
                Z=(cc-n)*(Rdummy(m+1,n)-Rdummy(m+1,n+1))+(rr-m)*((cc-n)*(Rdummy(m+1,n+1)-Rdummy(m+1,n))+Rdummy(m+1,n))/(rr-m-1);
            elseif ( rr-m==cc-n ) % "NW-diagonal"
                Z=Rdummy(m+1,n+1)*(rr-m)/(rr-m-1);
            end
            if isfinite(demdummy(m,n))
                if demdummy(m,n)>Z
                    azidummy=abs(atan((rr-m)/(cc-n))+secdummy);
                    vertdummy=atan(demdummy(m,n)/sqrt((rr-m)^2+(cc-n)^2));
                    if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                        viewdummy(m,n)=-1; % invisible
                    else
                        viewdummy(m,n)=azidummy; % vertically visible
                    end
                    Rdummy(m,n)=demdummy(m,n);
                else
                    viewdummy(m,n)=-1; % invisible
                    Rdummy(m,n)=Z;
                end
            else
                viewdummy(m,n)=NaN; % like invisible, but NaN
                Rdummy(m,n)=Z;
            end
            m=m-1;
        else
            m=rr-1;
            n=n-1;
        end
    end
    clear m n
    viewW(indexdummy)=viewdummy;
    R(indexdummy)=Rdummy;
    clear demdummy viewdummy Rdummy indexdummy N0_rc N0_rc_hor ...
        N0_rc_vert
end
clear rr cc secdummy azidummy vertdummy Z dir
clear R
%       visibility decision using azimuth of vertically visible pixels
if Azi(1,1)>Azi(2,1)
    viewW(viewW<Azi(1,1) & viewW>Azi(2,1))=-1;
else
    viewW(viewW<Azi(1,1) | viewW>Azi(2,1))=-1;
end
viewW(buffer_cam==1)=-1;
viewW(viewW>=0)=1;
viewW(viewW<0)=0;
viewAzi=Azi;
viewAziSec=azi_sec;
clear demWvs R index buffer_cam Azi azi_sec N

%%%%%%%%%%%%%%%%%%%%%%% projection %%%%%%%%%%%%%%%%%%%%%%%
toc
disp('Starting projection')
% prepare DEM dependent on the viewshed
%   initialization 
demW=NaN(7, headerv_W(2,1)*headerv_W(1,1)); % if row, col included: 6
k=0;
%   DEM grid point coordinates (middle of pixel)
for i=1:headerv_W(1,1) % ncols
    for j=1:headerv_W(2,1) % nrows
        if viewW(j,i)==1 % visible
            k=k+1;
            demW(1,k)=headerv_W(3,1)-(0.5*headerv_W(5,1))+(i*headerv_W(5,1)); % longitude
            demW(2,k)=headerv_W(4,1)+(headerv_W(2,1)*headerv_W(5,1)) ...
                +(0.5*headerv_W(5,1))-(j*headerv_W(5,1)); % latitude
            demW(3,k)=demWrcz(j,i); %altitude
            demW(4,k)=j; %ascii row
            demW(5,k)=i; %ascii col
        end
    end
end
demW=demW(1:6,isfinite(demW(1,:))); % if row, col included: 1:6
clear i j k 
% prepare projection parameters
%   in case that optimisation off
toc
Xopt=X0;

% project and scale the world coordinates (xyzW in m) to the photo plane 
% (crP in pixels)
%   1) translate world reference system W to reference system T 
%      -> set coordinate system origin = viewpoint 
%   2) rotate reference system T to camera reference system C 
%      -> set coordinate system axis in the directions right (U), up (V) 
%      and depth (N) of image
%   3) project and scale camera reference system C to perspective
%   projection system P
%      -> project to metric cam image using central projection 
%      (at depth of focal length)
%      -> scale metric to pixel based image size
%      -> translate to the origin of the 2D coordinate system of
%      the photograph
%      -> ceil to integer pixel values
%   4) extract RGB pixel values and classify
%      -> keep rows and columns of DEM (Arc/Info ASCII Grid)   
%      -> keep pixel rows and columns of photograph
%      -> extract RGB values for pixel rows and columns
%      -> classify with one of the available methods
%
%   call steps 1) to 3) of DEM projection 
%       

[demP]=Proj_PRACTISE_spher(demW, Xopt, pix_c, pix_r, demWrcz, headerv_W);
% loop for classification of multiple photographs
demW(6:7,:)=demP;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save & write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
disp('start saving and writing files')
%   define output file names
%       variables (MAT-file)
fout_mat=[fin_folder,'projection.mat'];

%   save variables
%       tidy up
clear ans c demWrcz dummy f_name fout_figname1 fout_figname2 ...
 fout_figname4 fout_figname5 l r

save(fout_mat, 'demW'); 
toc

disp('PRACTISE was executed successfully.')                                         