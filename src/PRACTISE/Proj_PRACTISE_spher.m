% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rpSP] = Proj_PRACTISE_spher(xyzW, X0, pix_c, pix_r, dem, header) 
%   Name:       Proj_PRACTISE
%   Purpose:    Project and scale the world coordinates (xyzW in m) to the photo plane 
%               (crP in pixels)
%               1) translate world reference system W to reference system T 
%                  -> set coordinate system origin = viewpoint 
%               2) rotate reference system T to camera reference system C 
%                  -> set coordinate system axis in the directions right (U), up (V) 
%                  and depth (N) of image
%               3) project and scale camera reference system C to perspective
%               projection system P
%                  -> project to metric cam image using central projection 
%                  (at depth of focal length)
%                  -> scale metric to pixel based image size
%                  -> translate to the origin of the 2D coordinate system of
%                  the photograph
%                  -> ceil to integer pixel values
%   
%   Output:     crP (projected and scaled photo coordinates in pixel; 2D)
%   Input:      xyzW (world coordinates in m; 3D)
%               cam (longitude, latitude and altitude of viewpoint and targetpoint; altitude including offset)
%               cam_rol
%               cam_wid (camera sensor width)
%               cam_hei (camera sensor height)
%               pix_c (number of pixel columns of photograph)
%               pix_r (number of pixel rows of photograph)
%               dem (digital elevation model; grid)
%               header (header of dem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cam(1,1)=X0(1,1);
cam(2,1)=X0(1,2);
cam_off(1)=X0(1,3);
cam_rol=X0(1,4);
cam(1,2)=X0(1,5);
cam(2,2)=X0(1,6);
cam_off(2)=X0(1,7);
cam_hei=X0(1,8);
cam_wid=X0(1,9);

[cam]=Cam_PRACTISE(cam, header, dem, cam_off);
% prepare new projected coordinate system axis (use middle of pixels)   
%   vector N (cam viewing direction)
%       depth in image is positive
N0=cam(:,2)-cam(:,1); %cam target - cam position
N=N0/norm(N0);
%   vector U (horizontal image direction)
%       right is positive; left is negative)
Nxy=[N(1:2,1);0];
if N(3,1) > 0                                                       
    U0=cross(N,(Nxy/norm(Nxy)));
elseif N(3,1) < 0                                                   
    U0=cross((Nxy/norm(Nxy)),N);
elseif N(3,1) == 0 
    V=[0;0;1]; % theoretical V
    U0=cross(V,N);
end
U=U0/norm(U0);
%   vector V (vertical image direction)
%       up is positive; down is negative
V=cross(U,N);
V=V/norm(V);
% prepare scaling from metric to pixel based image size   
scale_c=pix_c/cam_wid;
scale_r=pix_r/cam_hei;

% project
%   initialize for 1)
trans_WtoT=[-cam(1,1); -cam(2,1); -cam(3,1)];
xyzT=NaN(3,size(xyzW,2));
%   1) translate W to T
for i=1:3
    xyzT(i,:)=xyzW(i,:)+trans_WtoT(i,1);
end
clear i
%   initialize for 2)
trans_TtoC=[U(1,1), U(2,1), U(3,1); V(1,1), V(2,1), V(3,1); ...
            N(1,1), N(2,1), N(3,1)];
xyzC=NaN(3,size(xyzW,2));
aerC=NaN(3,size(xyzW,2));
%   2) rotate T to C
for i=1:size(xyzW,2)
    xyzC(1:3,i)=trans_TtoC*xyzT(1:3,i);
    if abs(cam_rol)>0 && abs(cam_rol)<=90 % 0 to 90: clockwise; -90 to 0:anti-clockwise
        if i==1
            roldummy=cam_rol*pi/180;
            trans_Crol=[cos(roldummy), sin(roldummy), 0; ...
                        -sin(roldummy), cos(roldummy), 0; ...
                        0, 0, 1];
            clear roldummy
        end
        xyzC(1:3,i)=trans_Crol*xyzC(1:3,i);
    elseif abs(cam_rol)>90 && i==1
        disp('cam_rol>90deg not defined, use like 0deg')
        cam_rol=0;
    end
end
clear i
%   initialize for 3)
crP=NaN(2,size(xyzW,2));
rpSP=NaN(2,size(xyzW,2));
%   3) ceil([+/-]scale*(project from C to P)+translate); [no mirroring/mirroring]  
%       indexing: no calculation if xyzC(3,:)=0
%aerC(xyzC(3,:)~=0) = cart2sph(xyzC(3,xyzC(3,:)~=0), -xyzC(1,xyzC(3,:)~=0), -xyzC(2,xyzC(3,:)~=0));
rpSP(1,xyzC(3,:)~=0) = atan2(-xyzC(1,xyzC(3,:)~=0), xyzC(3,xyzC(3,:)~=0));
rpSP(2,xyzC(3,:)~=0) = atan2(-xyzC(2,xyzC(3,:)~=0),sqrt(xyzC(1,xyzC(3,:)~=0).^2 + xyzC(3,xyzC(3,:)~=0).^2));
%r = sqrt(xyzC(1,xyzC(3,:)~=0).^2 + xyzC(3,xyzC(3,:)~=0).^2);
%r = sqrt(xyzC(2,xyzC(3,:)~=0).^2 + xyzC(2,xyzC(3,:)~=0).^2 + xyzC(3,xyzC(3,:)~=0).^2);
%rpSP(:,xyzC(3,:)~=0) = aerC(1:2,xyzC(3,:)~=0);

%crP(1,xyzC(3,:)~=0)=ceil(-scale_c*(cam_foc./xyzC(3,xyzC(3,:)~=0).*xyzC(1,xyzC(3,:)~=0))+pix_c/2);
%crP(2,xyzC(3,:)~=0)=ceil(+scale_r*(cam_foc./xyzC(3,xyzC(3,:)~=0).*xyzC(2,xyzC(3,:)~=0))+pix_r/2);