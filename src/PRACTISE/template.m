vs=1;
cs=1;
is=1;
fin_demW='..\dem\DEM_southtyrol.asc';
fin_folder='..\17112011_configHintereis\';
fin_im='slr_ufs\'; 
if is==1                                                                   %KEEP
    fin_imformat='.png';
end                                          %KEEP
cam(:,1)=[6.142488e+05, 5.163441e+06]; % altitude will be calculated!
if vs==1                                                                   %KEEP
    buffer_radius=100; % in  [m]    
end                                                                        %KEEP
cam(:,2)=[6.184486e+05, 5.157440e+06]; % altitude will be calculated!
cam_off=[6, 40]; 
cam_rol=0; 
cam_hFOV=0.018;
cam_vFOV=0.0148; 