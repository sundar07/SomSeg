function Map_projection_pixel_map_SN()
% This function determines the pixel size in the projected map
% Pre-requisite: You must have already performed sphere fit and map projection and therefore you must have the following mat files:
% Sphere_params.mat and Map_parameters.mat
% Created by Sundar Naganathan - July 2018

Folder_path = cd;
% Load the following mat files
load(strcat(Folder_path,'/Sphere_fit/Sphere_params.mat'))
load(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/Map_parameters.mat'))

Map_pixel.Downsample = 2; % By how much images were downsampled when performing fusion of multiview images
Map_pixel.Original_px_size = 0.228; % um per pixel
Map_pixel.Original_px_size_down = Map_pixel.Original_px_size * Map_pixel.Downsample; % two times downsampling
Map_pixel.Scale_factor = Map_parameters.Scale_factor; 
Map_pixel.Recenter_map = Map_parameters.Recenter_map; 

% The different iterations of sphere radii used
Map_pixel.Radius_sphere = Map_parameters.Radius_iterate;

% Latitude longitude limits
Map_pixel.Lat_start = Map_parameters.Lat_start;
Map_pixel.Lat_end = Map_parameters.Lat_end;

Map_pixel.Long_start = Map_parameters.Long_start;
Map_pixel.Long_end = Map_parameters.Long_end;

% x_map,y_map are the x,y-coordinates in the projected map respectively
Start_x = Map_pixel.Long_start*Map_pixel.Scale_factor; % In my case, the imaging was not performed the full 360 deg. ~160 deg projection was enough.
End_x = (Map_pixel.Long_end*Map_pixel.Scale_factor)-1;
Start_y = Map_pixel.Lat_start*Map_pixel.Scale_factor; % The regions of interest in my case does not extend beyond ~45 deg from the equator.
End_y = (Map_pixel.Lat_end*Map_pixel.Scale_factor)-1;

x_map = (Start_x:1:End_x);
y_map = (Start_y:-1:End_y);

%% Equidistant cylindrical projection
% x_map1 = repmat(x_map,length(y_map),1);
y_map1 = repmat(y_map',1,length(x_map));

% Lambda = Map_pixel.Recenter_map+x_map1./Map_pixel.Scale_factor;
Phi = y_map1./Map_pixel.Scale_factor;

% Lambda_rad = deg2rad(Lambda);
Phi_rad = deg2rad(Phi);

% Map_pixel.Circle_radii = nan(size(Phi_rad,1),size(Phi_rad,2),numel(Map_pixel.Radius_sphere));
Map_pixel.Pixel_size = nan(size(Phi_rad,1),size(Phi_rad,2),numel(Map_pixel.Radius_sphere));

for i = 1:numel(Map_pixel.Radius_sphere)
%     Map_pixel.Circle_radii(:,:,i) = Map_pixel.Radius_sphere(i) .* cos(Phi_rad);
    Map_pixel.Pixel_size(:,:,i) = Map_pixel.Circle_radii(:,:,i) .* (2*pi*Map_pixel.Original_px_size_down/(360*Map_pixel.Scale_factor)); 
end

% Get the pixel size at the equator across maps
Pixel_size.Equator = Map_pixel.Pixel_size(round(size(Phi_rad,1)/2),1,:);
Pixel_size.Equator = Pixel_size.Equator(:);

save(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/Map_pixel.mat'),'Map_pixel','-v7.3')
save(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/Pixel_size.mat'),'Pixel_size')


% Plot pixel size for full extent and mark pixel size at 45 deg
% figure,plot(Map_pixel.Pixel_size(901:end,1))
% set(gcf,'Color','w')
% set(gca,'XLim',[0 901])
% hold on
% line([450,450],[0.1 0.6],'Color','red','LineStyle','--')
% line([0,900],[0.4187 0.4187],'Color','red','LineStyle','--')
% saveas(strcat(Folder_path,'/Sphere_fit/Maps_scalefactor10/PixelSize_variation_45deg_mark'),'tif')