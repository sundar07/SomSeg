function Map_projection_SN()
% This code peforms map projection (equidistant cylindrical or Mercator) of fused spherical data
% Pre-requisites - Image files in ics/ids format and sphere parameters (radius and center) that best fits the sample
% Created by Sundar Naganathan - July 2018

% Go to the folder, where all images are stored in ics/ids format
% You should have already obtained a point cloud representing your sample and performed a sphere fit

Map_parameters.mfile = 'Map_projection_SN'; % Record the name of the mfile

Folder_path = cd; % get the Folder_path of the current folder

% Make a new folder called 'Maps_eqdcylin' under 'Sphere_fit' if it does not exist already
if exist(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin'))
end

load(strcat(Folder_path,'/Sphere_fit/Sphere_params.mat')) % Load sphere parameters mat file

% Centre and radius of sphere obtained from 'Get_embryo_surface_SN' code
Map_parameters.Centre_sphere = round(Sphere_params.Center(1,:));
Map_parameters.Radius_sphere = round(Sphere_params.Radius(1));

Map_parameters.Recenter_map = -120; % Depending on region of interest, recenter the map

% By how much the map has to be scaled
Map_parameters.Scale_factor = 15; % For images downsampled 2 times

% Latitude longitude limits. These values have to be changed according to your specific region of interest in the sample
% The regions of interest in our case does not extend beyond 35 deg in latitude from the equator.
Map_parameters.Lat_start = 35;
Map_parameters.Lat_end = -35;

% In our case, the imaging was not performed the full 360 deg. 160 deg projection was enough.
Map_parameters.Long_start = 0;
Map_parameters.Long_end = 160;

% Projections will be performed for spheres of different radii. Each sphere will be 'radius_step' units apart from each other
Map_parameters.Radius_step = 2;

% Usually I run a dummy map projection by putting arbitrary values in radius_iterate surrounding the sphere radius to know the extent
% to which I want to project and then enter the desired values
Map_parameters.Radius_iterate = 676:2:814;

% Read ics/ids files from both channels, get their names and sort them in alphabetical order
Ch1 = dir(strcat(Folder_path,'/*ch_1.ids'));
Ch1_names = {Ch1.name}';
for i = 1:numel(Ch1)
    Num_chars(i) = numel(Ch1_names{i});
end
[~,idx] = sort(Num_chars);
Ch1_names = Ch1_names(idx);

clear Num_chars
% Do the same as above for channel 0
Ch0 = dir(strcat(Folder_path,'/*ch_0.ids'));
Ch0_names = {Ch0.name}';
for i = 1:numel(Ch0)
    Num_chars(i) = numel(Ch0_names{i});
end
[~,idx] = sort(Num_chars);
Ch0_names = Ch0_names(idx);

%% Map projection

% x_map,y_map are the x,y-coordinates in the projected map respectively
Start_x = Map_parameters.Long_start*Map_parameters.Scale_factor; % Multiply start/end of latitudes and longitudes by scale factor to get the size
End_x = (Map_parameters.Long_end*Map_parameters.Scale_factor)-1;
Start_y = Map_parameters.Lat_start*Map_parameters.Scale_factor; 
End_y = (Map_parameters.Lat_end*Map_parameters.Scale_factor)-1;

x_map = (Start_x:1:End_x);
y_map = (Start_y:-1:End_y);

% Get lambda (longitude) and phi (latitude) that corresponds to each position in the projected map. The formula will vary depending on the desired projection. 
%% Equidistant cylindrical projection
% Inverse formulas for this projection were obtained here: https://mathworld.wolfram.com/CylindricalEquidistantProjection.html

x_map1 = repmat(x_map,length(y_map),1);
y_map1 = repmat(y_map',1,length(x_map));

Lambda = Map_parameters.Recenter_map+x_map1./Map_parameters.Scale_factor;
Phi = y_map1./Map_parameters.Scale_factor;

% Convert to radian
Lambda_rad = deg2rad(Lambda);
Phi_rad = deg2rad(Phi);

%% Mercator projection
% Inverse formulas were obtained here: https://mathworld.wolfram.com/MercatorProjection.html

% x_map1 = repmat(x_map,length(y_map),1);
% y_map1 = repmat(y_map',1,length(x_map));
% 
% y_map1 = deg2rad(y_map1);
% Lambda_rad = deg2rad(Map_parameters.Recenter_map+x_map1./Map_parameters.Scale_factor);
% Phi_rad = (atan(sinh(y_map1./Map_parameters.Scale_factor)));

%% Loop through each time point

for frame = 1:numel(Ch1_names)
    
    % I_map will host the projected images. Clear previous instances of I_map and initialize with NaN.
    clear I_map_ch0 I_map_ch1
    I_map_ch0 = nan(size(Lambda_rad,1),size(Lambda_rad,2),numel(Map_parameters.Radius_iterate));
    I_map_ch1 = nan(size(Lambda_rad,1),size(Lambda_rad,2),numel(Map_parameters.Radius_iterate));
    
    % Download bfopen from Mathworks to open ics/ids files
    data_ch0 = bfopen(strcat(Folder_path,'/',Ch0_names{frame}));
    data_ch1 = bfopen(strcat(Folder_path,'/',Ch1_names{frame}));

    data_dummy = cell2mat(data_ch1{1,1}(1));
    % Initialize I_ch0 and I_ch1, which will contain images to be map projected
    I_ch0 = nan(size(data_dummy,1),size(data_dummy,2),size(data_ch0{1,1},1));
    I_ch1 = nan(size(data_dummy,1),size(data_dummy,2),size(data_ch1{1,1},1));
    clear data_dummy
    % Store all images in mat format and double class
    for slice = 1:size(data_ch1{1,1},1)
        I_ch0(:,:,slice) = im2double(cell2mat(data_ch0{1,1}(slice)));
        I_ch1(:,:,slice) = im2double(cell2mat(data_ch1{1,1}(slice)));
    end
    clear data_ch0 data_ch1

    % Loop through each radius iteration and obtain projections
    for Np = 1:numel(Map_parameters.Radius_iterate)

        % Get x,y,z coordinates that correspond to a certain lambda and phi
        % These are standard spherical to cartesian coordinate conversions
        x = round(Map_parameters.Radius_iterate(Np) .* sin(Lambda_rad) .* cos(Phi_rad) + Map_parameters.Centre_sphere(1));
        y = round(Map_parameters.Centre_sphere(2) - Map_parameters.Radius_iterate(Np) .* sin(Phi_rad));
        z = round(Map_parameters.Centre_sphere(3) - Map_parameters.Radius_iterate(Np) .* cos(Lambda_rad) .* cos(Phi_rad));
        
        % Ensure that there are no indices less than zero or greater than the size of the desired image
        xidx = find(x <= 0 | x > size(I_ch1,2));
        x(xidx) = 1;
        y(xidx) = 1;
        z(xidx) = 1;

        yidx = find(y <= 0 | y > size(I_ch1,1));
        y(yidx) = 1;
        x(yidx) = 1;
        z(yidx) = 1;

        zidx = find(z <= 0 | z > size(I_ch1,3));
        z(zidx) = 1;
        x(zidx) = 1;
        y(zidx) = 1;

        % Get pixel values of the x,y,z coordinates and store them in map projection variables
        I_map_ch0(:,:,Np) = I_ch0(sub2ind(size(I_ch0),y,x,z));
        I_map_ch1(:,:,Np) = I_ch1(sub2ind(size(I_ch1),y,x,z));
        
        % In case there were indices less than zero or greater than the size of the desired image, convert them back to zero
        I_map_ch0_dummy = I_map_ch0(:,:,Np);
        I_map_ch0_dummy([xidx;yidx;zidx]) = 0;
        I_map_ch0(:,:,Np) = I_map_ch0_dummy;

        I_map_ch1_dummy = I_map_ch1(:,:,Np);
        I_map_ch1_dummy([xidx;yidx;zidx]) = 0;
        I_map_ch1(:,:,Np) = I_map_ch1_dummy;
        
        % You need the 'Mapping' toolbox to run the following 5 lines. If you do not have this toolbox, you can comment it out
        % Visualize the map in classic map style
        figure,axesm('eqdcylin','frame','on','grid','on','glinestyle','-','gcolor',[0 1 1],'ParallelLabel','on','MeridianLabel','on','FLatLimit',[-45 45],'FLonLimit',[20 180])
        geoshow(Phi,Lambda,I_map_ch1(:,:,26),'DisplayType','texturemap')
        set(gcf,'Position',[10 10 2000 1000])
        set(gcf,'Color','w')
        colormap(hot(100))
        
        % Save maps in tif format
        A = I_map_ch0(:,:,Np);
        imwrite(A,strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/I_map_ch0_fr_',num2str(frame),'_proj_',num2str(Np),'.tif'),'Compression','none')

        B = I_map_ch1(:,:,Np);
        imwrite(B,strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/I_map_ch1_fr_',num2str(frame),'_proj_',num2str(Np),'.tif'),'Compression','none')
    end
    % Save maps in mat format
    save(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/I_map_ch0_fr_',num2str(frame),'.mat'),'I_map_ch0','-v7.3')
    save(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/I_map_ch1_fr_',num2str(frame),'.mat'),'I_map_ch1','-v7.3')

    cprintf('comment',['Map projection for frame ' num2str(frame) ' completed\n'])
    clear I_ch0 I_ch1
end
% Save parameters that were used for performing map projection
save(strcat(Folder_path,'/Sphere_fit/Maps_eqdcylin/Map_parameters.mat'),'Map_parameters')

