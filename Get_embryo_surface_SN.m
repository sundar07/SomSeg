function Get_embryo_surface_SN()
% This code detects nuclei in fused images and performs a sphere fit of the point cloud
% Created by Sundar Naganathan - June 2018

% Go to the folder, where all images are stored
% Pre-requisites: All images must be in ics/ids format. Use FIJI for converting tiff to ics format. 
% The reason for using this format is because it allows files of larger size to be opened my matlab.

Folder_path = cd; % get the Folder_path of the current folder

% The following parameters should be changed according to imaging conditions
% For an image with a pixel size of 0.46 um, the following parameters work well.
Nuclei_radius = 9; 
Nuclei_out_radius = 19; % The number of px around nuclei centroid - to be considered for signal to noise calculation.
Nuclei_out_radius_half = floor(Nuclei_out_radius/2);

d0 = 10; % Minimum distance between nuclei

% Make a new folder called 'Sphere_fit' if it does not exist already
if exist(strcat(Folder_path,'/Sphere_fit'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Sphere_fit'))
end

% Read all ics/ids files, get their names and sort them in alphabetical order
Files = dir(strcat(Folder_path,'/*ch_0.ids'));
File_names = {Files.name}';
for i = 1:numel(Files)
    Num_chars(i) = numel(File_names{i});
end
[~,idx] = sort(Num_chars);
    
File_names = File_names(idx);

for frame = 1:numel(File_names)
    
    % Download bfopen from Mathworks to open ics/ids files
    data = bfopen(strcat(Folder_path,'/',File_names{frame}));
    
    % Dummy initialization
    Nuclei = [];
            
    for z = 142%1:5:size(data{1,1},1) % Perform nuclei detection in every 5th slice. Usually this is more than enough to get a fit sphere.
        
        I = cell2mat(data{1,1}(z));
        
        % Go through this slice only if the fluorescence intensities are high. 1200 is an arbitrary number. 
        % Open an image and check the gray value of a nucleus and decide this number accordingly.
        % For my images, if the max intensity is lower than this threshold, usually there are no nuclei in the image.
        if max(I(:)) > 1200 
            
            % Convert to 8-bit and double class
            I = im2double(255 * mat2gray(I));
            
            % Filter the image using a laplacian of Gaussian
            h = fspecial('log',15,5); 
            Ifilt = -imfilter(I,h);
            Ifilt(Ifilt < 0) = 0;
            
            Ifilt_dil = imdilate(Ifilt,strel('disk',3)); % Dilate the filtered image
            I_points = (Ifilt == Ifilt_dil) .* Ifilt; % Identify points of local maxima

            I_points(I_points > 0) = 1;

            % Remove points from the border
            I_points(size(I,1)-Nuclei_out_radius_half+1:end,:) = 0; % bottom
            I_points(1:Nuclei_out_radius_half,:) = 0; % top
            I_points(:,size(I,2)-Nuclei_out_radius_half+1:end) = 0; % right
            I_points(:,1:Nuclei_out_radius_half) = 0; % left
            
            [row,col] = find(I_points == 1);
            
            h2 = fspecial('average',Nuclei_radius);
            I_points_avg_h2 = imfilter(I_points,h2); % Apply the average filter with nuclei radius on local maxima image, which will be used as a mask
            I_points_avg_h2(I_points_avg_h2 > 0) = 1; 
            I_mult_h2 = I .* I_points_avg_h2; % Multiply the mask with original image
            I_mult_h2_avg = imfilter(I_mult_h2,h2); % Apply the average filter on masked original image

            I_std = stdfilt(I_mult_h2,true(Nuclei_radius)); % Determine local std of image based on nuclei radius 

            % Do the same steps as above for nuclei out radius
            h3 = fspecial('average',Nuclei_out_radius);
            I_points_avg_h3 = imfilter(I_points,h3);
            I_points_avg_h3(I_points_avg_h3 > 0) = 1;

            I_mult_h3 = I_points_avg_h3 .* I;
            % Multiply the complement to remove those pixels, which were considered for nuclei radius and that leaves us with only those pixels surrounding a nucleus
            I_mult_h3_update = I_mult_h3 .* imcomplement(I_points_avg_h2); 
            I_mult_h3_update(I_mult_h3_update == 0) = nan;

            points_int = [];
            for points = 1:numel(row)
                
                % The following numbers are determined based on trial and error. A std of 5 and an average pixel intensity of greater than 60 signified nuclei regions
                if I_std(row(points),col(points)) > 5 && I_mult_h2_avg(row(points),col(points)) > 60
                    % Get the list of pixels surrounding a potential nucleus
                    I_out = I_mult_h3_update(row(points)-Nuclei_out_radius_half:row(points)+Nuclei_out_radius_half,col(points)-Nuclei_out_radius_half:col(points)+Nuclei_out_radius_half);
                    I_out_mean = nanmean(I_out(:));
                    
                    % Signal to noise ratio - SNR
                    SNR = (I_mult_h2_avg(row(points),col(points)) - I_out_mean)/I_std(row(points),col(points));

                    if ~isinf(SNR) && SNR > 2 % SNR > 2 worked well for our nuclei images. Try out different values to see if you get better results
                        points_int = [points_int;points];
                    end
                end
            end

            row = row(points_int);
            col = col(points_int);
            
            
            
            %%%%%%%%% Remove too close nuclei (likely to be noise) using delaunay triangulation %%%%%%%%%
            
            Points_to_delete = [];
            if numel(row) >= 20
                dt = delaunayTriangulation(col, row);
                
%                 IC = incenter(dt);
%                 triplot(dt)
%                 hold on
%                 plot(IC(:,1),IC(:,2),'*r')
                
                edge = edges(dt);  % vertex ids in pairs
                pwise = reshape(dt.Points(edge', :), 2, size(edge,1), 2);

                difference = pwise(1,:,:) - pwise(2,:,:);
                edge_lengths = sqrt(difference(1,:,1).^2 + difference(1,:,2).^2); % Get all edge lengths

                idx = find(edge_lengths < d0); % Find the idx of those edge lengths that is less than minimum distance between nuclei
                
                % Loop through each idx and determine the points that need to be deleted
                for ed = 1:numel(idx)
                    P1 = dt.Points(edge(idx(ed),1),:);
                    P2 = dt.Points(edge(idx(ed),2),:);

                    [~,idx_pt] = min([I(P1(2),P1(1)) I(P2(2),P2(1))]);
                    Points_to_delete = [Points_to_delete;edge(idx(ed),idx_pt)];
                end

                Points_to_delete = unique(Points_to_delete);
                dt.Points(Points_to_delete,:) = [];
                
                if ~isempty(dt.Points)
                    I1 = false(size(I)); % Initialize a dummy image of the same size as I
                    I1(sub2ind(size(I1),dt.Points(:,2),dt.Points(:,1))) = 1; % Detected nuclei positions are represented as true values and the rest as false

                    [Nuc_row,Nuc_col] = find(I1==1);
                    Nuclei = [Nuclei; Nuc_col Nuc_row repmat(z,size(Nuc_row,1),1)];

                    % Visualize and save detected points superimposed on original image
                    figure(1);clf;imagesc(I,[50 255])
                    hold on
                    plot(Nuc_col,Nuc_row,'Marker','o','Markersize',1,'LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','none','LineWidth',0.1)
                    saveas(gcf,strcat(Folder_path,'/Sphere_fit/Detected_nuclei_fr_',num2str(frame),'_z',num2str(z)),'tif')
                end
            end % if numel(Nuclei_row) >= 20
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end % if max(I(:)) > 1200
    end % for z
    
    Pt_Cloud.Nuclei = Nuclei;
    
    Nuclei = Pt_Cloud.Nuclei;
    if ~isempty(Nuclei)        
        
        % Fit a sphere to the point cloud
        % ellipsoid_fit was downloaded from Mathworks developed by Yury Petrov (Sep 2015)
        [center,radii,evecs,v,chi2] = ellipsoid_fit(Nuclei,'xyz'); 
        Pt_Cloud.Center = center;
        Pt_Cloud.Radii = radii;
        Pt_Cloud.evecs = evecs;
        Pt_Cloud.v = v;
        Pt_Cloud.chi2 = chi2;
        
        %%%%%%%%% Plot the point cloud and the fit sphere %%%%%%%%%
        
        x = Nuclei(:,1);
        y = Nuclei(:,2);
        z = Nuclei(:,3);
        figure(2);clf;
        plot3(x,y,z,'.r' );
        hold on;
        %draw fit
        mind = min([x y z])-200; % for having a closed sphere, decrease the bounds
        maxd = max([x y z])+200; % for having a closed sphere, increase the bounds
        nsteps = 50;
        step = (maxd - mind) / nsteps;
        [x,y,z] = meshgrid(linspace(mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps));
        Ellipsoid = v(1)*x.*x + v(2)*y.*y + v(3)*z.*z + ...
        2*v(4)*x.*y + 2*v(5)*x.*z + 2*v(6)*y.*z + ...
        2*v(7)*x + 2*v(8)*y + 2*v(9)*z;
        p = patch( isosurface(x,y,z,Ellipsoid,-v(10) ) );
        hold off;
        set(p,'FaceColor','g','EdgeColor','none' );
        view(34.4,14.4);
        axis vis3d equal;
        camlight;
        lighting phong;
        
        % Save image of sphere fit
        saveas(gcf,strcat(Folder_path,'/Sphere_fit/Frame_',num2str(frame)),'tif')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        cprintf('comment',['Frame ' num2str(frame) ' completed\n'])
    else
        cprintf('comment',['Frame ' num2str(frame) ' - Embryo image not available or unable to segment\n'])
    end
    
    % Save the point cloud
    save(strcat(Folder_path,'/Sphere_fit/PointCloud_frame_',num2str(frame),'.mat'),'Pt_Cloud')
    clear Pt_Cloud Nuclei
end % for frame

%% Get sphere parameters from each of the sphere fits

Sphere_params.Frame_rate = 5; % min

Matfiles = dir(strcat(Folder_path,'/Sphere_fit/PointCloud_frame_*.mat'));

for i = 1:numel(Matfiles)
    load(strcat(Folder_path,'/Sphere_fit/PointCloud_frame_',num2str(i),'.mat'))
    
    Sphere_params.Radius(i) = Pt_Cloud.Radii(1);
    Sphere_params.Center(i,:) = Pt_Cloud.Center;
    
end
save(strcat(Folder_path,'/Sphere_fit/Sphere_params.mat'),'Sphere_params')

Sphere_params.Time = 0:Sphere_params.Frame_rate:(numel(Matfiles)-1)*Sphere_params.Frame_rate;
figure;clf;plot(Sphere_params.Time,Sphere_params.Radius,'Marker','o','LineStyle','--','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0.75 0.75])
hold on
plot(Sphere_params.Time,repmat(nanmean(Sphere_params.Radius),1,length(Sphere_params.Time)),'--r')
set(gcf,'Color','w')
saveas(gcf,strcat(Folder_path,'/Sphere_fit/Fit_sphere_radii'),'tif')
