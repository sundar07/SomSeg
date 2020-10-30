function Segment_somites_notochord_SN()
% This code segments notochord and somite boundaries in map projected images
% This code should perform equally well in non-map-projected images as well
% Created by Sundar Naganathan - August 2018

% Go to the folder with images to be segmented
Folder_path = cd; % get the Folder_path of the current folder

Num_projections = 70; % Alternatively number of z-slices in a single time-point
Project = 2; % If you want to perform maximum intensity projection every n-slices and then segment, choose a number here. If you want to segment every slice, then Project = 1
Num_slices = ceil(Num_projections/Project);

% Make a new folder called 'Detected_boundaries' if it does not exist already
if exist(strcat(Folder_path,'/Detected_boundaries'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Detected_boundaries'))
end

% Get info of all mat files that contain images
DirOutput = dir(strcat(Folder_path,'/I_map_ch1*.mat'));
Mat_names = {DirOutput.name};
for i = 1:numel(Mat_names)
    Num_chars(i) = numel(Mat_names{i});
end
[~,idx] = sort(Num_chars);
Mat_names = Mat_names(idx);

% If Seg_slices_info.mat exists, then load it. This means you have already attempted segmentation in this movie and you can get info from this mat file
if exist(strcat(Folder_path,'/Detected_boundaries/Seg_slices_info.mat'),'file')
    load(strcat(Folder_path,'/Detected_boundaries/Seg_slices_info.mat'))
else
    Seg_slices_info.Noto_slices_segmented = false(Num_slices,numel(Mat_names));
    Seg_slices_info.Som_slices_segmented = false(Num_slices,numel(Mat_names));
end

% For the first frame, go through all slices for segmentation. Alternatively, open your image in FIJI and enter the slice numbers where the notochord and somites are first visible
Noto_slice_first = 1;
Noto_slice_end = Num_slices;

Som_slice_first = 1;
Som_slice_end = Num_slices;

% Run through each frame
for frame = 1:numel(Mat_names)
    
    clear rect
    load(strcat(Folder_path,'/',Mat_names{frame}))
    
    %% Perform maximum intensity projection
    
    % Instead of segmenting every slice, you can perform maximum intensity projection of every n slices and segment those instead
    % project every n slices
    I_map = nan(size(I_map_ch1,1),size(I_map_ch1,2),Num_slices);
    I_allbdy = false(size(I_map_ch1,1),size(I_map_ch1,2));
    
    for i = 1:ceil(Num_projections/Project)
        if i == ceil(Num_projections/Project)
            I_map(:,:,i) = max(I_map_ch1(:,:,(i*Project)-(Project-1):end),[],3);
        else
            I_map(:,:,i) = max(I_map_ch1(:,:,(i*Project)-(Project-1):i*Project),[],3);
        end
    end
    
    %% Initialize variables
    
    Somite_seg.Boundary_img = nan(size(I_map,1),size(I_map,2));
    Somite_seg.Noto = nan(size(I_map,1),size(I_map,2));
    Somite_seg.Boundary_identity = nan(50,Num_slices);
    Somite_seg.Noto_up_identity = nan(50,Num_slices);
    Somite_seg.Noto_down_identity = nan(50,Num_slices);
    Somite_seg.Rectangle_analysis = nan(7,Num_slices);
    
    Somite_seg.mfile = 'Segment_somites_notochord_SN';
    
    %% Determine first and last slice to be segmented based on segmentation performed on previous frames
    
    if frame > 1
        
        Noto_slice_first = min(find(Seg_slices_info.Noto_slices_segmented(:,frame-1)))-2;
        if isempty(Noto_slice_first) || Noto_slice_first < 1 
            Noto_slice_first = 1;
        end
        Noto_slice_end = max(find(Seg_slices_info.Noto_slices_segmented(:,frame-1)))+2;
        if isempty(Noto_slice_end) || Noto_slice_end > size(I_map,3)
            Noto_slice_end = size(I_map,3);
        end

        Som_slice_first = min(find(Seg_slices_info.Som_slices_segmented(:,frame-1)))-2;
        if isempty(Som_slice_first) || Som_slice_first < 1
            Som_slice_first = 1;
        end
        Som_slice_end = max(find(Seg_slices_info.Som_slices_segmented(:,frame-1)))+2;
        if isempty(Som_slice_end) || Som_slice_end > size(I_map,3)
            Som_slice_end = size(I_map,3);
        end
    end
    
    %% Segment each slice
    
    Noto_inc = 0;
    Vertical_extent = false;
    for slice = Noto_slice_first:Noto_slice_end

        I = I_map(:,:,slice); % Get the first slice

        I8bit = im2double(255 * mat2gray(I)); % Convert to 8-bit and double class
        figure(1);clf;imagesc(I8bit)
        set(gcf,'Position',[102 294 1770 990]) % these parameters will most likely have to be changed based on the screen you use
        
        % If you want to detect notochord in this slice, press y else just press enter and it proceeds to the next slice
        YN = input(strcat('Do you want to detect notochord in slice ',num2str(slice),'?: \n'),'s');
        if strcmp(YN,'y') 

            Noto_inc = Noto_inc + 1;
            disp('Draw a rectangle where segmentation will be performed')
            figure(1); rect = round(getrect);
            I8bit = I8bit(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));

            Somite_seg.Rectangle_analysis(1:4,slice) = rect;
            Somite_seg.Rectangle_analysis(5,slice) = rect(1)+rect(3)-1;
            Somite_seg.Rectangle_analysis(6,slice) = rect(2)+rect(4)-1;
            
            
            YN1 = input(strcat('Do you want to determine vertical extent over which somites will be segmented later?: \n'),'s');
            if ~Vertical_extent && strcmp(YN1,'y')
                disp('Draw a rectangle where segmentation will be performed')
                figure(1); rect1 = round(getrect);
                Somite_seg.Rectangle_analysis(2,slice) = rect1(2);
                Somite_seg.Rectangle_analysis(6,slice) = rect1(2)+rect1(4)-1;
                Vertical_extent = true;
            end
            
            % Apply Frangi filter, which identifies tube-like structures in the image
            % c = 2*options.FrangiBetaTwo^2; Value of 'c' has been recommended to be a quarter of the value of the maximum intensity at the vessels of interest
            options = struct('FrangiScaleRange', [1 3], 'FrangiScaleRatio', 2, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 12, 'verbose',true,'BlackWhite',false);
            [outIm,~,~,Dxn,~,~] = FrangiFilter2D_sundar(I8bit,options);

            Dxn = rad2deg(Dxn);

            % If you want to display the eigenvectors uncomment the following lines
%             [outIm,~,~,Dxn,I_x,I_y] = FrangiFilter2D_sundar(I8bit,options);
%             figure,imshow(outIm,[])
%             hold on
%             quiver(I_x,I_y,'m')
%             set(gcf,'Color','w')
%             % If you want to zoom into the figure and save the zoomed figure, run the following command before saving
%             set(gcf,'CreateFcn','z = zoom(gcf);set(z,''ActionPostCallback'',''TMWtest'')')

            %% Notochord segmentation

            outIm_noto = outIm;
            % Discard all points with vectors that are not horizontal (notochord is usually horizontal in our images)
            outIm_noto(Dxn > -180 & Dxn < -120) = nan;
            outIm_noto(Dxn > -60 & Dxn < 60) = nan;
            outIm_noto(Dxn > 120 & Dxn < 180) = nan;
            
            % Discard border on the left and right
            outIm_noto([1:20 end-20:end],:) = nan;

            Filtered_image_noto = false(size(outIm_noto));
            % Run through the image column by column and identify peaks in fluorescence intensity
            for i = 1:size(outIm_noto,2)
                Line_data = outIm_noto(:,i);
                [~,Int_points] = findpeaks(Line_data,'minPeakHeight',nanmean(Line_data)); % Only bright pixels with a minimum intensity that equals mean intensity are chosen
                Filtered_image_noto(Int_points,i) = 1;
            end

            Filtered_image_open = bwareaopen(Filtered_image_noto,2); % Discard isolated points

            E = bwmorph(Filtered_image_open,'EndPoints'); % Find end points of continuous lines
            [row,col] = find(E);
            Ecr = [col row]; % in xy
            Ecr_pdist = squareform(pdist(Ecr)); % Get distance between end points of neighbouring lines
            Ecr_pdist(Ecr_pdist == 0) = nan; % Remove diagonal elements, which will have distance between a point and itself

            figure(2);clf;imshow(Filtered_image_open)

            R = regionprops(Filtered_image_open,'PixelList'); % Get pixel coordinates (in xy) for each isolated structure
            PixelList = struct2cell(R);

            for i = 1:numel(R)

                Segment = PixelList{i}; % Choose first line
                [~,sort_segment_idx] = sort(Segment(:,1)); % Arrange the pixels from left to right
                Segment = Segment(sort_segment_idx,:);

                [~,Loc_right] = ismember(Segment(end,:),Ecr,'rows'); % Identify the right-most poistion of the line in Ecr
                [~,Loc_left] = ismember(Segment(1,:),Ecr,'rows'); % Identify the left-most position of the line in Ecr

                if Loc_left ~=0 && Loc_right ~= 0 % This 'if' statement is not really needed
                    % Remove distance between end points from the same line
                    Ecr_pdist(Loc_right,Loc_left) = nan;
                    Ecr_pdist(Loc_left,Loc_right) = nan;
                    
                    % The idea is to join lines from left to right. Therefore find the closest line to the right-most end point of the current line
                    [Min_dist,idx1] = min(Ecr_pdist(Loc_right,:));

                    while Min_dist < 15 % This is an arbitrary number, which works best for our images. Connect those lines that are closer than 15 pixels.

                        Ecr_int = Ecr(idx1,:); % Choose the closest line and get pixel coordinates
                        Segment_near_idx = find(cellfun(@(x) ismember(Ecr_int, x,'rows'),PixelList));
                        Segment_near = PixelList{Segment_near_idx};
                        [~,sort_segment_idx] = sort(Segment_near(:,1));
                        Segment_near = Segment_near(sort_segment_idx,:);
                        
                        Segment_combine = [Segment;Segment_near]; % Combine the pixel coord of the two lines
                        
                        % Check for a few cases and choose whether to consider this line or not
                        % 1. If the end point (rather than the first point) of the chosen line is more closer to the end point of the previous line, then skip
                        % 2. If the new line is greater than five pixels away in y, then skip
                        % 3. If the two lines are overlapping in x, then skip
                        if ismember(Ecr_int,Segment_near(end,:),'rows') || abs(Segment(end,2) - Ecr_int(1,2)) > 5 || any(sum(Segment_combine(:,1) == Segment_combine(:,1)') > 1)
                            Ecr_pdist([Loc_left Loc_right],idx1) = nan;
                            Ecr_pdist(idx1,[Loc_left Loc_right]) = nan;
                            [Min_dist,idx1] = min(Ecr_pdist(Loc_right,:)); % Go to the next closest line.
                            continue;
                        else
                            figure(2);
                            hline = imline(gca,[Segment(end,:);Ecr_int]); % Link the two lines
                            binaryImage2 = hline.createMask();
                            Filtered_image_open(binaryImage2) = 1; % Store the lines (including the link) in a new image

                            Ecr_pdist(:,[Loc_left Loc_right]) = nan;  % For the next iteration, remove the end points of the current line
                            Ecr_pdist([Loc_left Loc_right],:) = nan;
                            Ecr_pdist(:,idx1) = nan; % For the next iteration, remove the point from new line that was used for linking
                            Ecr_pdist(idx1,:) = nan;
                            break;
                        end
                    end
                end
            end

            R_label = regionprops(Filtered_image_open,'extrema'); % Get the extrema coord of each line

            A = outIm;
            % Some normalization to better display lines
            A = A - 1E-8;
            A = A / 0.5;
            Fperim = bwperim(Filtered_image_open);
            Ioverlay = imoverlay(A,Fperim,[1 0 0]);

            figure(3);clf;
            imshow(Ioverlay);set(gcf,'Position',[36 212 1490 1100])
            for k = 1:numel(R_label)
               e = R_label(k).Extrema; % Use the extrema coord to label the different lines
               text(e(1,1), e(1,2)-3, sprintf('%d', k),'Color','red');
            end
            
            % Get the user to enter which lines constitute the upper and lower parts of the notochord. Enter in square brackets
            Boundary_img_temp = zeros(size(outIm));
            prompt = 'Which lines constitute the upper part of the notochord?';
            IP1 = input(prompt);
            Somite_seg.Noto_up_identity(1:numel(IP1),slice) = IP1;
            if ~isempty(IP1)
                Boundary_img_temp1 = ismember(bwlabel(Filtered_image_open),IP1);
            end

            prompt = 'Which lines constitute the lower part of the notochord?';
            IP2 = input(prompt);
            Somite_seg.Noto_down_identity(1:numel(IP2),slice) = IP2;
            if ~isempty(IP2)
                Boundary_img_temp2 = ismember(bwlabel(Filtered_image_open),IP2);
            end

            figure(2);clf;imshow(Boundary_img_temp1)
            R_noto = regionprops(logical(Boundary_img_temp1),'PixelList');
            PixelList_noto = struct2cell(R_noto);

            % Link all the lines in the upper part of the notochord one after the other
            if length(R_noto) > 1
                for i = 1:length(R_noto)-1

                    Segment = PixelList_noto{i};
                    [~,sort_segment_idx] = sort(Segment(:,1));
                    Segment = Segment(sort_segment_idx,:);

                    Segment_next = PixelList_noto{i+1};
                    [~,sort_segment_idx] = sort(Segment_next(:,1));
                    Segment_next = Segment_next(sort_segment_idx,:);

                    figure(2);
                    hline = imline(gca,[Segment(end,:);Segment_next(1,:)]);
                    binaryImage2 = hline.createMask();
                    Boundary_img_temp(binaryImage2) = 1;
                end
            end

            figure(2);clf;imshow(Boundary_img_temp2)
            R_noto = regionprops(logical(Boundary_img_temp2),'PixelList');
            PixelList_noto = struct2cell(R_noto);

            % Link all the lines in the lower part of the notochord one after the other
            if length(R_noto) > 1
                for i = 1:length(R_noto)-1

                    Segment = PixelList_noto{i};
                    [~,sort_segment_idx] = sort(Segment(:,1));
                    Segment = Segment(sort_segment_idx,:);


                    Segment_next = PixelList_noto{i+1};
                    [~,sort_segment_idx] = sort(Segment_next(:,1));
                    Segment_next = Segment_next(sort_segment_idx,:);

                    figure(2);
                    hline = imline(gca,[Segment(end,:);Segment_next(1,:)]);
                    binaryImage2 = hline.createMask();
                    Boundary_img_temp(binaryImage2) = 1;
                end
            end
            Noto = false(size(I8bit));
            Noto(Boundary_img_temp == 1 | Boundary_img_temp1 == 1 | Boundary_img_temp2 == 1) = 1; % Get all the lines into one image

            Somite_seg.Noto(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),Noto_inc) = Noto; % Bring the segmented image to original size
            Seg_slices_info.Noto_slices_segmented(slice,frame) = true;
        end
    end
    
    % Combine notochord segmentation from different slices of the same frame
    Noto_all = false(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    for ii = 1:size(Somite_seg.Noto,3)
        Noto_all(Somite_seg.Noto(:,:,ii) == 1) = 1;
    end
    [~,xmin] = min(Somite_seg.Rectangle_analysis(1,:)); xmin = Somite_seg.Rectangle_analysis(1,xmin);
    [~,ymin] = min(Somite_seg.Rectangle_analysis(2,:)); ymin = Somite_seg.Rectangle_analysis(2,ymin);
    [~,xmax] = max(Somite_seg.Rectangle_analysis(5,:)); xmax = Somite_seg.Rectangle_analysis(5,xmax);
    [~,ymax] = max(Somite_seg.Rectangle_analysis(6,:)); ymax = Somite_seg.Rectangle_analysis(6,ymax);
    Noto_all = Noto_all(ymin:ymax,xmin:xmax);

    A = bwdist(Noto_all) <= 5;
    A_label = bwlabel(A);
    
    % Extend line from first column of image to closest segmented line and from last column of image to closest segmented line
    for ii = 1:2
        B = ismember(A_label,ii); % Choose the first set of lines
        B(~Noto_all) = 0; % Get rid of added pixels

        [row,col] = find(B == 1,1);
        [row1,col1] = find(B == 1,1,'last');

        figure(2);clf;imshow(Noto_all)
        hline = imline(gca,[1 row; col row]);
        binaryImage2 = hline.createMask();
        Noto_all(binaryImage2 == 1) = 1;

        hline = imline(gca,[size(Noto_all,2) row1; col1 row1]);
        binaryImage2 = hline.createMask();
        Noto_all(binaryImage2 == 1) = 1;
    end
    
    % Add vertical lines in the first and last column of image, whose vertical extents are determined by the segmented notochord
    Noto_all(find(Noto_all(:,1),1):find(Noto_all(:,1),1,'last'),1) = 1;
    Noto_all(find(Noto_all(:,end),1):find(Noto_all(:,end),1,'last'),end) = 1;

    Noto_all = bwfill(Noto_all,'holes');

    % Ask the user if the notochord is fine, else dilate it up a bit and ask again.
    % If you are happy with the notochord, just press enter, else press n
    figure(2);clf;imshow(Noto_all)
    Q = input(strcat('Are you happy with the notochord?: \n'),'s'); 
    Add_bwdist = 1;
    while strcmp(Q,'n')
        Add_bwdist = Add_bwdist + 1;
        Noto_all = bwdist(Noto_all) < Add_bwdist;
        Noto_all = bwfill(Noto_all,'holes');
        figure(2);clf;imshow(Noto_all)
        Q = input(strcat('Are you happy with the notochord?: \n'),'s'); 
    end
    
    %% Somite segmentation
    
    % Get the median position of the notochord - between the upper and lower notochord positions
    Noto_median = nan(1,size(Noto_all,2));
    for ii = 1:size(Noto_all,2)
        Noto_median(ii) = round(nanmedian(find(Noto_all(:,ii) == 1)));
    end
    
    Slice_inc = 0;
    for slice = Som_slice_first:Som_slice_end
        
        I = I_map(:,:,slice); % Get the slice that will be segmented
        
        I8bit = im2double(255 * mat2gray(I)); % Convert to 8-bit and double class
        figure(1);clf;imagesc(I8bit(ymin:ymax,xmin:xmax))
        daspect([1 1 1])
%         set(gcf,'Position',[764 358 1770 990]) % Uncomment this line and change these parameters based on the screen you use if the image is too small or large
            
        % If you want to detect somite boundaries in this slice, press y
        YN = input(strcat('Do you want to detect boundaries in slice ',num2str(slice),'?: \n'),'s');
        if strcmp(YN,'y') 
            
            Slice_inc = Slice_inc + 1;
            
            I8bit = I8bit(ymin:ymax,xmin:xmax); % Choose a subset of the image that contains regions of interest already set through notochord segmentation section
            
            figure(9);clf;imagesc(I8bit)
            daspect([1 1 1])
            disp('Draw a rectangle for segmentation to be performed')
            figure(9); rect = round(getrect);
            
            % Apply Frangi filter, which identifies tube-like structures in the image
            % c = 2*options.FrangiBetaTwo^2; Value of 'c' has been recommended to be a quarter of the value of the maximum intensity at the vessels of interest
            options = struct('FrangiScaleRange', [1 3], 'FrangiScaleRatio', 2, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 12, 'verbose',true,'BlackWhite',false);
            [outIm,~,~,Dxn,~,~] = FrangiFilter2D_sundar(I8bit,options);

            Dxn = rad2deg(Dxn);
            
            outIm_somite = outIm;
            
            % Discard all points with vectors that are horizontal (somite boundaries are more or less vertical in our images)
            outIm_somite(Dxn < -70 & Dxn > -110) = nan;
            outIm_somite(Dxn < 110 & Dxn > 70) = nan;
            
            outIm_somite(:,[1:20 end-20:end]) = nan;
            outIm_somite(Noto_all == 1) = nan;
            outIm_somite([1:20 end-20:end],:) = nan;
            
            % Remove regions where segmentation is not needed. This saves time.
            outIm_somite(:,1:rect(1)) = nan;
            outIm_somite(:,rect(1)+rect(3):end) = nan;
            outIm_somite(1:rect(2),:) = nan;
            outIm_somite(rect(2)+rect(4):end,:) = nan;
            
            Filtered_image = false(size(outIm_somite));
            Filtered_image2 = false(size(outIm_somite));
            % Run through the image row by row and identify peaks in fluorescence intensity
            for i = 1:size(outIm_somite,1)

                Line_data = outIm_somite(i,:);
                [~,Int_points] = findpeaks(Line_data,'minPeakHeight',nanmean(Line_data)+nanstd(Line_data));
                [~,Int_points2] = findpeaks(Line_data,'minPeakHeight',nanmean(Line_data)-nanstd(Line_data));
                
                Filtered_image(i,Int_points) = 1;
                Filtered_image2(i,Int_points2) = 1;
            end
            
            Filtered_image_open = bwareaopen(Filtered_image,5);
            
            E = bwmorph(Filtered_image_open,'EndPoints'); % Find end points of continuous lines
            [row,col] = find(E);
            Ecr = [col row]; % in xy
            Ecr_pdist = squareform(pdist(Ecr)); % Get distance between end points of neighbouring lines
            Ecr_pdist(Ecr_pdist == 0) = nan; % Remove diagonal elements, which will have distance between a point and itself

            figure(2);clf;imshow(Filtered_image_open)

            R = regionprops(Filtered_image_open,'PixelList'); % Get pixel coordinates (in xy) for each isolated structure
            PixelList = struct2cell(R);
            Line_bet_segs = [nan nan];
            
            for j = 1:length(Ecr)
                Ecrj_pdist = Ecr_pdist;
                
                Segment_idx = find(cellfun(@(x) ismember(Ecr(j,:), x,'rows'),PixelList)); % Get the index of end point of interest
                
                if ~isempty(Segment_idx)
                    Segment = PixelList{Segment_idx}; % Get the pixel list of chosen segment
                    [~,sort_segment_idx] = sort(Segment(:,2)); % sort the pixel indices in ascending order
                    Segment = Segment(sort_segment_idx,:);
                    [~,Loc] = ismember(Ecr(j,:),Segment,'rows'); % Find the location of the end point in Segment. Is it at the start or end?
                    
                    if Loc == 1
                        Ecr_same_seg = Segment(end,:); % If it is at the start, then choose the end point of the same segment
                    else
                        Ecr_same_seg = Segment(1,:); % If it is at the end, choose the start point of the same segment
                    end
                    Seg_int = Segment;
                end
                
                % Find the location of the old end point in Ecr_same_seg and make it nan, so that it will not be used further
                Ecr_same_seg_row = find(ismember(Ecr,Ecr_same_seg,'rows')); 
                Ecrj_pdist(Ecr_same_seg_row,j) = nan;

                [Min_dist,idx1] = min(Ecrj_pdist(:,j)); % Find the closest end point 

                if Min_dist < 15 % This is an arbitrary number, which worked well for our images. Feel free to play around this number

                    Secondary_dist = Ecrj_pdist(Ecr_same_seg_row,idx1);
                    if Secondary_dist < Min_dist
                        Ecr1 = Ecr(Ecr_same_seg_row,:);
                    else
                        Ecr1 = Ecr(j,:);
                    end

                    Ecr2 = Ecr(idx1,:); 
                    Segment2_idx = find(cellfun(@(x) ismember(Ecr2, x,'rows'),PixelList));
                    if ~isempty(Segment2_idx)
                        Segment = PixelList{Segment2_idx};
                        Seg_int = [Seg_int;Segment];
                    end
                     
                    Closest_noto = median(find(Noto_all(:,round(nanmean(Seg_int(:,1)))) == 1)) + ymin;
                    
                    Angle_line = atan2(Ecr1(2)-Ecr2(2),Ecr1(1)-Ecr2(1))*180/pi; % Angle between two lines that will be connected
                    
                    % if the angle is within a certain range, combine the two lines. Check if the lines of interest are below or above the notochord as the angle will change accordingly
                    if nanmean(Seg_int(:,2)) >= Closest_noto && ((Angle_line >= -70 && Angle_line <= 30) || (Angle_line >= 110 && Angle_line <= 180) || (Angle_line >= -180 && Angle_line <= -165))
                         
                        % Check if the two lines are at overlapping levels, if yes, skip the line and go to next 'for' loop
                        N = sum(Seg_int(:,2) == Seg_int(:,2)');
                        if any(N > 1)
                            continue;
                        end
                        
                        % Combine the two lines
                        xx = round(linspace(Ecr1(1),Ecr2(1),ceil(Min_dist)));
                        yy = round(linspace(Ecr1(2),Ecr2(2),ceil(Min_dist)));

                        Filtered_image_open(sub2ind(size(Filtered_image_open),yy,xx)) = 1;
                        
                    elseif nanmean(Seg_int(:,2)) <= Closest_noto && ((Angle_line >= -30 && Angle_line <= 70) || (Angle_line >= -180 && Angle_line <= -110) || (Angle_line >= 165 && Angle_line <= 180))
                        
                        % Check if the two lines are at overlapping levels, if yes, skip the line and go to next 'for' loop
                        N = sum(Seg_int(:,2) == Seg_int(:,2)');
                        if any(N > 1)
                            continue;
                        end
                        
                        % Combine the two lines
                        xx = round(linspace(Ecr1(1),Ecr2(1),ceil(Min_dist)));
                        yy = round(linspace(Ecr1(2),Ecr2(2),ceil(Min_dist)));

                        Filtered_image_open(sub2ind(size(Filtered_image_open),yy,xx)) = 1;
                    end
                end
            end
            Filtered_image_open2 = bwareaopen(Filtered_image_open,25);
            
            R_label = regionprops(Filtered_image_open2,'extrema');
            
            A = outIm;
            % Normalize the image a little bit for bettter visualization
            A = A - 1E-8;
            A = A / 0.5;
            Fperim = bwperim(Filtered_image_open2);
            Ioverlay = imoverlay(A,Fperim,[1 0 0]);

            figure(8);clf;imshow(I8bit,[])
            for k = 1:numel(R_label)
               e = R_label(k).Extrema;
               text(e(1,1)-3, e(1,2), sprintf('%d', k),'Color','red');
            end
            
            figure(1);shg;
            figure(3);clf;imshow(Ioverlay)
            for k = 1:numel(R_label)
               e = R_label(k).Extrema;
               text(e(1,1)-3, e(1,2), sprintf('%d', k),'Color','red');
            end
            
            % You have to decide which lines constitute a boundary and enter
            prompt = 'Which lines constitute a boundary?';
            IP = input(prompt);
            Somite_seg.Boundary_identity(1:numel(IP),slice) = IP;
            if ~isempty(IP)
                Boundary_img_temp = ismember(bwlabel(Filtered_image_open2),IP);
            end
            
            Boundary_img_temp2 = Boundary_img_temp;
           
            % Go to filtered_image2, which will have more detected lines - Check line 431
            Filtered_image2_open = bwareaopen(Filtered_image2,5);
            R3 = regionprops(logical(Filtered_image2_open),'PixelList');

            PixelList_R3 = struct2cell(R3);
            E3 = bwmorph(Filtered_image2_open,'EndPoints');

            [row3,col3] = find(E3);
            Ecr_R3 = [col3 row3]; % in xy
            
            R4 = regionprops(logical(Boundary_img_temp2),'PixelList');
            PixelList_R4 = struct2cell(R4);
            E4 = bwmorph(Boundary_img_temp2,'EndPoints');
            [row4,col4] = find(E4);
            Ecr_R4 = [col4 row4]; % in xy
            
            % Find all end points that are already part of PixelList_R4 and make them invalid
            for i = 1:size(Ecr_R3,1)
                Common = find(cellfun(@(x) ismember(Ecr_R3(i,:), x,'rows'),PixelList_R4));
                if ~isempty(Common)
                    Ecr_R3(i,:) = nan;
                end
            end
            
            Not_common = Ecr_R4(~ismember(Ecr_R4,Ecr_R3,'rows'),:);
            Ecr_R3(end+1:end+size(Not_common,1),:) = Not_common;
            
            Ecr_combine = Ecr_R3;
            Ecr_combine_pdist = squareform(pdist(Ecr_combine));
            Ecr_combine_pdist(Ecr_combine_pdist == 0) = nan;
            figure(6);clf;imshow(Boundary_img_temp2)
            
            % Extend boundaries in decreasing 'y' direction
            for bdy = 1:numel(R4)
                Ecr_combine_pdist_temp = Ecr_combine_pdist;
                
                % Choose the somite and sort it in ascending 'y'
                Somite = PixelList_R4{bdy};
                [~,sort_idx] = sort(Somite(:,2));
                Somite = Somite(sort_idx,:);
                
                % Choose top point of boundary and find its position in Ecr_combine
                Line_pt1 = Somite(1,:);
                [~,Position_top] = ismember(Line_pt1,Ecr_combine,'rows');
                
                % Find position of bot point of boundary in Ecr_combine and 
                % make invalid the distance between Position_top and Position_bot in Ecr_combine_pdist_temp
                [~,Position_bot] = ismember(Somite(end,:),Ecr_combine,'rows');
                
                if Position_top == 0 || Position_bot == 0
                    continue;
                else
                    Ecr_combine_pdist_temp(Position_top,Position_bot) = nan;

                    % Find closest end point to Position_top
                    [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position_top,:));
                    % Store positions of starting and ending point of boundary, which will be later used to prevent drawing further lines to the same positions
                    For_nan = [Position_top Position_bot];

                    % As long as min_dist is less than 5 continue connecting lines
                    while Min_dist < 5
                        % Get (x,y) of closest point and find which section in PixelList_R3 corresponds to this point
                        Ecr_new2 = Ecr_combine(idx,:);
                        Section_idx = find(cellfun(@(x) ismember(Ecr_new2, x,'rows'),PixelList_R3));

                        % If the closest point is above or at the same position as end point (Position_top) and 
                        % if it is indeed present in one of the sections, then proceed
                        if Ecr_new2(1,2) <= Line_pt1(1,2) && ~isempty(Section_idx)

                            Section = PixelList_R3{Section_idx};
                            [~,sort_idx1] = sort(Section(:,2));
                            Section = Section(sort_idx1,:);

                            xx = round(linspace(Line_pt1(1),Ecr_new2(1),ceil(Min_dist)));
                            yy = round(linspace(Line_pt1(2),Ecr_new2(2),ceil(Min_dist)));

                            Boundary_img_temp2(sub2ind(size(Boundary_img_temp2),yy,xx)) = 1;

                            % Also, take the section points and include it in Boundary_img_temp2
                            Boundary_img_temp2(sub2ind(size(Boundary_img_temp2),Section(:,2),Section(:,1))) = 1;

                            % Time to proceed to next point, which will be top point of section
                            Line_pt1 = Section(1,:);
                            [~,Position1] = ismember(Line_pt1,Ecr_combine,'rows');

                            % If the top point is already part of PixelList_R4, Position1 will be zero
                            % In that scenario, extension of this somite boundary is over and therefore break the 'while'
                            % statement and proceed to next boundary
                            if Position1 == 0
                                break;
                            else
                                % If the line can still be extended, before finding closest point make invalid all points that were already considered
                                For_nan = [For_nan idx Position1];
                                Ecr_combine_pdist_temp(Position1,For_nan) = nan;
                                % Find closest point
                                [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position1,:));
                            end
                        else
                            % If closest point is below the end point (Position_top), then make it invalid and proceed to next closest point
                            For_nan = [For_nan idx];
                            Ecr_combine_pdist_temp(Position_top,For_nan) = nan;
                            % Position_top is still interesting and has to be considered until no points less than 5 pixels remain
                            Position1 = Position_top;
                            % Find closest point
                            [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position1,:));
                        end
                    end
                end
            end % for bdy = 1:numel(R4)
            
            % Go in opposite direction now
            % Extend boundaries in increasing 'y' direction
            for bdy = 1:numel(R4)
                Ecr_combine_pdist_temp = Ecr_combine_pdist;
                
                % Choose the somite and sort it in ascending 'y'
                Somite = PixelList_R4{bdy};
                [~,sort_idx] = sort(Somite(:,2));
                Somite = Somite(sort_idx,:);
                
                % Choose bottom point of boundary and find its position in Ecr_combine
                Line_pt1 = Somite(end,:);
                [~,Position_bot] = ismember(Line_pt1,Ecr_combine,'rows');
                
                % Find position of top point of boundary in Ecr_combine and 
                % make invalid the distance between Position_bot and Position_top in Ecr_combine_pdist_temp
                [~,Position_top] = ismember(Somite(1,:),Ecr_combine,'rows'); 
                
                if Position_top == 0 || Position_bot == 0
                    continue;
                else
                    Ecr_combine_pdist_temp(Position_bot,Position_top) = nan;

                    % Find closest end point to Position_bot
                    [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position_bot,:));
                    % Store positions of starting and ending point of boundary, which will be later used to prevent drawing further lines to the same positions
                    For_nan = [Position_bot Position_top];

                    % As long as min_dist is less than 5 continue connecting lines
                    while Min_dist < 5 
                        % Get (x,y) of closest point and find which section in PixelList_R3 corresponds to this point
                        Ecr_new2 = Ecr_combine(idx,:); 
                        Section_idx = find(cellfun(@(x) ismember(Ecr_new2, x,'rows'),PixelList_R3));

                        % If the closest point is below or at the same position as end point (Position_bot) and 
                        % if it is indeed present in one of the sections, then proceed
                        if Ecr_new2(1,2) >= Line_pt1(1,2) && ~isempty(Section_idx)

                            Section = PixelList_R3{Section_idx};
                            [~,sort_idx1] = sort(Section(:,2));
                            Section = Section(sort_idx1,:);

                            xx = round(linspace(Line_pt1(1),Ecr_new2(1),ceil(Min_dist)));
                            yy = round(linspace(Line_pt1(2),Ecr_new2(2),ceil(Min_dist)));

                            Boundary_img_temp2(sub2ind(size(Boundary_img_temp2),yy,xx)) = 1;

                            % Also, take the section points and include it in Boundary_img_temp2
                            Boundary_img_temp2(sub2ind(size(Boundary_img_temp2),Section(:,2),Section(:,1))) = 1;

                            % Time to proceed to next point, which will be bottom point of section
                            Line_pt1 = Section(end,:);
                            [~,Position1] = ismember(Line_pt1,Ecr_R3,'rows');

                            % If the bottom point is already part of PixelList_R4, Position1 will be zero
                            % In that scenario, extension of this somite boundary is over and therefore break the 'while'
                            % statement and proceed to next boundary
                            if Position1 == 0
                                break;
                            else
                                % If the line can still be extended, before finding closest point make invalid all points that were already considered
                                For_nan = [For_nan idx Position1];
                                Ecr_combine_pdist_temp(Position1,For_nan) = nan;
                                % Find closest point
                                [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position1,:));
                            end
                        else
                            % If closest point is above the end point (Position_bot), then make it invalid and proceed to next closest point
                            For_nan = [For_nan idx];
                            Ecr_combine_pdist_temp(Position_bot,For_nan) = nan;
                            % Position_bot is still interesting and has to be considered until no points less than 5 pixels remain
                            Position1 = Position_bot;
                            % Find closest point
                            [Min_dist,idx] = min(Ecr_combine_pdist_temp(Position1,:));
                        end
                    end
                end
            end % for bdy = 1:numel(R4)
            
            [row,col] = find(Boundary_img_temp2 == 1);
            
            figure(7);clf;subplot(1,2,1),imagesc(I8bit)
            hold on
            plot(col,row,'LineStyle','none','Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',1)
            set(gca,'XTick','','YTick','')
            subplot(1,2,2),imagesc(I8bit)
            set(gca,'XTick','','YTick','')
            set(gcf,'Position',[150 435 2300 900],'Color','w')
            saveas(gcf,strcat(Folder_path,'/Detected_boundaries/MaxProject_Ch1_fr_',num2str(frame),'_proj_',num2str(slice)),'tif')
            
            % Assign detected boundaries in selected region of interest back to bigger image
            Somite_seg.Boundary_img(ymin:ymax,xmin:xmax,Slice_inc) = Boundary_img_temp2;
            
            % Keep track of which slices were segmented
            Seg_slices_info.Som_slices_segmented(slice,frame) = true;            
            
            I_allbdy(Somite_seg.Boundary_img(:,:,Slice_inc) == 1) = 1;
            figure(4);clf;imshow(I_allbdy)
%             set(gcf,'Position',[2894 299 1774 795]) % for two screen case
            set(gcf,'Position',[1668 551 1774 795]) % for one big screen case
            
%             % for overlaying segmentation on top of original image
%             I = mat2gray(I_map_ch1(:,:,i));
%             Seg = Somite_seg.Noto_somite(:,:,i);
%             Seg(isnan(Seg)) = 0;
%             Ioverlay = imoverlay(I,bwperim(Seg),[0 1 1]);
%             figure(1);clf;imshow(Ioverlay)
%             export_fig(gcf,strcat('Fr_26_proj_',num2str(i),'.tif'),'-q101')
            
        end % if strcmp(YN,'y')
    end % for slice
    
    save(strcat(Folder_path,'/Detected_boundaries/Somite_seg_fr_',num2str(frame),'.mat'),'Somite_seg','-v7.3')
    save(strcat(Folder_path,'/Detected_boundaries/Seg_slices_info.mat'),'Seg_slices_info')

end % for frame
