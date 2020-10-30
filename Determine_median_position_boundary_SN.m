function Determine_median_position_boundary_SN()
% This code determines median positions of segmented somite boundaries and notochord
% Pre-requisites: You should have already used the Segment_somites_notochord_SN.m code and completed segmentation
% Created by Sundar Naganathan - August 2018

% Go to the folder that contains segmented mat files
Folder_path = cd;

Somite_pixels_int = 100; % Number of pixels along a boundary to be considered

% Make a new folder called 'Median_bdy' if it does not exist already
if exist(strcat(Folder_path,'/Median_bdy'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Median_bdy'))
end

% Get info of all segmented mat files and sort them in alphabetical order
DirOutput = dir(strcat(Folder_path,'/Somite_seg_fr_*.mat'));
Mat_names = {DirOutput.name};
for i = 1:numel(Mat_names)
    Num_chars(i) = numel(Mat_names{i});
end
[~,idx] = sort(Num_chars);
Mat_names = Mat_names(idx);

% Pixel_size = 0.35; % in um
% If working on non-map projected images, comment the following two lines and uncomment the previous line
load(strcat(Folder_path,'/../Pixel_size.mat'))
load(strcat(Folder_path,'/Seg_slices_info.mat'))


for frame = 1:numel(Mat_names)
    
    load(strcat(Folder_path,'/',Mat_names{frame}))
    
%     Median_bdy.Pixel_size = Pixel_size;
    % If working on non-map projected images, comment the following 5 lines and uncomment the previous line
    Seg_slices = find(Seg_slices_info.Som_slices_segmented(:,frame));
    Median_bdy.Pixel_sizes = Pixel_size.Equator(Seg_slices(1)*2-1:Seg_slices(end)*2);
    Median_bdy.Pixel_size_mean = nanmean(Median_bdy.Pixel_sizes);
    
    Frame_identity = strrep(Mat_names{frame},'Somite_seg_fr_','');    
    Frame_identity = str2num(strrep(Frame_identity,'.mat',''));
        
    
    % Initialize variables
    Median_bdy.Somite_median_left = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Somite_median_right = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Projected_somite = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Projected_noto = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Noto_mid = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Noto_top = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy.Noto_bot = nan(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));

    Projected_somite = zeros(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Projected_noto = zeros(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    
    Median_bdy_left_temp = zeros(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    Median_bdy_right_temp = zeros(size(Somite_seg.Noto,1),size(Somite_seg.Noto,2));
    
    Median_bdy.mfile = 'Determine_median_position_boundary_SN'; % Record the name of the mfile
    
    % Project all segmented slices of notochord from a single frame on top of each other
    for slice = 1:size(Somite_seg.Noto,3)        
        Noto_slice = Somite_seg.Noto(:,:,slice);
        Projected_noto(Noto_slice == 1) = 1;
    end
    Median_bdy.Projected_noto = Projected_noto;
    
    % Project all segmented slices of somite boundaries from a single frame on top of each other
    for slice = 1:size(Somite_seg.Boundary_img,3)
        Boundary_img_slice = Somite_seg.Boundary_img(:,:,slice);
        Projected_somite(Boundary_img_slice == 1) = 1;
    end
    Median_bdy.Projected_somite = Projected_somite;

    figure(1);clf;imshow(Median_bdy.Projected_noto)
%     set(gcf,'Position',[50 50 974 774]) % change these parameters depending on the screen you use
    
    % Check whether the segmented notochord is fine i.e no obscure lines
    % If there are indeed portions of notochord to be deleted, type n and then draw a rectangle where poritions of lines have to be deleted
    YN = input(strcat('Are you happy with the notochord in frame ',num2str(frame),'?: \n'),'s');
    while strcmp(YN,'n')
        disp('Draw a rectangle where segmented notochord can be deleted')
        figure(1); rect = round(getrect);
        Median_bdy.Projected_noto(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3)) = 0;
        
        figure(1);clf;imshow(Median_bdy.Projected_noto)
%         set(gcf,'Position',[50 50 974 774]) % change these parameters depending on the screen you use
        
        YN = input(strcat('Are you happy with the notochord in frame ',num2str(frame),'?: \n'),'s');
    end
    
    % Scan through each column of the image and determine mean and mid positions of notochord
    for i = 1:size(Median_bdy_left_temp,2)
        [row_noto,~] = find(Median_bdy.Projected_noto(:,i) == 1);
        
        % If number of notochord positions is greater than 1, then proceed. 
        % If it is empty or equal to one (i.e only top of bottom half of notochord segmented), then proceed to next column
        if numel(row_noto) > 1
            row_noto_uniq = unique(row_noto);
            
            % Get all top and bottom positions based on mean of row of segmented notochord
            row_noto_top = round(mean(row_noto_uniq(row_noto_uniq < mean(row_noto_uniq))));
            row_noto_bot = round(mean(row_noto_uniq(row_noto_uniq > mean(row_noto_uniq))));
            
            % If the difference between top and bottom is less than 30 then probably only top or bottom half of notochord is segmented
            % Proceed to next column
            if abs(row_noto_bot - row_noto_top) < 30
                continue;
            end
            Median_bdy.Noto_top(row_noto_top,i) = 1;
            Median_bdy.Noto_bot(row_noto_bot,i) = 1;
            Median_bdy.Noto_mid(round(mean([row_noto_top row_noto_bot])),i) = 1; % Get middle position of notochord
        end
    end
    
    % Find all connected points that are within 5 pixel distance. This makes sure that a single boundary is not labelled multiple times.
    A = double(bwdist(Projected_somite) <= 5);
    
    [A_label,num] = bwlabel(A); % Label the swollen boundaries
    A_regionprops = regionprops(A_label,'extrema'); % Find extrema of each labelled boundary and use it for labelling.
    Identity_A = 1:numel(A_regionprops);
    
    B = labeloverlay(A,A_label);
    figure(2);clf;imshow(B)
%     set(gcf,'Position',[50 50 974 774]) % change these parameters depending on the screen you use
    hold on
    for k = 1:numel(A_regionprops)
        e = A_regionprops(k).Extrema;
        text(e(1,1)-8, e(1,2), sprintf('%d', k),'Color','red');
    end
    
    % There might be lines which do not constitute a boundary. Very rare, but better to check before proceeding.
    prompt = ['Which lines should be deleted in frame ' num2str(frame) ' ?'];
    IP = input(prompt);
    
    % If any line has to be deleted, go ahead
    if ~isempty(IP)
        for i = 1:numel(IP)
            Identity_A(Identity_A == IP(i)) = [];
        end
        A = double(ismember(A_label,Identity_A));
        [A_label,num] = bwlabel(A); % Re-label boundaries
        A_regionprops = regionprops(A_label,'extrema');
        B = labeloverlay(A,A_label);
        figure(2);clf;imshow(B)
%         set(gcf,'Position',[29 582 974 774]) % change these parameters depending on the screen you use
        hold on
        for k = 1:numel(A_regionprops)
            e = A_regionprops(k).Extrema;
            text(e(1,1)-8, e(1,2), sprintf('%d', k),'Color','red');
        end
    end
    Median_bdy.Projected_somite = Median_bdy.Projected_somite .* A; % Remove deleted lines from the to-be-saved projected_somite variable as well
        
    % There might be lines which need to be considered as a single boundary
    prompt = ['Which lines should be combined in frame ' num2str(frame) ' ?'];
    IP = input(prompt);
    Removed = [];
    
    % Link lines until no more lines have to be dealt with
    while(~isempty(IP))
        for i = 2:numel(IP)
            A_label(A_label == IP(i)) = IP(1);
            Removed = [Removed;IP(i)];
        end
        prompt = ['Which lines should be combined in frame ' num2str(frame) ' ?'];
        IP = input(prompt);
    end
    
    % Loop through each boundary and determine median position
    for bdy = 1:num
        if ~ismember(bdy,Removed)
            B = ismember(A_label,bdy); % Choose a single boundary
            B(~Projected_somite) = 0; % Get rid of added pixels
            
            % Multiply with projected image as the projected image contains weights
            B = B .* Projected_somite;
            [row,col] = find(B > 0); % Find all positions which correspond to a boundary
            
            % In case of a weight of more than 1, represent these positions more times
            for i = 1:numel(row)
                Value = B(row(i),col(i));
                if Value > 1
                    for j = 2:Value
                        row(end+1) = row(i);
                        col(end+1) = col(i);
                    end
                end
            end
            
            % Determine median row, which will be later used to determine whether the somite boundary is to the right (top) or left (bottom) of notochord
            row_uniq = unique(row); 
            row_mid = round(median(row_uniq));
            col_mid = round(median(col));
            
            % Determine notochord position along the median column of a boundary
            Noto_row = find(Median_bdy.Noto_mid(:,col_mid) == 1);

            % If no notochord has been segmented in the same column as col_mid, then find the nearest notochord segmented position
            if isempty(Noto_row)
                % Search for closest notochord to the right of col_mid
                for ii = col_mid+1:size(Median_bdy.Noto_mid,2)
                    [Noto_row1,~] = find(Median_bdy.Noto_mid(:,ii) == 1);
                    if ~isempty(Noto_row1)
                        break;
                    end
                end

                % Search for closest notochord to the left of col_mid
                for jj = 1:col_mid-1
                    jj_change = col_mid - jj;
                    [Noto_row2,~] = find(Median_bdy.Noto_mid(:,jj_change) == 1);
                    if ~isempty(Noto_row2)
                        break;
                    end
                end

                % Find which is more closer and use that position
                if min(abs(col_mid - mean(Noto_row1)),abs(col_mid - mean(Noto_row2))) == abs(col_mid - mean(Noto_row1))
                    Noto_row = Noto_row1;
                else
                    Noto_row = Noto_row2;
                end
            end
            
            % Check if somite boundary is above or below and proceed
            % If median row of boundary is less than notochord position, then it is above (left), else it is below (right)
            % Above is left in our case as we have a dorsal view
            if row_mid < Noto_row
                
                % Find the end row of boundary towards the notochord
                Bot_most_row_bdy = find(row == max(row));
                Noto_top_row = find(Median_bdy.Noto_top(:,col_mid) == 1);
                
                if row(Bot_most_row_bdy(1)) >= Noto_top_row
                    row_idx = find(row < Noto_top_row & row > (Noto_top_row - Somite_pixels_int));
                else
                    row_idx = find(row < row(Bot_most_row_bdy(1)) & row > (row(Bot_most_row_bdy(1)) - Somite_pixels_int));
                end
                
                % Loop through each row and determine mean position of boundary
                for i = min(row(row_idx)):max(row(row_idx))

                    col_bdy = round(nanmean(col(row == i)));
                    if ~isnan(col_bdy) && ~isempty(col_bdy)
                        Median_bdy_left_temp(i,col_bdy) = 1;
                    end
                end
             
            % The same set of considerations as above
            elseif row_mid > Noto_row    
                
                Top_most_row_bdy = find(row == min(row));
                Noto_bot_row = find(Median_bdy.Noto_bot(:,col_mid) == 1);
                
                if row(Top_most_row_bdy(1)) <= Noto_bot_row
                    row_idx = find(row > Noto_bot_row & row < (Noto_bot_row + Somite_pixels_int));
                else
                    row_idx = find(row > row(Top_most_row_bdy(1)) & row < (row(Top_most_row_bdy(1)) + Somite_pixels_int));
                end
                
                for i = min(row(row_idx)):max(row(row_idx))

                    col_bdy = round(nanmean(col(row == i)));
                    if ~isnan(col_bdy) && ~isempty(col_bdy)
                        Median_bdy_right_temp(i,col_bdy) = 1;
                    end
                end
            end
        end
    end
    Median_bdy.Somite_median_left = Median_bdy_left_temp;
    Median_bdy.Somite_median_right = Median_bdy_right_temp;
    
    % Plot all the lines in the same image
    [rowsl,colsl] = find(Median_bdy.Somite_median_left == 1);
    [rowsr,colsr] = find(Median_bdy.Somite_median_right == 1);
    [rownt,colnt] = find(Median_bdy.Noto_top == 1);
    [rownb,colnb] = find(Median_bdy.Noto_bot == 1);
    [rownm,colnm] = find(Median_bdy.Noto_mid == 1);
    
    figure(3);clf;imshow(Median_bdy.Projected_somite)
%     set(gcf,'Position',[50 50 974 774]) % change these parameters depending on the screen you use
    hold on
    plot(colsl,rowsl,'or','MarkerSize',1)
    plot(colsr,rowsr,'or','MarkerSize',1)
    plot(colnt,rownt,'oc','MarkerSize',1)
    plot(colnb,rownb,'oc','MarkerSize',1)
    plot(colnm,rownm,'or','MarkerSize',1)
    
    saveas(gcf,strcat(Folder_path,'/Median_bdy/Median_bdy_',num2str(Frame_identity)),'tif')
    
    save(strcat(Folder_path,'/Median_bdy/Median_bdy_',num2str(Frame_identity),'.mat'),'Median_bdy')
    clear Median_bdy
end
