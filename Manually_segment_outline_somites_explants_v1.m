function Manually_segment_outline_somites_explants_v1()
% This code was written by Sundar Naganathan in June 2021
% This is a very simple code that takes in the frame number, z-slice and name and intimates the user to perform a freehand outline of the explant

close all
Folder_path = cd; % Keep a note of the current folder

Frame = 38; % Type in frame number
Tif_name = "t0038_ch2.tif"; % Type in name of tif file

Somite_seg.Slice_int = 45; % Type in z-slice that has to be outlined

% Make a folder called 'Segmentation_jove' if it is not already available
if exist(strcat(Folder_path,'/Segmentation_jove'),'dir') ~= 7
    mkdir(strcat(Folder_path,'/Segmentation_jove'))
end

Somite_seg.mfile = 'Manually_segment_outline_somites_explants_v1'; % Store the name of the mfile
    
I = imread(strcat(Folder_path,'/',Tif_name),Somite_seg.Slice_int); % Read the tif file
I8bit = im2double(255 * mat2gray(I)); % Convert to 8 bit

%     figure(1);clf;imagesc(I8bit)
%     set(gcf,'Position',[390 330 1770 990]) % with two screens
%     set(gcf,'Position',[1650 310 1770 990]) % with a single big screen
%     set(gcf,'Position',[102 294 1770 990]) % with a single screen

figure(1);clf;imagesc(I8bit,[10 130])
%     Pos = get(gcf,'Position');
%     set(gcf,'Position',[469 514 Pos(3) Pos(4)])
daspect([1 1 1])
            
%         Q = input(strcat('Can you observe the features properly?: \n'),'s');
%         upd = 1;
%         while strcmp(Q,'n')
%             Q1 = input(strcat('Do you want to increase or decrease dynamic range?: \n'),'s');
%             if strcmp(Q1,'inc')
%                 figure(1);clf;imagesc(I8bit,[0 175+25*upd])
%             else
%                 figure(1);clf;imagesc(I8bit,[0 175-25*upd])
%             end
%             upd = upd + 1;
%             Q = input(strcat('Can you observe the features properly?: \n'),'s');
%         end
    
M = imfreehand(gca); % Intimate the user to draw an outline of the explant

% F = false(size(M.createMask));
Somite_seg.P0 = M.getPosition; % get xy positions
% D = round([0; cumsum(sum(abs(diff(P0)),2))]); % Measure the distance between points
% [D, idx] = unique(D);
% P = interp1(D,P0(idx),D(1):.5:D(end)); % Close the gaps
% P = unique(round(P),'rows');
% S = sub2ind(size(I),P(:,2),P(:,1));
% F(S) = true;
% 
% Somite_seg.Boundary_img = F;

% dlmwrite(strcat(Folder_path,'/Segmentation_jove/Somite_seg_fr_',num2str(Frame),'.txt'),Somite_seg.P0,'delimiter','\t','newline','pc','precision','%2.2f')
save(strcat(Folder_path,'/Segmentation_jove/Somite_seg_fr_',num2str(Frame),'_slice',num2str(Somite_seg.Slice_int),'.mat'),'Somite_seg')
