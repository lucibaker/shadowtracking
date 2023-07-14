% calibration validation 
% do calibration with half the points (training points) and use that to
% predict the location of the other half (test points)

%% SETUP
clear
close all

gdrive_path = 'G:\My Drive\';  % 'C:\Users\ljbak\My Drive\';  %  
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set

% load experiment params
warning off
cal_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\cal_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

% image parameters 
cams = cell2mat(cal_params.Cam)';

img_shifts = zeros(4,2); 


%% LOOP OVER CAMERAS

rectify_quad = cell(length(cams),1);      % image rectification function handles for each camera view

for cam = 1%:length(cams)   
    
    cam_left = cam <= 2; % 'true' for the left two cameras, 'false' for right two cameras

    % load background and calibration images
    bkgd = cam_imread(sprintf('Cam%s-bkgd.tif',cams(cam)), cam_left);
    cal = cam_imread(sprintf('Cam%s-cal.tif',cams(cam)), cam_left);
    mask = cam_imread(sprintf('Cam%s-mask.tif',cams(cam)), cam_left);

    bkgd = double(bkgd);
    img_ix = size(bkgd,2); img_iy = size(bkgd,1);

    % crop to correct for camera shift
    bkgd_crop_rect = [abs(min([1,img_shifts(cam,1)])), max([1,img_shifts(cam,2)]), ...
        img_ix - max([0,img_shifts(cam,1)]), img_iy - abs(min([0,img_shifts(cam,2)]))];
    img_crop_rect = [max([1,img_shifts(cam,1)]), abs(min([1,img_shifts(cam,2)])), ...
        img_ix - abs(min([0,img_shifts(cam,1)])), img_iy - max([0,img_shifts(cam,2)])];
    
    bkgd = imcrop(bkgd, bkgd_crop_rect);
    cal = imcrop(cal, bkgd_crop_rect);
    mask = imcrop(mask, bkgd_crop_rect);

    img_ix = size(bkgd,2); img_iy = size(bkgd,1);



    %% CALIBRATION: GET MAPPING FUNCTION FROM IMAGE TO WORLD COORDS
    
    cal = double(cal) - bkgd; % subtract background
    cal = uint8(cal - min(cal(:)));  % shift intensities so that all are positive 
    cal = 255 - cal;  % invert image
    
    % binarize
    B = imbinarize(cal,'adaptive','Sensitivity',cal_params.calAdBinSens(cam)); % Adaptative binarization   
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); % remove false positive pixels by eroding and dilating twice
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 

    % number of rows and columns of dots in the calibration grid 
    n_rows_cal = cal_params.nRowsCal(cam); n_cols_cal = cal_params.nColsCal(cam);

    % mask top and edges (dots not part of the n_rows_cal x n_cols_cal grid)
    B = B.*logical(mask);

    figure; subplot(161); pcolor_img(cal); title('bkgd sub and inverted')
    subplot(162); pcolor_img(B); title('masked and binarized')

    % detect dots
    CC = bwconncomp(B);
    cal_dots = regionprops('table',CC,cal,'Centroid','Area','EquivDiameter');
    idx = cal_dots.Area > cal_params.calAreaThres1_px(cam) & cal_dots.Area < cal_params.calAreaThres2_px(cam);
    cal_dots = cal_dots(idx,:);    

    subplot(163); pcolor_img(cal); hold on
    viscircles(cal_dots.Centroid,cal_dots.EquivDiameter/2); title('dots detected')

    % known coords of the dots in physical space (world coords)
    W = generateCircleGridPoints([n_rows_cal, n_cols_cal], cal_params.spacing_m(cam), "PatternType","symmetric") + ...
        [cal_params.worldOffset_x(cam), cal_params.worldOffset_y(cam)]*cal_params.spacing_m(cam);

    % coords of the dots in image coords: sort dots into ascending order
    % from lower left corner of image by separating the dots by
    % x-coordinate into vertical bins   
    [~,top_row] = sort(cal_dots.Centroid(:,2),'descend'); 
    top_row = top_row(1:n_cols_cal);  % top row of dots
    bin_lim = cal_dots.Centroid(top_row,1); 
    bin_lim = sort(bin_lim,'ascend');
    bin_lim = [bin_lim - diff(bin_lim(1:2))/2; inf];  % bin limits
    
    I = zeros(n_rows_cal*n_cols_cal,2); 
    for j = 1:n_cols_cal    
        line([bin_lim(j) bin_lim(j)],[0 img_iy]) % plot bin limits
        cal_dots0 = cal_dots(cal_dots.Centroid(:,1) > bin_lim(j) & cal_dots.Centroid(:,1) < bin_lim(j+1),:);  % dots in a vertical bin
        [~,sort_idx] = sort(cal_dots0.Centroid(:,2),'ascend');  % sort by y-coord
        I(n_rows_cal*(j-1)+1 : n_rows_cal*j,:) = cal_dots0.Centroid(sort_idx,:);  % image points
    end

    point_colors = jet(size(I,1));
    subplot(164); scatter(I(:,1),I(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('detected points')
    subplot(165); scatter(W(:,1),W(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('known coords')


    % separate into training points and test points
    train_idx = false(size(I,1),1);
    train_idx(1:round(end/2)) = true;

    I_train = I(train_idx,:);
    W_train = W(train_idx,:);
    I_test = I(~train_idx,:);
    W_test = W(~train_idx,:);

    % undistortion and rectification function on training points
    rectify_quad{cam} = calibrate_camera(I_train,W_train,2);    
    
    % confirm that calibration dots are mapped to the correct world coords
    imagePoints2 = rectify_quad{cam}(I); 
    subplot(166); scatter(imagePoints2(:,1),imagePoints2(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('corrected points')

    % error on prediction of test point position
    imagePoints_test = rectify_quad{cam}(I_test);
    rms_error = rms(W_test - imagePoints_test);
    disp('rms calibration error x and y - test points [mm]')
    disp(rms_error*1000)

    % error on prediction when all points are used
    rectify_quad{cam} = calibrate_camera(I,W,2);
    imagePoints2 = rectify_quad{cam}(I);
    rms_error = rms(W - imagePoints2);
    disp('rms calibration error x and y - test points [mm]')
    disp(rms_error*1000)

end