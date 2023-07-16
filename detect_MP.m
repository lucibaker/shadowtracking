function detect_MP(n)
% plastic particle detection and tracking
% input: n = run number
%
% images for each run must stored be in folders named 'run1', 'run2', etc
% image filenames must start with the camera name, e.g. 'CamA_0001.tif'
% experiment parameters must be stored in a spreadsheet starting with
%   'run_parameters', e.g. 'run_parameters_demo.xlsx'
% calibration parameters must be stored in a spreadsheet starting with
%   'cal_parameters', e.g. 'cal_parameters_demo.xlsx'
%
% Luci Baker 2023 (ljbaker36@gmail.com)

%% SETUP - set these variables first 

code_path = 'H:\My Drive\MATLAB\shadowtracking\';  % path to shadowtracking functions
cal_path = 'H:\My Drive\MATLAB\shadowtracking\calibration_and_bkgd\';   % path to calibration, background, and mask images
expt_name = 'demo';  % experiment or dataset name 

% show plots or not
plot_calibration_on = true; % calibration
plot_detection_on = true; % detection
plot_merge_on = true; % merged views

% save outputs or not
save_on = true;


%% get experiment and calibration parameters 

% add path to shadowtracking functions
addpath([code_path 'calibration functions'])       
addpath([code_path 'track functions'])

% load experiment params
warning off
run_params = readtable(sprintf('%s\\data_%s\\run_parameters_%s.xlsx', code_path, expt_name, expt_name));
cal_params = readtable(sprintf('%s\\data_%s\\cal_parameters_%s.xlsx', code_path, expt_name, expt_name));
warning on

% display the current experimental run
fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), run_params.ParticleType{n});

% are the current particles spheres or nonspheres?
nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

% image parameters 
cams = cell2mat(cal_params.Cam)';
dir_name = sprintf('run%i\\', run_params.Run(n));
imgset = dir(sprintf('%sCam%s*',dir_name,cams(1))); 
img_nt = length(imgset);

img_shifts = floor([run_params.camA_offset_x(n), run_params.camA_offset_y(n); run_params.camC_offset_x(n), run_params.camC_offset_y(n); 
    run_params.camB_offset_x(n), run_params.camB_offset_y(n); run_params.camD_offset_x(n), run_params.camD_offset_y(n)]./2);

% adaptive binarization sensitivity
ABsens = run_params.adBinSensitivity(n); % 'Sensitivity' (range [0, 1]). High value thresholds more pixels as foregrnd, at risk of including some backgrnd pixels

% particle detection parameters
A_thres = [run_params.areaThreshold1_px(n), run_params.areaThreshold2_px(n)];
maj_thres = [run_params.majAxThres1(n), run_params.majAxThres2(n)];
min_thres = [run_params.minAxThres1(n), run_params.minAxThres2(n)];

% quiescent free surface water level
z_freesurf_m = zeros(8,2);

% set up figure
if plot_calibration_on
    Bfig = figure; 
    set(Bfig,'Units', 'Pixels','Position',[103  80  1192  568]); 
    set(gcf,'color','w');
end



%% Preallocate variables

tic  % time the computation

% Find coordinates of the key geometric points of the particle shadows for each camera (in px): 
% -endpoints of major and minor axes for disks
% -endpoints of minor axis for rods
% -centroid for spheres

% preallocate keypoints arrays:
% {xp_maj1 yp_maj1}, {xp_maj2 yp_maj2}, {xp_min1 yp_min1}, {xp_min2 yp_min2}, {xp, yp}
keypts1 = cell(img_nt, length(cams));  
keypts2 = cell(img_nt, length(cams));
keypts3 = cell(img_nt, length(cams));
keypts4 = cell(img_nt, length(cams));
keypts5 = cell(img_nt, length(cams));

% preallocate free surface array
z_freesurf_inst = cell(img_nt,length(cams));

% preallocate image rectification function handles for each camera view
rectify_quad = cell(length(cams),1);      


%% LOOP OVER CAMERAS
for cam = 1:length(cams)   
    
    cam_left = (cam <= 2); % 'true' for the left two cameras, 'false' for right two cameras
    
    % read image files from the directory they are in 
    imgset = dir(sprintf('%sCam%s*', dir_name, cams(cam))); 

    % load background and calibration images
    bkgd = cam_imread(sprintf('%s\\Cam%s-bkgd.tif', cal_path, cams(cam)), cam_left);
    cal = cam_imread(sprintf('%s\\Cam%s-cal.tif', cal_path, cams(cam)), cam_left);
    mask = cam_imread(sprintf('%s\\Cam%s-mask.tif', cal_path, cams(cam)), cam_left);

    % background image must be "double" type
    bkgd = double(bkgd);

    % image dimensions
    img_ix = size(bkgd,2); img_iy = size(bkgd,1);

    % crop to correct for camera shift (if camera shifted during experiment)
    bkgd_crop_rect = [abs(min([1,img_shifts(cam,1)])), max([1,img_shifts(cam,2)]), ...
        img_ix - max([0,img_shifts(cam,1)]), img_iy - abs(min([0,img_shifts(cam,2)]))];
    img_crop_rect = [max([1,img_shifts(cam,1)]), abs(min([1,img_shifts(cam,2)])), ...
        img_ix - abs(min([0,img_shifts(cam,1)])), img_iy - max([0,img_shifts(cam,2)])];
    
    % crop background, calibration, and mask images to account for shift
    bkgd = imcrop(bkgd, bkgd_crop_rect);
    cal = imcrop(cal, bkgd_crop_rect);
    mask = imcrop(mask, bkgd_crop_rect);

    % cropped image dimensions
    img_ix = size(bkgd,2); img_iy = size(bkgd,1);



    %% CALIBRATION: GET MAPPING FUNCTION FROM IMAGE TO WORLD COORDINATES
    
    % subtract background from calibration image
    cal = double(cal) - bkgd; 

    % shift image intensities so that all are positive
    cal = uint8(cal - min(cal(:)));

    % invert image intensities
    cal = 255 - cal;  
    
    % binarize calibration image (adaptative binarization method)
    B = imbinarize(cal,'adaptive','Sensitivity', cal_params.calAdBinSens(cam));   

    % remove false positive pixels by eroding and dilating twice
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 
    B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]); 

    % number of rows and columns of dots in the calibration grid 
    n_rows_cal = cal_params.nRowsCal(cam); 
    n_cols_cal = cal_params.nColsCal(cam);

    % mask top and edges of calibration image (to remove dots that are not 
    % part of the n_rows_cal x n_cols_cal grid)
    B = B.*logical(mask);

    % view calibration images
    if plot_calibration_on
        figure; subplot(161); pcolor_img(cal); title('bkgd sub and inverted')
        subplot(162); pcolor_img(B); title('masked and binarized')
    end

    % detect calibration dots
    CC = bwconncomp(B);
    cal_dots = regionprops('table', CC, cal, 'Centroid', 'Area', 'EquivDiameter');

    % remove false detections based on dot size
    idx = cal_dots.Area > cal_params.calAreaThres1_px(cam) & cal_dots.Area < cal_params.calAreaThres2_px(cam);
    cal_dots = cal_dots(idx,:);    

    % view detected dots
    if plot_calibration_on
        subplot(163); pcolor_img(cal); hold on
        viscircles(cal_dots.Centroid,cal_dots.EquivDiameter/2); title('dots detected')
    end

    % known coords of the dots in physical space (world coordinates):
    % generate grid
    W = generateCircleGridPoints([n_rows_cal, n_cols_cal], cal_params.spacing_m(cam), "PatternType","symmetric") + ...
        [cal_params.worldOffset_x(cam), cal_params.worldOffset_y(cam)]*cal_params.spacing_m(cam);

    % corresponding coordinates of detected dots in the image (image coordinates): 
    % sort dots into ascending order from lower left corner of image by
    % separating the dots by x-coordinate into vertical bins
    [~,top_row] = sort(cal_dots.Centroid(:,2),'descend'); 
    top_row = top_row(1:n_cols_cal);  % top row of dots
    bin_lim = cal_dots.Centroid(top_row,1); 
    bin_lim = sort(bin_lim,'ascend');
    bin_lim = [bin_lim - diff(bin_lim(1:2))/2; inf];  % bin limits
    
    I = zeros(n_rows_cal*n_cols_cal,2); 
    for j = 1:n_cols_cal    
        if plot_calibration_on; line([bin_lim(j) bin_lim(j)],[0 img_iy]); end  % plot bin limits
        cal_dots0 = cal_dots(cal_dots.Centroid(:,1) > bin_lim(j) & cal_dots.Centroid(:,1) < bin_lim(j+1),:);  % dots in a vertical bin
        [~,sort_idx] = sort(cal_dots0.Centroid(:,2),'ascend');  % sort by y-coord
        I(n_rows_cal*(j-1)+1 : n_rows_cal*j,:) = cal_dots0.Centroid(sort_idx,:);  % image points
    end
    if plot_calibration_on  % plot detected image points
        point_colors = jet(size(I,1));
        subplot(164); scatter(I(:,1),I(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('detected points')
        subplot(165); scatter(W(:,1),W(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('known coords')
    end

    % CALIBRATION: compute the image rectification function for this camera
    % from the world coords and image coords of the calibration dots (using
    % a quadratic transformation: order = 2)
    trans_order = 2;
    rectify_quad{cam} = calibrate_camera(I,W,trans_order);    
    
    % view calibration to confirm that calibration dots are mapped to the correct world coords
    if plot_calibration_on
        imagePoints2 = rectify_quad{cam}(I); 
        subplot(166); scatter(imagePoints2(:,1),imagePoints2(:,2),15,point_colors,"filled"); 
        grid on; axis tight equal; title('corrected points')
    end

    % transform quiescent free surface from image to world coordinates
    z_freesurf_m(2*(cam-1)+1:2*cam,:) = rectify_quad{cam}([cal_params.xSurfStill1_px(cam), cal_params.ySurfStill1_px(cam); ...
        cal_params.xSurfStill2_px(cam), cal_params.ySurfStill2_px(cam)]);

    % view calibration
    if plot_calibration_on; figure(Bfig); end


    %% LOOP OVER FRAMES
    i0 = 1;             % starting frame 
    nframes = img_nt;   % number of frames to process
    
    for i = i0:i0+nframes-1 % can parallelize this loop using parfor instead of for

        % load the frame
        A = cam_imread([dir_name imgset(i).name], cam_left);

        A0 = imcrop(A, img_crop_rect);  % crop to correct for camera shift
        A0 = double(A0) - bkgd; % subtract background
        A0 = uint8(A0 - min(A0(:)));  % shift intensities so that all are positive 
        A0 = 255 - A0;  % invert image 
    
        % adaptive binarization
        B = imbinarize(A0,'adaptive','Sensitivity',ABsens); % Adaptative binarization   
        B = imerode(B,[0 1 0; 1 1 1; 0 1 0]); % remove false positive pixels
        B = imerode(B,[0 1 0; 1 1 1; 0 1 0]);
        B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]);
        B = imdilate(B,[0 1 0; 1 1 1; 0 1 0]);


        % FIND PARTICLES
        S1 = regionprops('table',B,A0,'Centroid','Area','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','PixelList');
        
        % free surface 
        [~,idx] = max(S1.Area);
        S_freesurf = S1(idx,:);
        z_freesurf_px = zeros(img_ix,1);
        for j = 1:img_ix
            freesurf_col_j = S_freesurf.PixelList{1}(:,1)==j;
            if any(freesurf_col_j)
                z_freesurf_px(j) = min(S_freesurf.PixelList{1}(freesurf_col_j,2));
            else
                z_freesurf_px(j) = inf;
            end
        end
        z_freesurf_inst{i,cam} = [(1:length(z_freesurf_px))', z_freesurf_px];
        
        xp = []; yp = [];
        if ~isempty(S1)
            % remove based on area and proximity to free surface
            idx = S1.Area > A_thres(1) & S1.Area < A_thres(2) & S1.Centroid(:,2) < z_freesurf_px(round(S1.Centroid(:,1))) & ...
                S1.MajorAxisLength > maj_thres(1) & S1.MajorAxisLength < maj_thres(2) & S1.MinorAxisLength > min_thres(1) & S1.MinorAxisLength < min_thres(2); 
            S = S1(idx,:);

            if ~isempty(S)
                % remove doubles (false shadows when the particle that shows up in addition to its shadow 
                % if particle is near the back wall) 
                dx_double = 300; % max horizontal distance between doubles
                dx_double_l = dx_double*~cam_left; % if particle is to the left of the shadow/right camera (px)
                dx_double_r = -dx_double*cam_left; % if particle is to the right of the shadow/left camera (px)
                dz_double = 20; % max vertical distance between doubles (px)
                
                Np = height(S);
                xp = S.Centroid(:,1);
                yp = S.Centroid(:,2);
                distx = repmat(xp,1,Np) - repmat(xp',Np,1);
                disty = repmat(yp,1,Np) - repmat(yp',Np,1);
                [~,double_idx] = find(distx ~= 0 & disty ~= 0 & distx < dx_double_l & distx > dx_double_r & disty > -dz_double & disty < dz_double);
                if ~isempty(double_idx)
                    S(double_idx,:) = [];
                end
                
                % particle centroids and orientations
                Np = height(S);
                xp = S.Centroid(:,1);
                yp = S.Centroid(:,2); 

                % rods
                if strncmp(run_params.ParticleType{n},'r',1)  
                    len1 = sqrt(S.BoundingBox(:,3).^2 + S.BoundingBox(:,4).^2);  % length
                    len2 = S.Area./len1; % thickness
                    len1 = len1 - 0.25*len2;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need the (-S.Orientation) because images are flipped up-down)
                    d_p = len1; % apparent major axis length [px]
                
                % disks and spheres
                else                                            
                    len1 = S.MajorAxisLength;
                    len2 = S.MinorAxisLength;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need the (-S.Orientation) because images are flipped up-down)
                    d_p = len2; % apparent minor axis length [px]
                end
            else
                Np = 0;
            end
        else
            Np = 0;
        end
    
        % accumulate centroids [xp, yp], angles, and keypoints (px)
        if Np > 0
            keypts1{i,cam} = [xp - len1/2.*cos(th_p), yp - len1/2.*sin(th_p)];
            keypts2{i,cam} = [xp + len1/2.*cos(th_p), yp + len1/2.*sin(th_p)];
            keypts3{i,cam} = [xp - len2/2.*cos(th_p-pi/2), yp - len2/2.*sin(th_p-pi/2)]; 
            keypts4{i,cam} = [xp + len2/2.*cos(th_p-pi/2), yp + len2/2.*sin(th_p-pi/2)];
            keypts5{i,cam} = [xp, yp];
        else
            keypts1{i,cam} = zeros(0,2);
            keypts2{i,cam} = zeros(0,2);
            keypts3{i,cam} = zeros(0,2); 
            keypts4{i,cam} = zeros(0,2);
            keypts5{i,cam} = zeros(0,2);
        end
    
        % print progress
        if ~mod(i-1,1000)
            fprintf([num2str(i-1) '/' num2str(img_nt) '\n'])
        end
    
        % make plots
        if plot_detection_on
            cla; pcolor_img(A0); hold on
            if ~isempty(xp)
                % plot centroid
                plot(xp,yp,'r.','markersize',4'); 

                % plot major and minor axes
                if strncmp(run_params.ParticleType{n},'d',1)  
                    line([xp - len1/2.*cos(th_p), xp + len1/2.*cos(th_p)]', ...
                        [yp - len1/2.*sin(th_p), yp + len1/2.*sin(th_p)]', ...
                        'color','r','linewidth',0.75,'linestyle','-');
                    line([xp - len2/2.*cos(th_p-pi/2), xp + len2/2.*cos(th_p-pi/2)]', ...
                        [yp - len2/2.*sin(th_p-pi/2), yp + len2/2.*sin(th_p-pi/2)]', ...
                        'color','b','linewidth',0.75,'linestyle','-'); 
                elseif strncmp(run_params.ParticleType{n},'r',1)
                    line([xp - len1/2.*cos(th_p), xp + len1/2.*cos(th_p)]', ...
                        [yp - len1/2.*sin(th_p), yp + len1/2.*sin(th_p)]', ...
                        'color','r','linewidth',0.75,'linestyle','-'); 
                end
            end 
            hold off;
            pause(1/100); 

        end
   
    end

end

% combine keypoints into one array
keypts = cat(3,keypts1,keypts2,keypts3,keypts4,keypts5);



%% MERGE CAMERA VIEWS

% set up figure
if plot_merge_on
    Mfig = figure; set(Mfig,'position',[103  80  1192  568]); 
    subplot(211); axis equal; axis([0 4*img_ix 0 img_iy]); colormap gray
    subplot(212); axis equal; axis([-.5 .5 -.45 .05]); grid on

    imgset = cell(length(cams),1);
    for cam = 1:4
        imgset{cam} = dir(sprintf('%sCam%s*',dir_name,cams(cam))); 
    end
end

% find mean free surface elevation
z_freesurf_m_mean = mean(z_freesurf_m(:,2));

% apply calibration: convert point coords in px into meters across all images
keypts_rect = cell(img_nt,5);
for i = i0:i0+nframes-1  
    for j = 1:5
        % coordinate transformation
        keypts_rect{i,j} = [rectify_quad{1}(keypts{i,1,j}); 
            rectify_quad{2}(keypts{i,2,j}); 
            rectify_quad{3}(keypts{i,3,j}); 
            rectify_quad{4}(keypts{i,4,j})];

        % offset z coordinates so that z=0 is at the mean free surface
        keypts_rect{i,j}(:,2) = keypts_rect{i,j}(:,2) - z_freesurf_m_mean;
    end
end

% transform free surface coordinates
z_freesurf_inst_rect = cell(img_nt,1);
for i = i0:i0+nframes-1  
    % coordinate transformation
    z_freesurf_inst_rect{i} = [rectify_quad{1}(z_freesurf_inst{i,1}); 
        rectify_quad{2}(z_freesurf_inst{i,2}); 
        rectify_quad{3}(z_freesurf_inst{i,3}); 
        rectify_quad{4}(z_freesurf_inst{i,4})];

    % offset z coordinates so that z=0 is at the mean free surface
    z_freesurf_inst_rect{i}(:,2) = z_freesurf_inst_rect{i}(:,2) - z_freesurf_m_mean;
end

% convert rectified keypoints back into centroids and angles
% preallocate variables
centers = cell(img_nt,1);  % particle centroids [xp, yp]
angles = cell(img_nt,1);   % particle orientation info [th_p, d_p]
errchk = cell(img_nt,1);   % error check quantity

% point overlap threshold [m]: if two points are separated by less than
% this threshold, they are determined to belong to the same particle (when
% one particle is imaged by two cameras at the same time)
overlap_thres = 10e-3; 

for i = i0:i0+nframes-1  
    
    % rectified keypoints
    maj1 = keypts_rect{i,1};
    maj2 = keypts_rect{i,2};
    min1 = keypts_rect{i,3};
    min2 = keypts_rect{i,4};
    xpyp = keypts_rect{i,5};

    % disks
    if strncmp(run_params.ParticleType{n},'d',1)
        % centroid
        cen = [ mean([maj1(:,1), maj2(:,1), min1(:,1), min2(:,1)],2), ...
             mean([maj1(:,2), maj2(:,2), min1(:,2), min2(:,2)],2) ];
        
        % orientation
        th_p_rect = atan2( (maj2(:,2) - maj1(:,2))/2, (maj2(:,1) - maj1(:,1))/2 );
        d_p_rect = sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 );
        
        % error check
        ec = sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 );
    
    % rods
    elseif strncmp(run_params.ParticleType{n},'r',1)
        % centroid
        cen = [ mean([maj1(:,1), maj2(:,1)],2), mean([maj1(:,2), maj2(:,2)],2) ];

        % orientation
        th_p_rect = atan2( (maj2(:,2) - maj1(:,2))/2, (maj2(:,1) - maj1(:,1))/2 );
        d_p_rect = sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 );
        
        % error check
        ec = sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 );
        
    % spheres
    else
        % centroid
        cen = xpyp;

        % error check
        ec = mean([ sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 ), ...
            sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 ) ], 2);         
    end

    % accumulate particle data
    centers{i} = cen;
    errchk{i} = ec;
    if nonsphere
        angles{i} = [th_p_rect, d_p_rect];
    end

    % combine particles that appear in more than one camera view:
    % find overlapping particles
    Np = size(cen,1);
    distx = repmat(cen(:,1),1,Np) - repmat(cen(:,1)',Np,1);
    disty = repmat(cen(:,2),1,Np) - repmat(cen(:,2)',Np,1);
    distp = sqrt(distx.^2 + disty.^2);
    [idx1,idx2] = find(distp > 0 & distp < overlap_thres);

    % remove the smaller/shorter overlapping particle
    if strncmp(run_params.ParticleType{n},'r',1)
        L_discrim = d_p_rect; % rod length
    else
        L_discrim = errchk{i}; % disk major axis length/sphere diameter
    end
    idx_remove = idx1( L_discrim(idx1) < L_discrim(idx2) );
    centers{i}(idx_remove,:) = [];
    errchk{i}(idx_remove) = [];
    if nonsphere
        angles{i}(idx_remove,:) = [];
    end

    % plot image array and final particle coordinates
    if plot_merge_on
        maj1(idx_remove,:) = [];
        maj2(idx_remove,:) = [];
        min1(idx_remove,:) = [];
        min2(idx_remove,:) = [];

        A = zeros(img_iy, 4*img_ix);
        for cam = 1:4
            cam_left = cam <= 2;
            A(:,img_ix*(cam-1)+1:img_ix*cam) = cam_imread([dir_name imgset{cam}(i).name], cam_left);
        end

        subplot(211); cla; imagesc(flipud(A),[0 80]); axis equal; axis([0 4*img_ix 0 img_iy]); 
        colormap gray; set(gca,'YTick',[],'XTick',[]); 
        subplot(212); cla;
        if ~isempty(maj1)
            % plot axes 
            line([maj1(:,1), maj2(:,1)]',[maj1(:,2), maj2(:,2)]', ...
                'color','r','linewidth',0.75,'linestyle','-'); hold on
            line([min1(:,1), min2(:,1)]',[min1(:,2), min2(:,2)]', ...
                'color','b','linewidth',0.75,'linestyle','-'); hold off
        end 
        pause(.1); 
    end

end



% save centers and angles
if save_on
    if nonsphere
        save(sprintf('%s\\outputs_%s\\centers_run%02d.mat', code_path, expt_name, run_params.Run(n)),'centers','angles','errchk','z_freesurf_inst_rect');
    else
        save(sprintf('%s\\outputs_%s\\centers_run%02d.mat', code_path, expt_name, run_params.Run(n)),'centers','errchk','z_freesurf_inst_rect');
    end
end

toc  % time the computation



%% TRACK PARTICLES
fs = run_params.imagingFreq_Hz(n);  % imaging freq (Hz)
searchrad = 10e-3; % search radius (m)

tic  % time the computation

% use nearest-neighbor algorithm to track particles
if nonsphere
    [tracks0,tracklength0] = get_particle_tracks(centers,fs,searchrad,angles,errchk);
else
    [tracks0,tracklength0] = get_particle_tracks(centers,fs,searchrad,errchk);
end
toc  % time the computation

% save tracks
if save_on
    save(sprintf('%s\\outputs_%s\\tracks_run%02d.mat', code_path, expt_name, run_params.Run(n)),'tracks0','tracklength0')
end

% preprocess tracks
avar_k = preprocess_MP(n);
