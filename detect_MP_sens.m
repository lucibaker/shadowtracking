function detect_MP_sens(n)
% input n = run number
% MP particle detection and coordinate mapping

%% SETUP

% clear
% close all

% n = 1;  % run number
gdrive_path = 'C:\Users\ljbak\My Drive\';  % 'G:\My Drive\';  %  
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\sens_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
cal_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\cal_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

% make plots or not
plot_on1 = 0; % calibration
plot_on2 = 1; % detection
plot_on3 = 0; % merged views

% save results or not
save_on = 1;

% image parameters 
cam0 = 1;
cams = cell2mat(cal_params.Cam)';
dir_name = sprintf('sens%i\\', run_params.Run(n));
imgset = dir(sprintf('%sCam%s_sens*.tif',dir_name,cams(cam0))); 
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
if plot_on1
    Bfig = figure; 
    set(Bfig,'Units', 'Pixels','Position',[0.2642    0.1274    1.1368    0.6200]*1000); 
    set(gcf,'color','w');
end



%% LOOP OVER CAMERAS

% tic

% get array of coords of relevant points in px for each cam: endpoints of maj/min axes for disks,
% endpoints of minor axis for rods, centroid for nurdles/wax pellets
% centers1 = cell(img_nt,length(cams));  % particle centroids [xp, yp] (distorted)
% angles1 = cell(img_nt,length(cams));   % particle orientation info [th_p, d_p] (distorted)

keypts1 = cell(img_nt,length(cams));  % particle key points (endpoints of major and minor axes and centroid 
                                       % {xp_maj1 yp_maj1}, {xp_maj2 yp_maj2}, {xp_min1 yp_min1}, {xp_min2 yp_min2}, {xp, yp})
keypts2 = cell(img_nt,length(cams));
keypts3 = cell(img_nt,length(cams));
keypts4 = cell(img_nt,length(cams));
keypts5 = cell(img_nt,length(cams));

z_freesurf_inst = cell(img_nt,length(cams));

rectify_quad = cell(length(cams),1);      % image rectification function handles for each camera view

for cam = 1:2%1:length(cams)   
    
    cam_left = cam <= 2; % 'true' for the left two cameras, 'false' for right two cameras
    
    % read image files from directory structure 
    imgset = dir(sprintf('%sCam%s_sens*.tif',dir_name,cams(cam))); 
    imgset_mask = dir(sprintf('%sCam%s_mask*.tif',dir_name,cams(cam)));

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
    if plot_on1
        figure; subplot(161); pcolor_img(cal); title('bkgd sub and inverted')
        subplot(162); pcolor_img(B); title('masked and binarized')
    end

    % detect dots
    CC = bwconncomp(B);
    cal_dots = regionprops('table',CC,cal,'Centroid','Area','EquivDiameter');
    idx = cal_dots.Area > cal_params.calAreaThres1_px(cam) & cal_dots.Area < cal_params.calAreaThres2_px(cam);
%     idx = cal_dots.Area > 300 & cal_dots.Area < 5000;
    cal_dots = cal_dots(idx,:);    
    if plot_on1
        subplot(163); pcolor_img(cal); hold on
        viscircles(cal_dots.Centroid,cal_dots.EquivDiameter/2); title('dots detected')
    end

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
        if plot_on1; line([bin_lim(j) bin_lim(j)],[0 img_iy]); end  % plot bin limits
        cal_dots0 = cal_dots(cal_dots.Centroid(:,1) > bin_lim(j) & cal_dots.Centroid(:,1) < bin_lim(j+1),:);  % dots in a vertical bin
        [~,sort_idx] = sort(cal_dots0.Centroid(:,2),'ascend');  % sort by y-coord
        I(n_rows_cal*(j-1)+1 : n_rows_cal*j,:) = cal_dots0.Centroid(sort_idx,:);  % image points
    end
    if plot_on1  % plot detected image points
        point_colors = jet(size(I,1));
        subplot(164); scatter(I(:,1),I(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('detected points')
        subplot(165); scatter(W(:,1),W(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('known coords')
    end

    % undistortion and rectification function
    rectify_quad{cam} = calibrate_camera(I,W,2);    
    
    % confirm that calibration dots are mapped to the correct world coords
    if plot_on1
        imagePoints2 = rectify_quad{cam}(I); 
        subplot(166); scatter(imagePoints2(:,1),imagePoints2(:,2),15,point_colors,"filled"); grid on; axis tight equal; title('corrected points')
    end

    % quiescent free surface coordinates
    z_freesurf_m(2*(cam-1)+1:2*cam,:) = rectify_quad{cam}([cal_params.xSurfStill1_px(cam), cal_params.ySurfStill1_px(cam); ...
        cal_params.xSurfStill2_px(cam), cal_params.ySurfStill2_px(cam)]);

    if plot_on1; figure(Bfig); end


    %% LOOP OVER FRAMES
    i0 = 1;
    nframes = img_nt;  
    
    for i = i0:i0+nframes-1 % parfor

        A = cam_imread([dir_name imgset(i).name], cam_left);
        M = cam_imread([dir_name imgset_mask(i).name], cam_left);

        A0 = imcrop(A,img_crop_rect);  % crop to correct for camera shift
        M = imcrop(M,img_crop_rect);
        A0 = double(A0) - bkgd; % subtract background
        A0 = uint8(A0 - min(A0(:)));  % shift intensities so that all are positive 
        A0 = 255 - A0;  % invert image 
    
        % adaptive binarization
        B = imbinarize(A0,'adaptive','Sensitivity',ABsens); % Adaptative binarization   
        B = logical(B.*logical(M));
        
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
%             freesurf_col_j = S_freesurf.PixelList{1}(:,1)==j;
%             if any(freesurf_col_j)
%                 z_freesurf_px(j) = min(S_freesurf.PixelList{1}(freesurf_col_j,2));
%             else
                z_freesurf_px(j) = inf;
%             end
        end
        z_freesurf_inst{i,cam} = [(1:length(z_freesurf_px))', z_freesurf_px];
        
        xp = []; yp = [];
        if ~isempty(S1)
            % remove based on area and proximity to free surface
            idx = S1.Area > A_thres(1) & S1.Area < A_thres(2) & S1.Centroid(:,2) < z_freesurf_px(round(S1.Centroid(:,1))) & ...
                S1.MajorAxisLength > maj_thres(1) & S1.MajorAxisLength < maj_thres(2) & S1.MinorAxisLength > min_thres(1) & S1.MinorAxisLength < min_thres(2); 
            S = S1(idx,:);

            if ~isempty(S)
%                 % remove doubles (false shadows when the particle that shows up in addition to its shadow 
%                 % if particle is near the back wall) 
%                 dx_double = 300; % max horizontal distance between doubles
%                 dx_double_l = dx_double*~cam_left; % if particle is to the left of the shadow/right camera (px)
%                 dx_double_r = -dx_double*cam_left; % if particle is to the right of the shadow/left camera (px)
%                 dz_double = 20; % max vertical distance between doubles (px)
%                 
%                 Np = height(S);
%                 xp = S.Centroid(:,1);
%                 yp = S.Centroid(:,2);
%                 distx = repmat(xp,1,Np) - repmat(xp',Np,1);
%                 disty = repmat(yp,1,Np) - repmat(yp',Np,1);
%                 [~,double_idx] = find(distx ~= 0 & disty ~= 0 & distx < dx_double_l & distx > dx_double_r & disty > -dz_double & disty < dz_double);
%                 if ~isempty(double_idx)
%                     S(double_idx,:) = [];
%                 end
                
                % particle centroids and orientations
                Np = height(S);
                xp = S.Centroid(:,1);
                yp = S.Centroid(:,2); 

                if strncmp(run_params.ParticleType{n},'r',1)        % rods
                    len1 = sqrt(S.BoundingBox(:,3).^2 + S.BoundingBox(:,4).^2);  % length
                    len2 = S.Area./len1; % thickness
                    len1 = len1 - 0.25*len2;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
                    d_p = len1; % apparent major axis length [px]
                else                                                % disks and spheres
                    len1 = S.MajorAxisLength;
                    len2 = S.MinorAxisLength;
                    th_p = -S.Orientation(:)*pi/180; % in-plane orientation, radians (need (-S.Orientation) because image is flipped up-down)
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
%             centers1{i,cam} = [xp,yp];
%             if nonsphere
%                 angles1{i,cam} = [th_p,d_p];
%             end
            keypts1{i,cam} = [xp - len1/2.*cos(th_p), yp - len1/2.*sin(th_p)];
            keypts2{i,cam} = [xp + len1/2.*cos(th_p), yp + len1/2.*sin(th_p)];
            keypts3{i,cam} = [xp - len2/2.*cos(th_p-pi/2), yp - len2/2.*sin(th_p-pi/2)]; 
            keypts4{i,cam} = [xp + len2/2.*cos(th_p-pi/2), yp + len2/2.*sin(th_p-pi/2)];
            keypts5{i,cam} = [xp, yp];
        else
%             centers1{i,cam} = zeros(0,2);
%             if nonsphere
%                 angles1{i,cam} = zeros(0,2);
%             end
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
    
        % MAKE PLOTS
        if plot_on2
            subplot(121); cla; pcolor_img(B); 
            subplot(122); cla; pcolor_img(A0); hold on
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
                        'color','y','linewidth',0.75,'linestyle','-'); 
                elseif strncmp(run_params.ParticleType{n},'r',1)
                    line([xp - len1/2.*cos(th_p), xp + len1/2.*cos(th_p)]', ...
                        [yp - len1/2.*sin(th_p), yp + len1/2.*sin(th_p)]', ...
                        'color','r','linewidth',0.75,'linestyle','-'); 
                end
            end 
            hold off;
            pause(1/100); 
            
%             % write to avi
%             if i == i0
%                 mname = 'test_disks_bkgd'; % movie file name 
%                 vid = VideoWriter(mname,'Grayscale AVI');
%                 vid.FrameRate = 10;
%                 open(vid)
%             end
%             writeVideo(vid,flipud(uint8(A0)));
%             if i == i0+nframes
%                 close(vid);
%             end
        
%             % write to gif
%             fn_gif = sprintf('C:\\Users\\ljbak\\My Drive\\MP in OSBL\\imaging expts\\%s-%i-%s.gif',run_params.ParticleType{n},run_params.WindSpeed_m_s(n),cams(cam));
%             fig_to_gif(fn_gif,0.1)
        end
   
    end

end
keypts = cat(3,keypts1,keypts2,keypts3,keypts4,keypts5);



%% MERGE CAMERA VIEWS

if plot_on3
%     Mfig1 = figure; set(Mfig1,'position',[0.0010    0.0410    1.5360    0.7488]*1000);
%     axis equal; axis([0 4*img_ix 0 img_iy]); colormap gray
%     Mfig2 = figure; set(Mfig2,'position',[-1.9190   -0.1750    1.9200    0.9648]*1000);
%     axis equal; axis([-.5 .5 -.45 .05]); grid on
    Mfig = figure; set(Mfig,'position',[1.1983    0.1397    1.3400    1.1000]*1000);
    subplot(211); axis equal; axis([0 4*img_ix 0 img_iy]); colormap gray
    subplot(212); axis equal; axis([-.5 .5 -.45 .05]); grid on

    imgset = cell(length(cams),1);
    for cam = 1:2%1:4
        imgset{cam} = dir(sprintf('%sCam%s*',dir_name,cams(cam))); 
    end

    if strcmp(fmt,'avi')
        vid = cell(length(cams),length(imgset));
        for cam = 1:2%1:4
            for j = 1:max(file_id)
                vid{cam,j} = VideoReader([dir_name imgset{cam}(j).name]);
            end
        end
    end
end

z_freesurf_m_mean = mean(z_freesurf_m(:,2));

% apply calibration: convert point coords in px into meters across all images
keypts_rect = cell(img_nt,5);
for i = i0:i0+nframes-1  
    for j = 1:5
        keypts_rect{i,j} = [rectify_quad{1}(keypts{i,1,j}); rectify_quad{2}(keypts{i,2,j}); ...
            rectify_quad{3}(keypts{i,3,j}); rectify_quad{4}(keypts{i,4,j})];
        keypts_rect{i,j}(:,2) = keypts_rect{i,j}(:,2) - z_freesurf_m_mean;
    end
end

% convert freesurf coords
z_freesurf_inst_rect = cell(img_nt,1);
for i = i0:i0+nframes-1  
    z_freesurf_inst_rect{i} = [rectify_quad{1}(z_freesurf_inst{i,1}); rectify_quad{2}(z_freesurf_inst{i,2}); ...
        rectify_quad{3}(z_freesurf_inst{i,3}); rectify_quad{4}(z_freesurf_inst{i,4})];
    z_freesurf_inst_rect{i}(:,2) = z_freesurf_inst_rect{i}(:,2) - z_freesurf_m_mean;
end

% convert rectified keypoints back into centroids and angles
centers = cell(img_nt,1);  % particle centroids [xp, yp]
angles = cell(img_nt,1);   % particle orientation info [th_p, d_p]
errchk = cell(img_nt,1);   % error check quantity

overlap_thres = (2e-3)^2; % m (squared)

for i = i0:i0+nframes-1  
    
    % rectified keypoints
    maj1 = keypts_rect{i,1};
    maj2 = keypts_rect{i,2};
    min1 = keypts_rect{i,3};
    min2 = keypts_rect{i,4};
    xpyp = keypts_rect{i,5};

    if strncmp(run_params.ParticleType{n},'d',1)
        cen = [ mean([maj1(:,1), maj2(:,1), min1(:,1), min2(:,1)],2), ...
             mean([maj1(:,2), maj2(:,2), min1(:,2), min2(:,2)],2) ];
        th_p_rect = atan2( (maj2(:,2) - maj1(:,2))/2, (maj2(:,1) - maj1(:,1))/2 );
        d_p_rect = sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 );
        ec = sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 );
    
    elseif strncmp(run_params.ParticleType{n},'r',1)
        cen = [ mean([maj1(:,1), maj2(:,1)],2), mean([maj1(:,2), maj2(:,2)],2) ];
        th_p_rect = atan2( (maj2(:,2) - maj1(:,2))/2, (maj2(:,1) - maj1(:,1))/2 );
        d_p_rect = sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 );
        ec = sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 );
        
    else
        cen = xpyp;
        ec = mean([ sqrt( (maj2(:,1) - maj1(:,1)).^2 + (maj2(:,2) - maj1(:,2)).^2 ), ...
            sqrt( (min2(:,1) - min1(:,1)).^2 + (min2(:,2) - min1(:,2)).^2 ) ], 2);         
    end

    % accumulate particle data
    centers{i} = cen;
    errchk{i} = ec;
    if nonsphere
        angles{i} = [th_p_rect, d_p_rect];
    end

    % combine particles that appear in more than one camera view
    Np = size(cen,1);
    distx = repmat(cen(:,1),1,Np) - repmat(cen(:,1)',Np,1);
    disty = repmat(cen(:,2),1,Np) - repmat(cen(:,2)',Np,1);
    distp = distx.^2 + disty.^2;
    [idx1,idx2] = find(distp > 0 & distp < overlap_thres);

    % remove the smaller/shorter particle
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

    % plot image array and rectified particles
    if plot_on3
%         maj1(idx_remove,:) = [];
%         maj2(idx_remove,:) = [];
%         min1(idx_remove,:) = [];
%         min2(idx_remove,:) = [];

        A = zeros(img_iy, 4*img_ix);
        for cam = 1:2%1:4
            cam_left = cam <= 2;
            if strcmp('fmt','tif')
                A(:,img_ix*(cam-1)+1:img_ix*cam) = cam_imread([dir_name imgset{cam}(i).name], cam_left);
            else
                A(:,img_ix*(cam-1)+1:img_ix*cam) = cam_aviread(vid{cam,file_id(i)}, i-start_idx(file_id(i))+1, cam_left);
            end
        end

%         figure(Mfig1); 
        subplot(211); cla; imagesc(flipud(A),[0 80]); axis equal; axis([0 4*img_ix 0 img_iy]); colormap gray
%         figure(Mfig2); 
        subplot(212); cla;
        if ~isempty(maj1)
            % plot axes 
            line([maj1(:,1), maj2(:,1)]',[maj1(:,2), maj2(:,2)]', ...
                'color','r','linewidth',0.75,'linestyle','-'); hold on
            line([min1(:,1), min2(:,1)]',[min1(:,2), min2(:,2)]', ...
                'color','b','linewidth',0.75,'linestyle','-'); hold off
        end 
        pause(.1); 
%         fig_to_gif('C:\Users\ljbak\My Drive\MP in OSBL\imaging expts\figs\d5-16-combined.gif',0.1)
    end

end



%% SAVE CENTERS AND ANGLES
if save_on
    if nonsphere
        save(sprintf('%sMP in OSBL\\imaging expts\\run%s\\outputs\\centers_sens%02d.mat', gdrive_path, expt_string, run_params.Run(n)),'centers','angles','errchk','z_freesurf_inst_rect');
    else
        save(sprintf('%sMP in OSBL\\imaging expts\\run%s\\outputs\\centers_sens%02d.mat', gdrive_path, expt_string, run_params.Run(n)),'centers','errchk','z_freesurf_inst_rect');
    end
end

% toc


%% TRACK PARTICLES
fs = run_params.imagingFreq_Hz(n);  % imaging freq (Hz)
searchrad = 10e-3; % search radius (m)

tic
if nonsphere
    [tracks0,tracklength0] = get_particle_tracks(centers,fs,searchrad,angles,errchk);
else
    [tracks0,tracklength0] = get_particle_tracks(centers,fs,searchrad,errchk);
end
toc

if save_on
    save(sprintf('%sMP in OSBL\\imaging expts\\run%s\\outputs\\tracks_sens%02d.mat', gdrive_path, expt_string, run_params.Run(n)),'tracks0','tracklength0')
end
