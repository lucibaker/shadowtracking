% plot images

% path containing images
dir_name = 'D:\luci\run220613\run4\';

% image range to show
% rng = 1:100;
rng = 490:530; %smtracks(smtracks(:,5)==45 | smtracks(:,5)==78, 7);
rng = min(rng):max(rng);

% read image set
cams = 'ACBD';
imgset = cell(size(cams));
for cam = 1:length(cams)
    imgset{cam} = dir(sprintf('%sCam%s*', img_path,cams(cam)));
end
B = imread(sprintf('%s%s',dir_name, imgset{1}(1).name));
[img_ix, img_iy] = size(B);

% figure
fig = figure; set(gcf,'position',[0.0010    0.0410    2.5600    1.3273]*1000);

% show images
A = zeros(img_iy, 4*img_ix);
for frameno = rng    
    for cam = 1:length(cams)
        cam_left = cam <= 2;
        A(:,img_ix*(cam-1)+1:img_ix*cam) = cam_imread([dir_name imgset{cam}(frameno).name], cam_left);
    end
    cla; imagesc(flipud(A),[0 80]); axis equal; axis([0 4*img_ix 0 img_iy]); colormap gray; title(num2str(frameno))
    pause(1/10)
end