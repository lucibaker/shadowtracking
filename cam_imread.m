function A = cam_imread(fn,cam_left)
% function to correctly load and rotate tif or bmp images from 4x1 camera array
%
% A = loaded and rotated image
% fn = image filename, e.g. 'CamA-0001.tif'
% cam_left: 'true' for the left two cameras (A and C), 'false' for right
% two cameras (B and D)

try    % try loading image assuming it is black and white
    if cam_left
        A = imread(fn)'; 
    else
        A = rot90(imread(fn),2)';  % images from cam B and D need to be rotated
    end

catch   % if that gives an error, load it as a color image
    A = rgb2gray(imread(fn));
    if cam_left
        A = A'; 
    else
        A = rot90(A,2)';  % images from cam B and D need to be rotated
    end
end

end