function A = cam_imread(fn,cam_left)
% function to load and rotate images from 4x1 camera array
%
% A = loaded and rotated image
% fn = image filename
% cam_left: 'true' for the left two cameras, 'false' for right two cameras

try
    if cam_left
        A = imread(fn)'; 
    else
        A = rot90(imread(fn),2)'; 
    end
catch
    A = rgb2gray(imread(fn));
    if cam_left
        A = A'; 
    else
        A = rot90(A,2)'; 
    end
end

end