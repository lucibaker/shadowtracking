function A = cam_aviread(vid,N,cam_left)
% function to load and rotate images from 4x1 camera array
%
% A = loaded and rotated image
% vid = avi VideoReader object
% N = frame number
% cam_left: 'true' for the left two cameras, 'false' for right two cameras

try
    if cam_left
        A = read(vid,N)'; 
    else
        A = rot90(read(vid,N),2)'; 
    end
catch
    A = rgb2gray(read(vid,N));
    if cam_left
        A = A'; 
    else
        A = rot90(A,2)'; 
    end
end

end