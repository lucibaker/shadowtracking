% DEMO for calibrate_camera.m

addpath('..\')

% calibrate the cameras
load('calib_data_C.mat')  % contains coordinates of calibration points in image coordinates (I matrix) and world coordinates (W matrix) 
    % (note: this is the calibration data from the WASIRF experiment)
rectify = calibrate_camera(I,W,2);  % get image-to-world transformation function

% get particle positions in world coordinates
points_px = (rand([1,2]) - 0.5)*2000;  % randomly generated particle position in image coordinates [units: pixels] for testing our rectify function
points_m = rectify(points_px);  % convert particle position from image coordinates to world coordinates [units: meters]
