# shadowtracking: Shadow tracking code for WASIRF plastic particle experiments

## References

Baker, L. J., & Coletti, F. (2022). [Experimental investigation of inertial fibres and disks in a turbulent boundary layer.](https://doi.org/10.1017/jfm.2022.438) Journal of Fluid Mechanics, 943, A27.

Baker, L., & DiBenedetto, M. (2023). [Large-scale particle shadow tracking and orientation measurement with collimated light.](https://doi.org/10.1007/s00348-023-03578-y) Experiments in Fluids, 64(3), 52.


## Instructions

1. Create the input files (see *Inputs*)

2. Create a small subset of images (~100 frames). 

3. Run main script (detect_MP) on the subset and tune the thresholds in the run_parameters spreadsheet to optimize particle detection and tracking.

4. Run detect_MP with the full image set.


## Inputs

**run_parameters_demo.xlsx**  
Spreadsheet containing the experimental parameters for each run. Use the demo spreadsheet as a template.

**cal_parameters_demo.xlsx**  
Spreadsheet containing the calibration parameters for each run. Use the demo spreadsheet as a template.

**calibration_and_bkgd**  
Folder containing calibration image, background image, and mask image for each camera.

**run1, run2, run3**  
Folders containing the images for each run. Image filenames must start with the camera name, e.g. 'CamA_0001.tif'


## Outputs

**centers_run01.mat**

Intermediate output. Contains cell arrays of particle centroids and orientations, organized by image frame (one cell per frame)  

**tracks_run01.mat**

Intermediate output. Contains unsmoothed particle tracks and apparent (projected) orientations

**smtracks_run01.mat**

Final output. Contains smoothed particle tracks and 3D orientations. 

The columns of the smtracks array are:

1. x [m]  
2. z [m]  
3. u [m/s]  
4. w [m/s]  
5. track ID  
6. time counter for each track  
7. frame number  
8. a_x [m/s^2]  
9. a_z [m/s^2]  

The columns of the smangles array are:

1. px  
2. py  
3. pz  
4. px_dot [s^-1]  
5. py_dot [s^-1]  
6. pz_dot [s^-1]  
7. px_dd [s^-2]  
8. py_dd [s^-2]  
9. pz_dd [s^-2]

The smangles_cont array is the same as smangles, but orientations don't wrap around from -1 to 1 as the particle rotates, so these values are continuous in time. They can be used for differentiation or autocorrelations. 

The kernel value is also saved in this file as a record of what smoothing kernel width was used. 



## Function descriptions

**detect_MP**  
main script for plastic particle detection and tracking 

**preprocess_MP**  
preprocess particle tracks: clean up data, compute orientation (p vector), apply temporal smoothing to tracks and orientations

**calibrate_camera**  
calculate the transformation function to convert image coordinates to world (physical) coordinates, including calibration (converting pixels to meters), undistortion, and rectification

**cam_imread**  
function to correctly load and rotate tif or bmp images from 4x1 camera array

**pcolor_img**  
helper function to display camera images

**track functions**  
functions to operate on particle track data

**calibration example**  
demonstrates how to use calibrate_camera as a standalone function


## Troubleshooting

**Code takes a long time to run**
  
Switch main for-loop (loop over frames) in detect_MP script to parfor-loop and run on parallel cluster. (reference: [[https://www.mathworks.com/help/parallel-computing/parfor.html]])

**Calibration/rectification failed**

Adjust adaptive binarization thresholds for dot detection in cal_parameters spreadsheet.  
Adjust area and axis length thresholds for dot detection in cal_parameters spreadsheet. 
Check that the order of the dots in the world coordinates array corresponds exactly to order of the dots in the image coordinates vector for each camera view.   

**No particles detected** 
 
Adjust adaptive binarization thresholds for particle detection in run_parameters spreadsheet.  
Adjust area and axis length thresholds for particle detection in run_parameters spreadsheet.  

**Incorrect merging of the four camera views**

Check world offsets of each camera view in cal_parameters spreadsheet.

**No tracks found**

Adjust search radius in run_parameters spreadsheet.

