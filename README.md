# shadowtracking: Shadow tracking code for WASIRF plastic particle experiments

## References
Baker, L. J., & Coletti, F. (2022). Experimental investigation of inertial fibres and disks in a turbulent boundary layer. Journal of Fluid Mechanics, 943, A27. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/experimental-investigation-of-inertial-fibres-and-disks-in-a-turbulent-boundary-layer/0261D874F98C7C0981D1B06ACB8223DD

Baker, L., & DiBenedetto, M. (2023). Large-scale particle shadow tracking and orientation measurement with collimated light. Experiments in Fluids, 64(3), 52. https://link.springer.com/article/10.1007/s00348-023-03578-y

## Inputs
### run_parameters_demo.xlsx
Spreadsheet containing the experimental parameters for each run. Use the demo spreadsheet as a template.

### cal_parameters_demo.xlsx
Spreadsheet containing the calibration parameters for each run. Use the demo spreadsheet as a template.

### calibration_and_bkgd
Folder containing calibration image, background image, and mask image for each camera.

### run1, run2, run3
Folders containing the images for each run

## Instructions


## Function descriptions
### detect_MP
main plastic particle detection and tracking script

### preprocess_MP
preprocess particle tracks:
 -clean up data
 -compute orientation (p vector)
 -apply temporal smoothing to tracks and orientations

### calibrate_camera
calculate the transformation function to convert image coordinates to world (physical) coordinates, including calibration (converting pixels to meters), undistortion, and rectification

### calibrate_camera_DEMO
demonstrates how to use calibrate_camera as a standalone function

### cam_imread
function to correctly load and rotate tif or bmp images from 4x1 camera array

### pcolor_img
helper function to display camera images

### track functions
functions to operate on particle track data

## Troubleshooting



