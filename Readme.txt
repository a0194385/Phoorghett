This is a ultrasound imaging simulator base on MATLAB GUI.
Run UltrasoundImagingSimulator_Zheheng.m in MATLAB to open the simulator interface.
It simulates the L1-regularization of the linear inversion problem for the computational ultrasound imaging.
With Born approximation, ultrasound imaging is approximated as a linear inversion problem.

In this simulator, Field II (version 3.24a windows) is used to generate ultrasound data.
http://www.field-ii.dk//

Top-left area of the simulator interface is to set the parameters.
Below it, "Display" buttom is to show the ultrasound beam in imaging domain as well as ultrasound pulse.
Below, "Display phantom" button is to show the linear-scale phantom, and "Generate RF data" is to generate RF data of phantom.
Bottom-left area  is to set the frequency band of the signals.
"Generate system matrix" button at top-right area is to generate system matrix A, which contains the echo data of one point scatterer at the center of each pixel.
Below it, by clicking "LSQR reconstruction", LSQR algorithm solves the linear equations y=Ax to reconstruct the image and reconstructed image is shown with log-compression.

LSQR_B (Per Christian Hansen, IMM, August 13, 2001.) is used.
https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/52/versions/2/previews/lsqr_b.m/index.html

Bottom-right area is to solve the linear equations with L1-regularization by YALL1 solver. 
Selecting multiple bases means using the conjugation of these bases for L1-regularization.

Here YALL1 v1.4 is used.
http://yall1.blogs.rice.edu/

The files do not include YALL1 and Field II, which have to be downloaded by users.

After changing the parameters, new system matrix and RF data should be generated.
Because the image reconstruction will work correctly only with proper system matrix and RF data.

Users can use the sparse transform defined by users.
To use a transform matrix, users can save the matrix as variable 'phi' in a .mat file, and then write the "filename.mat" in the edit box next to the "user-defined" button.
To use a transform function, users can write the transform and inverse transform handle in "transform_handles.m".

Zheheng Liu. Sept. 16, 2020.
