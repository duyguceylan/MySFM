The libraries "sfm-driver" and "sba-1.5" include our code for the symmetry-aware structure-from-motion (SfM) optimization (Section 4.3 of our paper). The implementation is done on top of the original versions of these libraries by adding additional functions. (The remaining code has not been modified.) "sfm-driver" library is linked with the "sba-1.5" library. The original "sfm-driver" library comes with the open source distribution of Bundler written by Noah Snavely (http://www.cs.cornell.edu/~snavely/bundler/). Since this library also includes the "matrix" and "imagelib" libraries that also come with this code, we also include them. We have not made any changes in these libraries.

In order to call the symmetry-aware SfM optimization routine from your work, you should link your code with the "sfm-driver" library. You can call the "run_sfm_grid_opt" function. This function assumes, grid and non-grid correspondences have been established. Initial estimates of the 
	1. 3D position of the non-grid correspondence points, 
	2. the 3D position of the reference grid point and the vertical and horizontal transformations for the detected grids, 
	3. camera parameters
should be provided. An explanation of the arguments the function takes is provided in "sfm.c".

You can use the makefile files included together with the libraries to compile the code. You should first compile the "matrix", "imagelib", and "sba" libraries and modify the makefile for the "sfm-driver" library with the correct path to these libraries. You can then compile the "sfm-driver" library. 

We would like to thank Noah Snavely and Manolis Lourakis for making their code online.

We have tested the code on Mac Os X environment and we expect that there should not be a problem with compilation on Linux environments as well. We haven't tested on Windows platform.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

If you use this code, please cite our paper:
Coupled Structure-from-motion and 3D Symmetry Detection for Urban Facades
Ceylan, Duygu and Mitra, Niloy J. and Zheng, Youyi and Pauly, Mark
ACM Trans. Graph., January 2014