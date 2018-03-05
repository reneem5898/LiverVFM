# LiverVFM

These scripts allow one to use the isotropic optimised virtual fields method (as outlined in Connesson 2015 in Strain) to estimate a single complex shear modulus in a region of interest. The main scripts to run are: main_liverVFM.m or main_liverVFM_multiSlice.m, depending on whether stiffness estimates are required for a single or multiple slices. 

Currently, this code is slow because a 3D mesh was used to represent the image. 

There are numerous element types available - varying the element integration method - C3D8 (selectively reduced integration elements), C3D8R (fully reduced integration) and C3D8F (fully integrated - giving the user the choice of the number of Gauss points).
