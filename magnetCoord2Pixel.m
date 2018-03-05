function pixCoord = magnetCoord2Pixel(img, coord_3d)

% ------------------------------------------------------------------------
% This function converts one cooridnate from the magnet coordinate system
% to pixel coordinates using the header information
% 
% Inputs: 1) image - filename
%         2) 3D coordinates in magnet coordinates
% Output  1) pixel coordinates for that point
%
% Written by: Renee Miller
% Last modified: 20 December 2017
%
% ------------------------------------------------------------------------

header = dicominfo(img); %Read dicom header

% Get the vectors corresponding to rows and columns of image -
% ImageOrientationPatient
orient = header.ImageOrientationPatient;

% Get resolution of images (mm/pixel)
res = header.PixelSpacing;
dc = res(1);
dr = res(2);

% Get vector corresponding to first row and column of image in
% cardiac coordinate system
rowVec = orient(1:3);
colVec = orient(4:6);

% Get the coordinates of the top left pixel - ImagePositionPatient
TL = header.ImagePositionPatient; 

% Initialise transformation matrix 
coord2pix = zeros(4,4);

% Construct the matrix - see documentation
% CoordinateSystemConversion.docx
coord2pix(1:3,1) = colVec*dc;
coord2pix(1:3,2) = rowVec*dr;
coord2pix(1:3,4) = TL;
coord2pix(4,4) = 1.0;

% RHS
coordinate_3d = [coord_3d(1); coord_3d(2); coord_3d(3); 1];

% Convert coordinate to pixel coordinates
pixCoord = coordinate_3d\coord2pix;
pixCoord = pixCoord(1:2);