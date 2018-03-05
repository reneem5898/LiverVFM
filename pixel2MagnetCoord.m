function magnetCoord = pixel2MagnetCoord(img, coord)

% ------------------------------------------------------------------------
% This function converts one pixel from an image into magnet coordinates
% using the header information
% 
% Inputs: 1) image
%         2) coordinates in pixel coordinates
% Output  1) set of magnet coordinates for that pixel
%
% Written by: Renee Miller
% Last modified: 3 November 2015
%
% Function called from: main_dti.m, getCardiacPoints.m
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

% Initialise transformation matrix (mag coordinate system to
% cardiac)
pix2coord = zeros(4,4);

% Construct the matrix - see documentation
% CoordinateSystemConversion.docx
pix2coord(1:3,1) = colVec*dc;
pix2coord(1:3,2) = rowVec*dr;
pix2coord(1:3,4) = TL;
pix2coord(4,4) = 1.0;

r = coord(2); % row value is the y coord
c = coord(1); % column value is the x coord

coordinate = [r-1; c-1; 0; 1];
magnetCoord = pix2coord*coordinate;
magnetCoord = magnetCoord(1:3);