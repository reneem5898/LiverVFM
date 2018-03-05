function liverContour = getLiverContour(directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns coordinates (in image) of the liver contour
% (made by LiverLab)
%
% Input: patient directory
%
% Output: contour coordinates of liver
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date: 12 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list of all directories
allDirs = dir(directory);

% Get directory with stiffness map
for i = 1:length(allDirs)
    if strfind(allDirs(i).name, 'StiffC95')
        stiffDir = allDirs(i).name;
    end
end

% Open an image (any image) in the stiffness map directory
AllFiles = dir(stiffDir);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Open an image and get the contour from the dicom header
info = dicominfo(imageFiles{1});
contourImg = info.OverlayData_0;

% Find indices of non-zero elements in contour map and save as liverContour
[x, y] = find(contourImg);
liverContour = [x y];