function liverMask = getLiverMask(directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a mask of the liver from the contour (made by LiverLab)
%
% Input: patient directory
%
% Output: mask of the liver
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
contour = info.OverlayData_0;

% Use region growing algorithm to convert contour to a mask
liverMask = regiongrowing(contour);