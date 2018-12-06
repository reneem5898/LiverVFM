function liverContour = getLiverContour(directory, seq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns coordinates (in image) of the liver contour
% (made by LiverLab)
%
% Input: 1) patient directory
%        2) seq - name of sequence
%
% Output: contour coordinates of liver
%
% Written by: Renee Miller (reneem5898@gmail.com)
% Date: 4 December 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list of all directories
d = sprintf('%s%sP_StiffC95', directory, seq);

% allDirs = dir(directory);
% 
% % Get directory with stiffness map
% for i = 1:length(allDirs)
%     if strfind(allDirs(i).name, 'StiffC95')
%         stiffDir = allDirs(i).name;
%     end
% end

% Open an image (any image) in the stiffness map directory
%AllFiles = dir(sprintf('%s/%s', directory, stiffDir));
AllFiles = dir(d);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Open an image and get the contour from the dicom header
info = dicominfo(sprintf('%s/%s', d, imageFiles{1}));
contourImg = info.OverlayData_0;

% Find indices of non-zero elements in contour map and save as liverContour
[x, y] = find(contourImg);
liverContour = [x y];