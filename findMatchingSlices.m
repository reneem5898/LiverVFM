function [fileListOut1, fileListOut2] = findMatchingSlices(dir1, dir2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findMatchingSlices - returns file list with common slice locations 
% between two folders
%
% Inputs: 1) dir1 - first directory
%         2) dir2 - second directory 
%
% Outputs: 1) fileListOut1 - list of files to use from first directory
%          2) fileListOut2 - list of files to use from second directory
%
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date modified: 19 February 2018
%
% Called from: getTransform.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List all files in two directories

% Get list of all files in first directory
AllFiles1 = dir(dir1);
fileList1 = {AllFiles1.name};
fileList1(ismember(fileList1,{'.','..'})) = [];

% Get list of all files in second directory
AllFiles2 = dir(dir2);
fileList2 = {AllFiles2.name};
fileList2(ismember(fileList2,{'.','..'})) = [];


%% Initialise variables
sliceLocList1 = zeros(length(fileList1),1);
sliceLocList2 = zeros(length(fileList2),1);
fileListOut1 = fileList1;
fileListOut2 = fileList2;


%% Get slice locations

% Get slice locations - list 1
for i = 1:length(fileList1)
    
    % Get current image files
    fn = sprintf('%s/%s', dir1, fileList1{i});
       
    % Save image info temporarily
    info = dicominfo(fn);
        
    % Get slice height
    sliceLocList1(i) = info.SliceLocation;

end

% Get slice locations - list 2
for i = 1:length(fileList2)
    
    % Get current image files
    fn = sprintf('%s/%s', dir2, fileList2{i});
       
    % Save image info temporarily
    info = dicominfo(fn);
        
    % Get slice height
    sliceLocList2(i) = info.SliceLocation;

end

%% Get intersection of slice location lists
slicesCommon = intersect(round(sliceLocList1), round(sliceLocList2));


%% Get final file lists to use

% Get slices to use in dir1
for i = 1:length(fileList1)
    if ~any(slicesCommon == round(sliceLocList1(i)))
        fileListOut1{i} = [];
    end
end

% Get slices to use in dir2
for i = 1:length(fileList2)
    if ~any(slicesCommon == round(sliceLocList2(i)))
        fileListOut2{i} = [];
    end
end

% Remove empty cells
fileListOut1 = fileListOut1(~cellfun('isempty',fileListOut1));
fileListOut2 = fileListOut2(~cellfun('isempty',fileListOut2));