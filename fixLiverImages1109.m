function fixLiverImages1109(fldr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix naming of folders for liver images solely in volunteer: 20171109
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date modified: 19 February 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List file names (without . and ..)
All = dir(fldr);
filenames = {All.name};
imageFolders = filenames([All.isdir]);
imageFolders(ismember(imageFolders,{'.','..'})) = [];

% Loop through file names
for i = 1:length(imageFolders)
    
    % Get current image file
    currentFolder = sprintf('%s/%s', fldr, imageFolders{i});
    disp(currentFolder)
    
    % Create new folder name - without number at the end
    newFolder = currentFolder(1:(end-5));
    if ~exist(newFolder, 'dir')
        mkdir(newFolder);
    end
    
    % Move files in currentFolder to newFolder
    movefile(sprintf('%s/*', currentFolder), newFolder);
    
    % Remove previous folder name
    rmdir(currentFolder);
end