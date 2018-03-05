function sortLiverImages(fldr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to sort images from liver scan by the name of the protocol
% 
% Renee Miller (rmil520@aucklanduni.ac.nz)
% 7 November 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List file names (without . and ..)
AllFiles = dir(fldr);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Loop through file names
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', fldr, imageFiles{i});
    
    % Save image info temporarily
    imgInfo = dicominfo(filename); 
    
    % Protocol name
    %protocolName = imgInfo.ProtocolName;
    protocolName = imgInfo.SeriesDescription;
    
    % Make a directory with the protocol name
    protocolDir = sprintf('%s/%s', fldr, protocolName);
    if ~exist(protocolDir, 'dir')
        mkdir(protocolDir)
        disp(protocolName)
    end
    
    % Move file to protocol directory
    %copyfile(filename, sprintf('%s/', protocolDir));
    movefile(filename, sprintf('%s/', protocolDir));
    
end