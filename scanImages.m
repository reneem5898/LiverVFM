%close all
folder = 'P:\Data\Liver\MRE_20171106_2\SIEMENS greMRE_tra_p2_bh_128_P_Wave';
%folder = 'P:\Data\Liver\MRE_20171106_1\greMRE_tra_p2_bh_128__Mag';

% Get list of all files in directory
AllFiles = dir(folder);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];



% Slice locations
sl = [];

% Save phase images in matrix
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', folder, imageFiles{i});
    
    % Save image temporarily
    img = dicomread(filename);
    
    % Slice height
    tmp = dicominfo(filename);
    sliceLoc = tmp.SliceLocation;
    
    figure(1)
    imagesc(img)
    hold on
    title(sprintf('Image: %d\nSlice Location: %.3f', i, sliceLoc))
    %colormap(gray)
    pause
    close(1)
    
end