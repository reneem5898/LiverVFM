function [contour2magnet, magnet2image] = getTransform(volunteerDir, magFldr, numSlices)

%% getTransform - gets transformation matrices between ROI contour and current image
% 
% Inputs: 1) volunteerDir - directory for current volunteer/patient
%         2) magFldr - magnitude folder for current sequence
%         3) Number of image slices
%
% Outputs: 1) contour2magnet - transformation matrix between product image coordinates
%             to patient (magnet) coordinates
%          2) magnet2image - transformation matrix between patient coordinates
%             to image coordinates for current sequence (not necessarily 
%             the product sequence) 
%
% Note: matrices are used for converting contours from one MRE sequence
% to another MRE sequence
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date modified: 19 February 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put MRE image data into one variable - ImgMat

% Get list of all files in magnitude directory
AllFiles = dir(magFldr);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Product sequence directory - to get contour2magnet transformation matrices
% Folder names hard coded for Siemens product sequence - bad
productDir = sprintf('%s/SIEMENS greMRE_tra_p2_bh_128_Mag', volunteerDir);
if ~exist(productDir, 'dir')
    productDir = sprintf('%s/SIEMENS_GREMRE_TRA_P2_BH_128_MAG', volunteerDir); %20171109
end

% Get list of all files in magnitude directory for product sequence
AllFilesProduct = dir(productDir);
imageFilesProduct = {AllFilesProduct.name};
imageFilesProduct(ismember(imageFilesProduct,{'.','..'})) = [];

% Create array to store data (phase offsets x 256 pixels x 256 pixels)
magnet2image = zeros(numSlices, 4, 4);
contour2magnet = zeros(numSlices, 4, 4);

% Slice locations
sl = zeros(length(imageFiles), 2);

% Check to see that there are the same number of images in both the product
% and wip folders
if ~(length(imageFilesProduct)==length(imageFiles))
    [imageFiles, imageFilesProduct] = findMatchingSlices(magFldr, productDir);
end

% Get slice locations
for i = 1:length(imageFiles)
    
    % Get current image files
    wip_filename = sprintf('%s/%s', magFldr, imageFiles{i});
    seq_filename = sprintf('%s/%s', productDir, imageFilesProduct{i});
       
    % Save image info temporarily
    info_wip = dicominfo(wip_filename);
    info_product = dicominfo(seq_filename);
    
    % Get slice height
    sl(i,:) = [info_wip.SliceLocation info_product.SliceLocation];

end

% Round slice heights
sl_round = round(sl);

% Check to see if images were acquired at the same heights for each sequence - approximately
err = sl_round(:,1) - sl_round(:,2);
if err
    disp('WARNING: Images were not acquired at the same slice location as the product sequence.');
end

% Get only unique slice locations
sliceHeights = unique(sl_round(:,1), 'stable');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save transformation matrices

for i = 1:length(imageFiles)
    
    % Get current image files
    wip_filename = sprintf('%s/%s', magFldr, imageFiles{i});
    product_filename = sprintf('%s/%s', productDir, imageFilesProduct{i});
       
    % Save image info temporarily
    info_wip = dicominfo(wip_filename);
    info_product = dicominfo(product_filename);
    
    % Get slice height
    sliceLoc_wip = round(info_wip.SliceLocation);
    sliceLoc_product = round(info_product.SliceLocation);
       
    % Save MRI data
    mri_data_wip = MRI(wip_filename);
    mri_data_product = MRI(product_filename);
    
    % Save transformation matrices
    contour2magnet(find(sliceHeights == sliceLoc_wip),:,:) = mri_data_product.TI2P;
    magnet2image(find(sliceHeights == sliceLoc_product),:,:) = mri_data_wip.TP2I;

end




