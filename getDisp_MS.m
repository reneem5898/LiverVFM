function [fha, maskAbdomen, pixelSpacing, avgMag, phaseMat, magMat, dc] = getDisp_MS(res, phaseFldr, magFldr, outDir)

%% getDisp - gets complex fha from liver images
% Inputs: 1) Image resolution
%         2) Folder containing phase images
%         3) Folder containing magnitude images
%
% Outputs: 1) Disp: displacement data
%          2) masks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put MRE image data into one variable - ImgMat

% Get list of all files in directory
AllFiles = dir(phaseFldr);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Create array to store data (phase offsets x 256 pixels x 256 pixels)
phaseMat = zeros(length(imageFiles), res, res);

% Slice locations
sl = [];

% Save phase images in matrix
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', phaseFldr, imageFiles{i});
    
    % Save image temporarily
    img = dicomread(filename);
    
    % Save image info temporarily
    info = dicominfo(filename);
    
    % Get slice height
    sl = [sl; info.SliceLocation];
    
    % Load image data
    phaseMat(i,1:size(img,1),1:size(img,2)) = img;

end

% Get only unique slice locations
sliceHeights = unique(sl, 'stable');

% Get number of slices
numSlices = length(sliceHeights);

% Get number of offsets
numOffsets = size(phaseMat,1)/numSlices;

% Reshape phaseMat = phase offsets x slices x res x res
phaseMat = reshape(phaseMat, [numOffsets, numSlices, res, res]); %% works because images were collected in order - one slice at a time

% Get pixel spacing
info = dicominfo(filename);
pixelSpacing = info.PixelSpacing;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate complex FHA from MRE images

disp('Taking fft of phase data...');

% FFT along first dimension, which are the phase offsets
fourier = fft(phaseMat,[],1); 

% 1st location of fft corresponds to 0th harmonic or DC, 2nd location in
% fft corresponds to 1st harmonic, and so on for n harmonics. we need the first harmonics
fha = squeeze(fourier(2,:,:,:)); 

% Save the DC component
dc = squeeze(fourier(1,:,:,:));

fha = fha*2/size(phaseMat,1); %DIFFERENT FROM ARUN - in FE model, need to apply median to peak amplitude 
%Disp_Ax = squeeze(fha)*2*pi*menc/4096; % convert from pixel units to displacement (um) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put magnitude image data into one variable for masking

% Get list of all files in directory
AllFiles = dir(magFldr);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Initialise slice positions - magnitude
sl_mag = [];

% Create array to store data (phase offsets x 256 pixels x 256 pixels)
magMat = zeros(length(imageFiles), res, res);

% Save phase images in matrix
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', magFldr, imageFiles{i});
    
    % Save image temporarily
    img = dicomread(filename);
    
    % Save image info temporarily
    info = dicominfo(filename);
    
    % Get slice height
    sl_mag = [sl_mag; info.SliceLocation];
    
    % Load image data
    magMat(i,1:size(img,1),1:size(img,2)) = img;

end

% Reshape magMat = phase offsets x slices x res x res
magMat = reshape(magMat, [length(imageFiles)/numSlices, numSlices, res, res]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask fha

% Take average of phase offset magnitude images
avgMag = squeeze(mean(magMat,1)); 

% Initialise mask
maskAbdomen = zeros(size(magMat,2), res, res);

% Loop through slices
for i = 1:size(magMat,2)
    
    % Create mask using region growing algorithm
    % First part starts at the beginning (1,1) and masks out everything outside of the abdomen
    % Second part starts at last pixel (256, 256) and masks out part of image that is "too big" for original image
    maskAbdomen(i,:,:) = ~regiongrowing(squeeze(avgMag(i,:,:)), 1, 1, 5.9).*~regiongrowing(squeeze(avgMag(i,:,:)), 256, 256);
    % Value of 5.9 controls the gradient between regions (found by trial and error) -- WARNING: MASKING NOT PERFECT
    
    % Plot mask
    FH = figure;
    imagesc(squeeze(maskAbdomen(i,:,:)))
    hold on
    colormap gray
    axis image
    title(sprintf('Slice %d', i));
    saveas(FH, sprintf('%s/mask_slice%d.png', outDir, i));
    
end

% Mask fha
fha = fha.*maskAbdomen;

