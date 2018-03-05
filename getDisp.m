function [fha, maskAbdomen, pixelSpacing, avgMag, phaseMat, dc] = getDisp(res, phaseFldr, magFldr)

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

% Save phase images in matrix
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', phaseFldr, imageFiles{i});
    
    % Save image temporarily
    img = dicomread(filename);
    
    % Load image data
    phaseMat(i,1:size(img,1),1:size(img,2)) = img;

end

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
fha = squeeze(fourier(2,:,:)); 

% Save the DC component
dc = squeeze(fourier(1,:,:,:,:));

fha = fha*2/size(phaseMat,1); %DIFFERENT FROM ARUN - in FE model, need to apply median to peak amplitude 
%Disp_Ax = squeeze(fha)*2*pi*menc/4096; % convert from pixel units to displacement (um) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put magnitude image data into one variable for masking

% Get list of all files in directory
AllFiles = dir(magFldr);
imageFiles = {AllFiles.name};
imageFiles(ismember(imageFiles,{'.','..'})) = [];

% Create array to store data (phase offsets x 256 pixels x 256 pixels)
magMat = zeros(length(imageFiles), res, res);

% Save phase images in matrix
for i = 1:length(imageFiles)
    
    % Get current image file
    filename = sprintf('%s/%s', magFldr, imageFiles{i});
    
    % Save image temporarily
    img = dicomread(filename);
    
    % Load image data
    magMat(i,1:size(img,1),1:size(img,2)) = img;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask fha

% Take average of phase offset magnitude images
avgMag = squeeze(mean(magMat,1)); 

% Create mask using region growing algorithm
% First part starts at the beginning (1,1) and masks out everything outside of the abdomen
% Second part starts at last pixel (256, 256) and masks out part of image that is "too big" for original image
maskAbdomen = ~regiongrowing(avgMag, 1, 1, 5.9).*~regiongrowing(avgMag, 256, 256);
% Value of 5.9 controls the gradient between regions (found by trial and error) -- WARNING: MASKING NOT PERFECT

% Mask fha
fha = fha.*maskAbdomen;

% Plot mask
figure
imagesc(maskAbdomen)
hold on
colormap gray
axis image

