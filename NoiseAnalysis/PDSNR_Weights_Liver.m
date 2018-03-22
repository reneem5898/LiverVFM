function [pdsnr, mask] = PDSNR_Weights_Liver(fha, AvgImgs, regionContour)

% ------------------------------------------------------------------------
% PDSNR - calculates Phase Difference SNR from the magnitude and phase (displacements) MRE images
% PDSNR = abs(FHA)/((1/sqrt(2)*asin(1/MagSNR))
%
% Inputs: 1) fha - first harmonic amplitudes
%         2) AvgImgs - magnitude images averaged over the phase offsets
%         3) masks - used to create folder for saving images
%
% Outputs: 1) pdsnr - matrix of phase difference snr values 
%
% Written by: Renee Miller
% Date: 21 March 2018
%
% Called in MRE_Noise_Liver.m
%
% ------------------------------------------------------------------------

% Create region mask from contour
m = poly2mask(regionContour(:,1), regionContour(:,2), size(AvgImgs,2), size(AvgImgs,3));
mask = zeros(size(AvgImgs));
for i = 1:size(AvgImgs,1);
    mask(i,:,:) = m;
end

% Mask average images
maskedAvgImgs = mask.*AvgImgs;

% Put image data into a vector
signal = maskedAvgImgs(:);

% remove all zero values from masking
signal(signal==0) = [];

% compute the mean of the signal values
MeanSignal = mean(signal);
% compute standard deviation of signal region
StdSignal = std(signal);
% compute the magnitude SNR
MagSNR = MeanSignal/StdSignal;

% compute phase error - See ppt on PDSNR calculation for explanation
PhaseError = (1.0/sqrt(2))*asin(1.0/MagSNR);

% Mask displacements
MaskFHA = (mask.*abs(fha))*2; %Multiply by 2 because displacements saved in Disp are mean to peak. For noise calculation, use peak to peak displacements.

% get pdsnr matrix
menc = 1;
pdsnr = MaskFHA/(menc*PhaseError)/2048;

