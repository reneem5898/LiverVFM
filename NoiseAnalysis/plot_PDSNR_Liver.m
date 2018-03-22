function [pdsnr1, meanPDSNR, stdevPDSNR] = plot_PDSNR_Liver(pdsnr, mask, figDir)

% ------------------------------------------------------------------------
% This function plots the pdsnr weighting values for each slice in each
% direction
%
% Inputs: 1) pdsnr matrix
%         2) Image masks
%         3) Folder where to save images
%
% Outputs: 1) pdsnr1 - pdsnr vector
%          2) meanPDSNR - mean pdsnr value for region
%          3) stdevPDSNR - standard deviation of pdsnr for region
%
% Written by: Renee Miller
% Date: 21 March 2018
% 
% Called from: runAll_LiverNoise.m
% ------------------------------------------------------------------------

% Directory where to save images
figDirWeight = strcat(figDir, '/Weighting/');
if ~exist(figDirWeight, 'dir')
    mkdir(figDirWeight);
end

minPDSNR = floor(min(pdsnr(:)/10))*10;
maxPDSNR = ceil(max(pdsnr(:)/10))*10;

sz = size(pdsnr);
res = sz(end);

% Crop the image around the myocardium
ROW_START = 1;
ROW_END = res;
COL_START = 1;
COL_END = res;

% Create plots of weights for each slice and direction
for i = 1:sz(1)
    
    % Get mask for the particular slice
    I = squeeze(~mask(i,ROW_START:ROW_END,COL_START:COL_END));
    % Create white background for images
    white = cat(3, ones(res,res), ones(res,res), ones(res,res));
    white = white(ROW_START:ROW_END,COL_START:COL_END,:);
    
    FH = figure('Position', [100, 100, 600, 400]);
    imagesc(squeeze(pdsnr(i,ROW_START:ROW_END,COL_START:COL_END)))
    hold on
    h = imshow(white);
    hold off
    set(h, 'AlphaData', I)
    title(sprintf('Slice %d', i), 'FontSize', 12)
    colorbar
    caxis([minPDSNR maxPDSNR])
    axis image
    
    figName = sprintf('%sPDSNR_weights_Slice%d', figDirWeight, i);
    saveas(FH,figName,'png')
end

% Put values of pdsnr from each independent direction into vectors to
% create histograms
pdsnr1 = mask.*pdsnr;
pdsnr1 = pdsnr1(:);
pdsnr1(pdsnr1==0) = [];
meanPDSNR = mean(pdsnr1);
stdevPDSNR = std(pdsnr1);
medianPDSNR = median(pdsnr1);

FH = figure('Position', [100, 100, 600, 400]);
hist(pdsnr1,20)
hold on
%axis([0 30 0 250])
title('PDSNR Histogram', 'FontSize', 12)
annotation('textbox', [0.3 0.75 0.25 0.125], 'String', ...
    {sprintf('Mean: %.2f\nMedian: %.2f', meanPDSNR, medianPDSNR)},'FontSize', 12)
set(gca, 'FontSize', 12)

figName = sprintf('%sPDSNR-Histogram', figDirWeight);
saveas(FH,strcat(figName, '.png'));

close all