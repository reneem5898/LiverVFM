function liverVFM_multiSlice_2D(volunteerDir, sequence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caluclate VFM isotropic stiffness for liver data
%
% Written by: Renee Miller
% 6 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define constants

% Level of Gauss Point "refinement" within each element (number of integration 
% points in each direction)
GaussPoints = 2;

% Density = kg/mm^3
rho = 1.0e-6; % Density of water

% Frequency (Hz)
f = 60.1; %%%%%%%%%%%%%%%%% change as needed %%%%%%%%%%%%%%%%%%%%%%
omega = f * 2 * pi; % angular frequency - used in wave equation

% Recalculate or not
RECALC = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define directories:

% Phase images
file_spec = sprintf('%s/*%s*P_Wave', volunteerDir, sequence);
tmp = dir(file_spec);
tmpName = {tmp.name};
phaseDir = sprintf('%s/%s', volunteerDir, tmpName{1});

% Magnitude images
file_spec = sprintf('%s/*%s*Mag', volunteerDir, sequence);
tmp = dir(file_spec);
tmpName = {tmp.name};
magDir = sprintf('%s/%s', volunteerDir, tmpName{1});

% Output directory
tmp = strsplit(magDir, '/');
tmp2 = strsplit(tmp{end}, '_');
tmp3 = strsplit(sequence, '*');
if length(tmp3) > 1
    outDir = sprintf('%s/%s_%s_results', volunteerDir, tmp2{1}, tmp3{end});
else
    outDir = sprintf('%s/%s_results', volunteerDir, tmp2{1});
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load image data and calculate complex FHA

% Defines image size
% Somewhat arbitrary - does not have to be the actual image size but must be
% AT LEAST as large as the largest image dimension. Currently: set to be the same size
% as the FE mesh (256 x 256). IF YOU NEED TO CHANGE THIS, LET ME KNOW. WE'LL NEED
% TO CREATE ANOTHER MESH.
res = 256;

% Results file
resultsFile = sprintf('%s/displacements.mat', outDir);

if ~exist(resultsFile, 'file') || RECALC
    
    % Get complex first harmonic amplitude and mask
    [fha, mask, pixelSpacing, avgMag, phaseMat, magMat, dc] = getDisp_MS(res, phaseDir, magDir, outDir);
    save(resultsFile, 'fha', 'mask', 'pixelSpacing', 'avgMag', 'phaseMat', 'magMat', 'dc');
    
    % Plot fft fit for one pixel (70,100) - just for visualisation
    plotLiverPhaseOffsets(squeeze(phaseMat(:,1,:,:)), squeeze(fha(1,:,:)), outDir, 70, 100);
    
    % Create gif of wave propagation
    makeWaveGif(phaseDir, outDir);
    
else
    load(resultsFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model data

disp('Loading model node and element data...');

% Node coordinates
%nodes = load('nodeCoords.txt');
load nodes_2D.mat;

% Each element = 1 row with eight node numbers
%elems = load('elems.txt');
load elems_2D.mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load liver contour and get elements and nodes inside the region 
% ROI - same for each slice and sequence

% Load ROI from text file
liverContourFile = sprintf('%s/RoIContour.txt', volunteerDir);
contourCoords = load(liverContourFile);

% Get transformation matrices to convert contour coordinates into coordinates
% for current image
numSlices = size(mask, 1);
[contour2magnet, magnet2image] = getTransform(volunteerDir, magDir, numSlices);
contour2magnet = squeeze(contour2magnet(1,:,:)); % Since contour is the same for each slice, we can just take the first set of transformation matrices
magnet2image = squeeze(magnet2image(1,:,:));

% Initialise new contour coordinates vector
imageContourCoords = zeros(size(contourCoords));

% Convert contours
for i = 1:length(contourCoords)
    ROI_mag_coord = contour2magnet * [contourCoords(i,1); contourCoords(i,2); 0; 1];
    tmp = round(magnet2image * ROI_mag_coord);
    imageContourCoords(i,:) = tmp(1:2);
end

% Plot contour on magnitude image of liver
FH2 = figure;
imagesc(double(squeeze(avgMag(1,:,:))))
caxis([0 80])
hold on
colormap(gray)
scatter(imageContourCoords(:,1), imageContourCoords(:,2), 'r.')
axis off
axis image
title('Region for stiffness estimation', 'FontSize', 12)
saveas(FH2, sprintf('%s/region-mag.png', outDir));
close(FH2)

% Plot contour on magnitude image of liver
FH3 = figure('position', [100, 100, 1800, 800]);

subplot(1,2,1)
imagesc(real(squeeze(fha(1,:,:))))
caxis([-100 100])
hold on
scatter(imageContourCoords(:,1), imageContourCoords(:,2), 'ko', 'filled')
axis off
axis image
title('Real(FHA)', 'FontSize', 16)

subplot(1,2,2)
imagesc(imag(squeeze(fha(1,:,:))))
caxis([-100 100])
hold on
scatter(imageContourCoords(:,1), imageContourCoords(:,2), 'ko', 'filled')
axis off
axis image
title('Imag(FHA)', 'FontSize', 16)
saveas(FH3, sprintf('%s/region-fha.png', outDir));
close(FH3)

% Get nodes and elements inside the mask - slow
regionMeshFile = sprintf('%s/regionMesh_2D.mat', outDir);

if ~exist(regionMeshFile, 'file') || RECALC
    
    % Get nodes which are in the mask
    [cn, ~] = inpoly(nodeCoords(:,2:3), imageContourCoords);
    liverNodes = [cn cn cn] .* nodeCoords;
    
    % Remove all rows with zeros
    liverNodes( ~any(liverNodes,2), : ) = [];
    
    % Get elements which contain liver nodes (slow) -- HOW TO MAKE FASTER?
    disp('Getting mesh in selected region...');
    elemsRemove = ones(size(elems));
    for i = 1:size(elems,1)
        for j = 1:(size(elems,2)-1)
            if ~any(liverNodes(:,1) == elems(i,j+1))
                elemsRemove(i,j+1) = 0;
                j = size(elems,2)-1;
            end
        end
    end
    
    tmp = elemsRemove .* elems;
    liverElems = tmp(all(tmp,2),:); % Get rid of any row with a zero in it
    
    % Put node coordinates into units of mm
    liverNodes(:,2) = liverNodes(:,2)*pixelSpacing(1); % x
    liverNodes(:,3) = liverNodes(:,3)*pixelSpacing(2); % y
    
    % Save element and node data in ROI
    save(regionMeshFile, 'liverElems', 'liverNodes', 'imageContourCoords');
    
else
    load(regionMeshFile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through image slices and calculate shear modulus

for slice = 1:size(phaseMat,2)
    
    % Display current slice number
    disp(sprintf('Slice #: %d', slice));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get image coordinates for liver MRE data
    
    disp('Interpolating image data at nodal coordinates...');
    
    % Create meshgrid of image pixel coordinates
    [xx, yy] = meshgrid(linspace(0.5,res-0.5,res), linspace(0.5,res-0.5,res));
    
    % Current image slice displacements
    uImg = squeeze(fha(slice,:,:));
        
    % Interpolate displacements at nodal coordinates
    Uy_r = interp2(xx, yy, real(uImg), nodeCoords(:,2), nodeCoords(:,3), 'linear');
    Uy_i = interp2(xx, yy, imag(uImg), nodeCoords(:,2), nodeCoords(:,3), 'linear');
    Uy = Uy_r + 1i*Uy_i;
    
    % Turn NaNs into 0's
    Uy(isnan(Uy)) = 0;
    
    % Create list of displacements with zeros for x direction
    U = zeros(length(Uy)*2,1);
    U(2:2:end) = Uy;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate shear stiffness for liver region only
    
    % Final nodes and elements to use:
    % nodes: liverNodes
    % elements: liverElems
    
    % Calculate the numeric virtual field
    disp('Calculating the numeric virtual field...')
    [uVF, eta, strain1] = numericVF_Iso2D(U, liverNodes, liverElems, nodeCoords(:,1), [], GaussPoints);
    
    % Calculate shear modulus
    disp('Calculating the shear modulus...')
    [fk, fg, b, strain2] = solveIsoVFM_C2D4(U, uVF, rho, omega, liverNodes, liverElems, nodeCoords(:,1), GaussPoints);
    
    
    % Calculate complex shear modulus
    shear = b/fg;
    magShear = abs(shear)
    
    % Calculate normalised sensitivity value = eta/G
    nG = calcNormSensitivity(eta, shear);
    
    % Save results
    save(sprintf('%s/shearResult_slice%d_2D.mat', outDir, slice), 'shear', 'uVF', 'eta', 'fk', 'fg', 'b', 'strain2', 'nG', 'magShear');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    plotVF_2D(uVF, U, liverNodes, outDir);
    
    
end


