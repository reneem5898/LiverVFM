function liverVFM_multiSlice(volunteerDir, sequence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caluclate VFM isotropic stiffness for liver data
%
% Written by: Renee Miller
% 9 October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define constants

% Element type used for VFM analysis
% Choices: C3D8 (most accurate), C3D8R (faster, less accurate) or C3D8F
elemType = 'C3D8';

% Level of Gauss Point "refinement" within each element
% If C3D8, GaussPoints = [] (defined within functions and differs for deviatoric versus dilatational stress calculations)
% If C3D8R, GaussPoints = [] (always 1, defined within functions)
% If C3D8F, GaussPoints = 1, 2, or 3. 1 = single Gauss Point at centroid, 2 = 8 Gauss points, 3 = 27 Gauss Points
GaussPoints = [];

% Density = kg/mm^3
rho = 1.0e-6; % Density of water

% Frequency (Hz)
f = 60.1; %%%%%%%%%%%%%%%%% change as needed %%%%%%%%%%%%%%%%%%%%%%
omega = f * 2 * pi; % angular frequency - used in wave equation

% Define directories:

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
tmp3 = strsplit(sequence, '_');
if length(tmp3) > 1
    outDir = sprintf('%s/%s_%s_results', volunteerDir, tmp2{1}, tmp3{end});
else
    outDir = sprintf('%s/%s_results', volunteerDir, tmp2{1});
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Recalculate or not
RECALC = 0;

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
load nodes.mat;

% Each element = 1 row with eight node numbers
%elems = load('elems.txt');
load elems.mat;


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
regionMeshFile = sprintf('%s/regionMesh.mat', outDir);

if ~exist(regionMeshFile, 'file') || RECALC
    
    % Get nodes which are in the mask
    [cn, ~] = inpoly(nodes(:,2:3), imageContourCoords);
    liverNodes = [cn cn cn cn] .* nodes;
    
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
    save(regionMeshFile, 'liverElems', 'liverNodes');
    
else
    load(regionMeshFile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through image slices and calculate shear modulus

for slice = 1:size(phaseMat,2)
    
    % Display current slice number
    fprintf('Slice #: %d', slice);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get image coordinates for liver MRE data
    
    disp('Interpolating image data at nodal coordinates...');
    
    % Create meshgrid of image pixel coordinates
    [xx, yy, zz] = meshgrid(linspace(0.5,res-0.5,res), linspace(0.5,res-0.5,res), linspace(0, 5, 6));
    
    % Displacements at each pixel (must be a faster/cleaner way of doing this)
    uImg = zeros(size(fha,2), size(fha,3), size(xx,3));
    for i = 1:size(xx,3) % x, y, z
        uImg(:,:,i) = squeeze(fha(slice,:,:));
    end
    
    % Interpolate displacements at nodal coordinates
    Uy_r = interp3(xx, yy, zz, real(uImg), nodes(:,2), nodes(:,3), nodes(:,4), 'linear');
    Uy_i = interp3(xx, yy, zz, imag(uImg), nodes(:,2), nodes(:,3), nodes(:,4), 'linear');
    Uy = Uy_r + 1i*Uy_i;
    
    % Turn NaNs into 0's
    Uy(isnan(Uy)) = 0;
    
    % Create list of displacements with zeros for x and y directions
    U = zeros(length(Uy)*3,1);
    U(2:3:end) = Uy;
    
    % Plot interpolated virtual displacement field - check
    %plotVF(U, U, nodes, outDir);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate shear stiffness for liver region only - no subzones
    
    % Final nodes and elements to use:
    % nodes: liverNodes
    % elements: liverElems
    
    % Calculate the numeric virtual field
    disp('Calculating the numeric virtual field...')
    [uVF, eta, strain1] = numericVF_Iso(U, liverNodes, liverElems, nodes(:,1), [], elemType, GaussPoints);
    
    % Calculate shear modulus
    disp('Calculating the shear modulus...')
    
    if strcmp(elemType, 'C3D8R')
        % Uniform strain elements - C3D8R
        [fk, fg, b, strain2] = solveIsoVFM_C3D8R(U, uVF, rho, omega, liverNodes, liverElems, nodes(:,1));
        
    elseif strcmp(elemType, 'C3D8')
        % Selectively reduced integration type element - C3D8
        [fk, fg, b, strain2] = solveIsoVFM_C3D8(U, uVF, rho, omega, liverNodes, liverElems, nodes(:,1));
        
    else
        % Fully integrated element
        [fk, fg, b, strain2] = solveIsoVFM_C3D8F(U, uVF, rho, omega, liverNodes, liverElems, nodes(:,1), GaussPoints);
        
    end
    
    % Calculate complex shear modulus
    shear = b/fg;
    magShear = abs(shear)
    
    % Calculate normalised sensitivity value = eta/G
    nG = calcNormSensitivity(eta, shear);
    
    % Save results
    save(sprintf('%s/shearResult_slice%d.mat', outDir, slice), 'shear', 'uVF', 'eta', 'fk', 'fg', 'b', 'strain2', 'nG', 'magShear');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    plotVF(uVF, U, liverNodes, outDir);
    
    
end


