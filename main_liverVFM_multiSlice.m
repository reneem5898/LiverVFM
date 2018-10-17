%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caluclate VFM isotropic stiffness for liver data
%
% Written by: Renee Miller
% 9 October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

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
f = 60.1;
omega = f * 2 * pi; % angular frequency - used in wave equation

% Define directories:
% Patient directory - CHANGE THIS FOR EACH NEW PATIENT
patientDir = 'P:\Data\Liver\MRE_20171106_1';

% Phase images
phaseFldr = 'epiMRE923_tra_bh_100_4slc_fract_P_Wave'; % may need to change the number at the end of this directory
phaseDir = sprintf('%s/', patientDir, phaseFldr);
% Magnitude images
magFldr = 'epiMRE923_tra_bh_100_4slc_fract_Mag'; % may need to change the number at the end of this directory
magDir = sprintf('%s/%s', patientDir, magFldr);
% Output directory
tmp = strsplit(magFldr, '_');
outDir = sprintf('%s/%s_results', patientDir, tmp{1});
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
resultsFile = sprintf('%s/results.mat', outDir);

if ~exist(resultsFile, 'file') || RECALC
    
    % Get complex first harmonic amplitude and mask
    [fha, mask, pixelSpacing, avgMag, phaseMat, dc] = getDisp_MS(res, phaseDir, magDir, outDir);
    save(sprintf('%s/results.mat', outDir), 'fha', 'mask', 'pixelSpacing', 'avgMag', 'phaseMat', 'dc');
    
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
%% Loop through  image slices

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
    for i = 1:size(xx,3)
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
    plotVF(U, U, nodes, outDir);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate shear stiffness for liver region only - no subzones
    
    % Manually segment a region within the liver for VFM stiffness estimation
    disp('Creating liver contour...');
    liverContourFile = sprintf('%s/liverRegionContour.txt', outDir); % Filename
    %contourCoords = get_contour(liverContourFile, squeeze(avgMag(slice,:,:)), 'Outline the region of the liver'); % If file already exists, read in contour. Else, ask user to contour region.  
    contourCoords = getLiverContour(patientDir);
    
    % Plot contour on magnitude image of liver
    FH2 = figure;
    imagesc(double(squeeze(avgMag(slice,:,:))))
    caxis([0 80])
    hold on
    colormap(gray)
    scatter(contourCoords(:,1), contourCoords(:,2), 'r.')
    axis off
    axis image
    title('Region for stiffness estimation', 'FontSize', 12)
    saveas(FH2, sprintf('%s/region-mag.png', outDir));
    close(FH2)
    
    % Plot contour on magnitude image of liver
    FH3 = figure('position', [100, 100, 1800, 800]);
    
    subplot(1,2,1)
    imagesc(real(squeeze(fha(slice,:,:))))
    caxis([-100 100])
    hold on
    scatter(contourCoords(:,1), contourCoords(:,2), 'ko', 'filled')
    axis off
    axis image
    title('Real(FHA)', 'FontSize', 16)
    
    subplot(1,2,2)
    imagesc(imag(squeeze(fha(slice,:,:))))
    caxis([-100 100])
    hold on
    scatter(contourCoords(:,1), contourCoords(:,2), 'ko', 'filled')
    axis off
    axis image
    title('Imag(FHA)', 'FontSize', 16)
    saveas(FH3, sprintf('%s/region-fha.png', outDir));
    close(FH3)
    
    % Get nodes which are in the mask
    [cn, ~] = inpoly(nodes(:,2:3), contourCoords);
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
    shear(slice) = b/fg;
    magShear(slice) = abs(shear)
    
    % Calculate normalised sensitivity value = eta/G
    nG(slice) = calcNormSensitivity(eta, shear);
    
    % Save results
    save(sprintf('%s/shearResult.mat', outDir), 'shear', 'uVF', 'eta', 'fk', 'fg', 'b', 'strain1', 'strain2', 'nG');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    plotVF(uVF, U, liverNodes, outDir);
    
    
end
