%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caluclate 2D VFM isotropic stiffness for model data
%
% Written by: Renee Miller
% 6 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all

addpath('D:\Renee\UoA\MATLAB\HyperelasticModels');
addpath('D:\Renee\UoA\MATLAB\VFM\Optimised_VFM');

modelDir = 'P:\UA - PhD\FEM\Frequency Analysis\PreStretch\beam_iso_loads\Comp_Y\SSD0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define constants

% Level of Gauss Point "refinement" within each element (number of integration 
% points in each direction)
GaussPoints = 2;

% Density = kg/mm^3
rho = 1.0e-6; % Density of water

% Frequency (Hz)
f = 80; %%%%%%%%%%%%%%%%% change as needed %%%%%%%%%%%%%%%%%%%%%%
omega = f * 2 * pi; % angular frequency - used in wave equation

% Recalculate or not
RECALC = 0;

% Number of short axis slices to create
numSlices = 6;

% Pixel spacing (mm/pixel)
pixelSpacing = 1;

% Image size
imgSize = 256; % 128 x 128


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make output directory

outDir = sprintf('%s/VFM/', modelDir);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load image matrix and masks

load(sprintf('%s/imgMatrix.mat', modelDir));
load(sprintf('%s/masks.mat', modelDir));


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
%% Create region contour

for i = 1:numSlices
    
    [c, ~] = contour(squeeze(masks(:,:,i)));
    indices = find(c(1,:) < 1);
    regionContour = c(:,indices(1)+1:indices(2)-1)';

end

% Get nodes and elements inside the mask - slow
regionMeshFile = sprintf('%s/regionMesh_2D.mat', outDir);

if ~exist(regionMeshFile, 'file') || RECALC
    
    % Get nodes which are in the mask
    %Shape = alphaShape(regionContour(:,2), regionContour(:,3), 'HoleThreshold', 100);
    [cn, ~] = inpoly(nodeCoords(:,2:3), regionContour);
    regionNodes = [cn cn cn] .* nodeCoords;
    
    % Remove all rows with zeros
    regionNodes( ~any(regionNodes,2), : ) = [];
    
    % Get elements which contain liver nodes (slow) -- HOW TO MAKE FASTER?
    disp('Getting mesh in selected region...');
    elemsRemove = ones(size(elems));
    for i = 1:size(elems,1)
        for j = 1:(size(elems,2)-1)
            if ~any(regionNodes(:,1) == elems(i,j+1))
                elemsRemove(i,j+1) = 0;
                j = size(elems,2)-1;
            end
        end
    end
    
    tmp = elemsRemove .* elems;
    regionElems = tmp(all(tmp,2),:); % Get rid of any row with a zero in it
    
    % Put node coordinates into units of mm
    regionNodes(:,2) = regionNodes(:,2)*pixelSpacing; % x
    regionNodes(:,3) = regionNodes(:,3)*pixelSpacing; % y
    
    % Save element and node data in ROI
    save(regionMeshFile, 'regionElems', 'regionNodes');
    
else
    load(regionMeshFile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate displacements from image matrices

[~, ~, fha] = calcDispModel(ImgMat, 1, 4); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through image slices and calculate shear modulus

for slice = 1:numSlices
    
    % Display current slice number
    fprintf('Slice #: %d', slice);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get image coordinates for liver MRE data
    
    disp('Interpolating image data at nodal coordinates...');
    
    % Create meshgrid of image pixel coordinates
    [xx, yy] = meshgrid(linspace(0.5,imgSize-0.5,imgSize), linspace(0.5,imgSize-0.5,imgSize));
    
    % Calculate shear modulus for displacements in each orthogonal
    % direction separately  
    
    % Initialise variables
    uVF = zeros(length(regionNodes)*2,3);
    eta = zeros(3,3);
    shear = zeros(3,1);
    nG = zeros(3,1);
    
    for dim = 1:3
        
        U = zeros(length(nodeCoords)*2,1);
        
        % Displacements for current slice
        u = squeeze(fha(dim,slice,:,:));
        
        % Interpolate displacements
        u_r = interp2(xx, yy, real(u), nodeCoords(:,2), nodeCoords(:,3), 'linear');
        u_i = interp2(xx, yy, imag(u), nodeCoords(:,2), nodeCoords(:,3), 'linear');
        
        U(1:2:end) = (u_r + 1i*u_i)*1e7;
        U(isnan(U)) = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calculate shear stiffness for region of interest only
        
        % Calculate the numeric virtual field
        disp('Calculating the numeric virtual field...')
        [uVF(:,dim), eta(dim,:), strain1] = numericVF_Iso2D(U, regionNodes, regionElems, nodeCoords(:,1), [], GaussPoints);
                
        % Calculate shear modulus
        disp('Calculating the shear modulus...')
        [fk, fg, b, strain2] = solveIsoVFM_C2D4(U, uVF, rho, omega, regionNodes, regionElems, nodeCoords(:,1), GaussPoints);
        
        % Calculate complex shear modulus
        shear(dim) = b/fg
        
        % Calculate normalised sensitivity value = eta/G
        nG(dim) = calcNormSensitivity(eta(dim,:), shear(dim));
        
    end
    
    % Save results
    save(sprintf('%s/shearResult_slice%d_2D.mat', outDir, slice), 'shear', 'uVF', 'eta', 'nG');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    %plotVF_2D(uVF, U, regionNodes, outDir);
    
    
end
