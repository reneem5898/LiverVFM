function VFM_multiSlice_3D(modelDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caluclate 3D VFM isotropic stiffness for model data
%
% Written by: Renee Miller
% 17 October 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all

addpath('D:\Renee\UoA\MATLAB\HyperelasticModels');
addpath('D:\Renee\UoA\MATLAB\VFM\Optimised_VFM');

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
RECALC = 1;

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
load nodes.mat;

% Each element = 1 row with eight node numbers
%elems = load('elems.txt');
load elems.mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create region contour

for i = 1:numSlices
    
    [c, ~] = contour(squeeze(masks(:,:,i)));
    indices = find(c(1,:) < 1);
    regionContour = c(:,indices(1)+1:indices(2)-1)';

end

% Get nodes and elements inside the mask - slow
regionMeshFile = sprintf('%s/regionMesh.mat', outDir);

if ~exist(regionMeshFile, 'file') || RECALC
    
    % Get nodes which are in the mask
    %Shape = alphaShape(regionContour(:,2), regionContour(:,3), 'HoleThreshold', 100);
    [cn, ~] = inpoly(nodes(:,2:3), regionContour);
    regionNodes = [cn cn cn cn] .* nodes;
    
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
    
    % Interpolate displacements at nodes
    U = zeros(length(nodes)*3,1);
    for dim = 1:3
   
        % Displacements for current slice
        u = squeeze(fha(dim,slice,:,:));
        
        % Interpolate displacements
        u_r = interp2(xx, yy, real(u), nodes(:,2), nodes(:,3), 'linear');
        u_i = interp2(xx, yy, imag(u), nodes(:,2), nodes(:,3), 'linear');
        
        U(dim:3:end) = (u_r + 1i*u_i)*1e7;
        U(isnan(U)) = 0;
    
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate shear stiffness for region of interest only
    
    % Calculate the numeric virtual field
    disp('Calculating the numeric virtual field...')
    [uVF, eta, strain1] = numericVF_Iso(U, regionNodes, regionElems, nodes(:,1), [], 'C3D8', GaussPoints);
    
    % Calculate shear modulus
    disp('Calculating the shear modulus...')
    [~, fg, b, ~] = solveIsoVFM_C3D8(U, uVF, rho, omega, regionNodes, regionElems, nodes(:,1));
    
    % Calculate complex shear modulus
    shear = b/fg
    
    % Calculate normalised sensitivity value = eta/G
    nG = calcNormSensitivity(eta, shear);
    
    % Save results
    save(sprintf('%s/shearResult_slice%d.mat', outDir, slice), 'shear', 'uVF', 'eta', 'nG');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    %plotVF_2D(uVF, U, regionNodes, outDir);
    
    
end
