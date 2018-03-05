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

% Size of subzone - # number of elements in x (and y) dimension per subzone
% bigger number = more free nodes (more accuracy)
% smaller number = faster
subZoneDim = 10;

% Define directories:
% Patient directory - CHANGE THIS FOR EACH NEW PATIENT
patientDir = 'C:/Users/rmil520/Documents/MATLAB/Liver/Volunteer_1/887/';

% Slice location
sliceLoc = 35;

% Phase images
%phaseFldr = 'epiMRE923_tra_bh_100_4slc_fract_P_Wave'; % may need to change the number at the end of this directory
phaseFldr = 'Phase Wave';
phaseDir = sprintf('%s/%s/%d/', patientDir, phaseFldr, sliceLoc); 
% Magnitude images
%magFldr = 'epiMRE923_tra_bh_100_4slc_fract_Mag'; % may need to change the number at the end of this directory
magFldr = 'Magnitude';
magDir = sprintf('%s/%s/%d/', patientDir, magFldr, sliceLoc); 
% Output directory
tmp = strsplit(magFldr, '_');
outDir = sprintf('%s/results/%d/', patientDir, sliceLoc);
%outDir = sprintf('%s/%s_results', patientDir, tmp{1});
if ~exist(outDir, 'dir')
    mkdir(outDir);
end


% % Phase images
% phaseDir = sprintf('%s\\Phase Wave\\%g', patientDir,Slice); % may need to change the number at the end of this directory
% % Magnitude images
% magDir = sprintf('%s\\Magnitude\\%g', patientDir,Slice); % may need to change the number at the end of this directory
% % Stiffness map
% stiffDir = sprintf('%s\\Stiffness Map\\%g', patientDir,Slice);
% % Output directory
% outDir = sprintf('C:\\Users\\jmon961\\Documents\\Joshua\\Inversion Algorithms\\Input_Image_Files\\PatientData\\Patient%g\\SavedFiles%g',Volunteer,Slice);
% if ~exist(outDir, 'dir')
%     mkdir(outDir);
% end

% Calculate heterogeneous (slow) or homogeneous (fast) shear stiffness over region of interest
%area = 'het'; 
area = 'hom'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load image data and calculate complex FHA

% Defines image size
% Somewhat arbitrary - does not have to be the actual image size but must be
% AT LEAST as large as the largest image dimension. Currently: set to be the same size 
% as the FE mesh (256 x 256). IF YOU NEED TO CHANGE THIS, LET ME KNOW. WE'LL NEED
% TO CREATE ANOTHER MESH.
res = 256;

% Get complex first harmonic amplitude and mask
[fha, mask, pixelSpacing, avgMag, phaseMat, dc] = getDisp(res, phaseDir, magDir);
save(sprintf('%s/results.mat', outDir), 'fha', 'mask', 'pixelSpacing', 'avgMag', 'phaseMat', 'dc');

% Plot fft fit for one pixel (70,100) - just for visualisation
plotLiverPhaseOffsets(phaseMat, fha, outDir, 70, 100);

% Create gif of wave propagation
makeWaveGif(phaseDir, outDir);


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
%% Get image coordinates for liver MRE data

disp('Interpolating image data at nodal coordinates...');

% Create meshgrid of image pixel coordinates
[xx, yy, zz] = meshgrid(linspace(0.5,res-0.5,res), linspace(0.5,res-0.5,res), linspace(0, 5, 6));

% Displacements at each pixel (must be a faster/cleaner way of doing this)
uImg = zeros(size(fha,1), size(fha,2), size(xx,3));
for i = 1:size(xx,3)
    uImg(:,:,i) = fha;
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
%% Check if the user wants to calculate the heterogeneous shear moduli over the entire image
% or just in a region of interest

if strcmp(area, 'het') %%% SLOW %%%
    
    % Manually segment a the liver for heterogeneous VFM stiffness estimation
    disp('Creating liver contour...'); 
    liverContourFile = sprintf('%s/liverContour.txt', outDir); % Filename
    contourCoords = get_contour(liverContourFile, avgMag, 'Outline the region of the liver'); % If file already exists, read in contour. Else, ask user to contour region.
    
    % Plot contour on magnitude image of liver
    FH2 = figure;
    imagesc(double(avgMag))
    caxis([0 80])
    hold on
    colormap(gray)
    scatter(contourCoords(:,1), contourCoords(:,2), 'r.')
    axis off
    axis image
    title('Region for stiffness estimation', 'FontSize', 12)
    saveas(FH2, sprintf('%s/region-mag.png', outDir));
    close(FH2)
    
    % Convert contour to mask
    liverMask = poly2mask(contourCoords(:,1), contourCoords(:,2), res, res);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop through subzones and calculate shear modulus for each pixel in the liver 
    
    % Put node coordinates into units of mm
    nodes(:,2) = nodes(:,2)*pixelSpacing(1);
    nodes(:,3) = nodes(:,3)*pixelSpacing(2);
    
    % Get list of x and y nodal coordinate locations
    xNodeCoords = unique(nodes(:,2));
    yNodeCoords = unique(nodes(:,3));
    
    % Initialise shear modulus map
    shearMap = zeros(res, res);
    
    % Total number of pixels to calculate shear modulus (only used for waitbar)
    tmp = liverMask(:);
    tmp(tmp==0) = [];
    totalPix = length(tmp);
    
    % Create waitbar
    WB = waitbar(0, 'Image analysis...');
    count = 0;
    
    for xx = (subZoneDim/2 + 1):(length(xNodeCoords)-subZoneDim/2)
        for yy = (subZoneDim/2 + 1):(length(yNodeCoords)-subZoneDim/2)
            
            if liverMask(xx, yy) % If central pixel of subzone is not zero
                
                count = count + 1;
                
                % Update waitbar to give user an indication of time
                percentComplete = count/totalPix;
                waitbar(percentComplete, WB, sprintf('%.2f%% complete...', 100*percentComplete))
                
                % Define range of coordinates in each direction to create subzones
                xRange = [xNodeCoords(xx-subZoneDim/2) xNodeCoords(xx+subZoneDim/2)];
                yRange = [yNodeCoords(yy-subZoneDim/2) yNodeCoords(yy+subZoneDim/2)];
                zRange = [0 5]; % Hard coded because image thickness is = 5 mm; Mesh was also made to be 5 mm thick
                
                % Creating subzone node and element lists
                [nodesSubZone, elemSubZone] = getSubZone(nodes, elems, xRange, yRange, zRange, 1, 1, 1);
                
                % Calculate the numeric virtual field
                disp('Calculating the numeric virtual field...')
                [uVF, eta, strain1] = numericVF_Iso(U, nodesSubZone, elemSubZone, nodes(:,1), [], elemType, GaussPoints);
                
                % Calculate shear modulus
                disp('Calculating the shear modulus...')
                
                if strcmp(elemType, 'C3D8R')
                    % Uniform strain elements - C3D8R
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8R(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1));
                    
                elseif strcmp(elemType, 'C3D8')
                    % Selectively reduced integration type element - C3D8
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1));
                    
                else
                    % Fully integrated element
                    [fk, fg, b, strain2] = solveIsoVFM_C3D8F(U, uVF, rho, omega, nodesSubZone, elemSubZone, nodes(:,1), GaussPoints);
                    
                end
                
                % Calculate complex shear modulus
                shearMap(xx,yy) = b/fg;
                
            end
            
        end
    end
    
    close(WB);
    
    % Save results in .mat file
    save(sprintf('%s/shear_subZoneSize%d.mat', outDir, subZoneDim), 'shearMap', 'liverMask');
    
    % Plot results
    FH = figure;
    imagesc(abs(shearMap))
    axis image
    axis off
    cbh = colorbar;
    caxis([0 8])
    xlabel(cbh, 'kPa', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    
elseif strcmp(area, 'hom')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate shear stiffness for liver region only - no subzones
    
    % Manually segment a region within the liver for VFM stiffness estimation
    disp('Creating liver contour...'); 
    liverContourFile = sprintf('%s/liverRegionContour.txt', outDir); % Filename
    contourCoords = get_contour(liverContourFile, avgMag, 'Outline the region of the liver'); % If file already exists, read in contour. Else, ask user to contour region.
    
    % Plot contour on magnitude image of liver
    FH2 = figure;
    imagesc(double(avgMag))
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
    imagesc(real(fha))
    caxis([-100 100])
    hold on
    scatter(contourCoords(:,1), contourCoords(:,2), 'ko', 'filled')
    axis off
    axis image
    title('Real(FHA)', 'FontSize', 16)
    
    subplot(1,2,2)
    imagesc(imag(fha))
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
    shear = b/fg;
    magShear = abs(shear)
    
    % Calculate normalised sensitivity value = eta/G
    nG = calcNormSensitivity(eta, shear);
       
    % Save results
    save(sprintf('%s/shearResult.mat', outDir), 'shear', 'uVF', 'eta', 'fk', 'fg', 'b', 'strain1', 'strain2', 'nG');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot measured displacement field and virtual displacement fields
    plotVF(uVF, U, liverNodes, outDir);
       
    
end
