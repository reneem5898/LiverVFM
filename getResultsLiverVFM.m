function [magG, G, sensitivity] = getResultsLiverVFM(volunteerDir, sequence, numSlices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get results from VFM isotropic stiffness for liver data
%
% Written by: Renee Miller
% 20 February 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define results directory
if strcmp(sequence, '923')
    file_spec = sprintf('%s/*%s_results', volunteerDir, sequence);
else
    file_spec = sprintf('%s/*%s*results', volunteerDir, sequence);
end
tmp = dir(file_spec);
tmpName = {tmp.name};
outDir = sprintf('%s/%s', volunteerDir, tmpName{1});
disp(outDir)

% Initialise result variables
magG = zeros(numSlices,1);
G = zeros(numSlices,1);
sensitivity = zeros(numSlices,1);

% Load results for each slice
for i = 1:numSlices
   
    load(sprintf('%s/shearResult_slice%d_2D.mat', outDir, i));
    
    magG(i) = magShear;
    G(i) = shear;
    sensitivity(i) = nG;
    
end