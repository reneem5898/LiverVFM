% Run all liver vfm analyses back to back

close all
clear all

addpath('C:\Users\rmil520\Documents\MATLAB\MRI_class');
import MRI.*

% Parent directory - where all data is located
parentDir = 'P:/Data/Liver';
AllFolders = dir(parentDir);
allFolders = {AllFolders.name};
allFolders(ismember(allFolders,{'.','..'})) = [];

% Scanning sequences
sequences = {'SIEMENS', '923', '887', '923*fract'};

% Number of slices
numSlices = 4;

% Counter
c = 0;

% Initialise results variables
magG = zeros(6, length(sequences), numSlices); % Abs value
G = zeros(6, length(sequences), numSlices); % Complex shear modulus
nG = zeros(6, length(sequences), numSlices); % Sensitivity value

% Loop through all folders
for i = 1:length(allFolders)
    
    if any(strfind(allFolders{i},'MRE_')) && ~any(strfind(allFolders{i}, '_F')) % Use only folders that have MRE in the name and not F (failed)
        
        % Counter
        c = c + 1;
        
        % Volunteer directory
        volunteerDir = sprintf('%s/%s', parentDir, allFolders{i});
        disp(volunteerDir);
        
        % Loop through all sequences
        for j = 1:length(sequences)
            
            % Scanning sequence
            seq = sequences{j};
            disp(seq)
            
            % Run analysis
%             liverVFM_multiSlice(volunteerDir, seq); %3D analysis
%             liverVFM_multiSlice_2D(volunteerDir, seq); %2D analysis
            
            % Collect results
            [mg, g, ng] = getResultsLiverVFM(volunteerDir, seq, numSlices);
            magG(c, j, :) = mg;
            G(c, j, :) = g;
            nG(c, j, :) = ng;
            
        end
    end
end

% volunteerDir = 'P:/Data/Liver/MRE_20171214_2';
% % Loop through all sequences
% for j = 1:length(sequences)
%
%     % Scanning sequence
%     seq = sequences{j};
%
%     % Run analysis
%     liverVFM_multiSlice(volunteerDir, seq);
% end

% Save matfile with all results
save(sprintf('%s/results-allLiver-2D.mat', parentDir), 'magG', 'G', 'nG')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results

% Plot all shear stiffness values and sensitivity values
figure(1)
scatter(abs(nG(:)), magG(:), 'ro', 'filled')
ylabel('Shear Modulus (kPa)', 'FontSize', 16)
xlabel('Normalised Sensitivity (\eta/G)', 'FontSize', 16)
set(gca, 'FontSize', 12)

% Plot shear stiffness values and sensitivities - by sequence
nG_product = squeeze(nG(:,1,:));
nG_923 = squeeze(nG(:,2,:));
nG_887 = squeeze(nG(:,3,:));
nG_923F = squeeze(nG(:,4,:));

magG_product = squeeze(magG(:,1,:));
magG_923 = squeeze(magG(:,2,:));
magG_887 = squeeze(magG(:,3,:));
magG_923F = squeeze(magG(:,4,:));

figure(2)
scatter(abs(nG_product(:)), magG_product(:), 'ro', 'filled')
hold on
scatter(abs(nG_923(:)), magG_923(:), 'bo', 'filled')
scatter(abs(nG_887(:)), magG_887(:), 'ko', 'filled')
scatter(abs(nG_923F(:)), magG_923F(:), 'mo', 'filled')
ylabel('Shear Modulus (kPa)', 'FontSize', 16)
xlabel('Normalised Sensitivity (\eta/G)', 'FontSize', 16)
set(gca, 'FontSize', 12)
legend('Product', '923', '887', '923F', 'Location', 'Best')
