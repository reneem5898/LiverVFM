% Run all liver noise analyses - signal to noise ratios

close all
clear all

% Parent directory - where all data is located
parentDir = 'P:/Data/Liver';
AllFolders = dir(parentDir);
allFolders = {AllFolders.name};
allFolders(ismember(allFolders,{'.','..'})) = [];

% Scanning sequences
sequences = {'SIEMENS', '923', '887', '923*fract'};
%sequences = {'923*fract'};

% Counter
c = 0;

% Create variables to store all mean and standard deviation of pdsnr values
% in region of interest
mean_PDSNR_vals = zeros(6, length(sequences));
stdev_PDSNR_vals = zeros(6, length(sequences));
all_PDSNR_vec = []; 

% Loop through all folders
for i = 1:length(allFolders)
    
    if any(strfind(allFolders{i},'MRE_')) && ~any(strfind(allFolders{i}, '_F')) % Use only folders that have MRE in the name and not F (failed)
        
        % Volunteer directory
        volunteerDir = sprintf('%s/%s', parentDir, allFolders{i});
        disp(volunteerDir);
        
        % Counter
        c = c + 1;
        
        % Loop through all sequences
        for j = 1:length(sequences)
            
            % Scanning sequence
            seq = sequences{j};
            disp(seq)
            
            % Magnitude images
            file_spec = sprintf('%s/*%s*Mag', volunteerDir, seq);
            tmp = dir(file_spec);
            tmpName = {tmp.name};
            magDir = sprintf('%s/%s', volunteerDir, tmpName{1});
            
            % Output directory
            tmp = strsplit(magDir, '/');
            tmp2 = strsplit(tmp{end}, '_');
            tmp3 = strsplit(seq, '*');
            if length(tmp3) > 1
                outDir = sprintf('%s/%s_%s_results', volunteerDir, tmp2{1}, tmp3{end});
            else
                outDir = sprintf('%s/%s_results', volunteerDir, tmp2{1});
            end
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            disp(outDir);
            
            % Load data needed for noise analysis
            load(strcat(outDir, '/displacements.mat')); % Load data from previous analysis
            load(strcat(outDir, '/regionMesh_2D.mat')); % Load mat file with region contour
            
            pdsnrFile = sprintf('%s/pdsnr.mat', outDir);
            if ~exist(pdsnrFile, 'file')
                % Calculate phase difference SNR
                disp('Calculating PDSNR weights...');
                [pdsnr, regionMask] = PDSNR_Weights_Liver(fha, avgMag, imageContourCoords); % Inputs = first harmonic amplitudes, average magnitude images and contour coordinates
                save(pdsnrFile,'pdsnr', 'regionMask');
            else
                load(pdsnrFile);
            end
            
            % Plot weights: output vector of PDSNR values, mean and
            % standard deviation in ROI
            [pdsnrVec, mean_PDSNR_vals(c, j), stdev_PDSNR_vals(c,j)] = plot_PDSNR_Liver(pdsnr, regionMask, outDir);
            
            % Save vector of PDSNR values
            seqVals = j * ones(length(pdsnrVec),1);
            all_PDSNR_vec = [all_PDSNR_vec; seqVals, pdsnrVec];
            
        end
    end
end

% Save combined results to mat file
save(sprintf('%s/pdsnr-all-results.mat', parentDir), 'mean_PDSNR_vals', 'stdev_PDSNR_vals', 'all_PDSNR_vec');

% Plot mean PDSNR values for each sequence type
[X,~] = meshgrid(1:length(sequences), 1:size(mean_PDSNR_vals,1));
FH = figure;
scatter(X(1,:), mean_PDSNR_vals(1,:), 60, 'bo', 'filled')
hold on
scatter(X(2,:), mean_PDSNR_vals(2,:), 60, 'ro', 'filled')
scatter(X(3,:), mean_PDSNR_vals(3,:), 60, 'mo', 'filled')
scatter(X(4,:), mean_PDSNR_vals(4,:), 60, 'yo', 'filled')
scatter(X(5,:), mean_PDSNR_vals(5,:), 60, 'go', 'filled')
scatter(X(6,:), mean_PDSNR_vals(6,:), 60, 'ko', 'filled')
xticklabels(sequences)
xticks([1 2 3 4])
xlim([0.5 4.5])
xlabel('Sequence', 'FontSize', 12)
ylabel('Mean PDSNR', 'FontSize', 12)
legend('1', '2', '3', '4', '5', '6', 'Location', 'Best')
saveas(FH, sprintf('%s/PDSNR-mean-volunteer.png', parentDir))

% Sort vector values by first column
allPDSNR = sortrows(all_PDSNR_vec,1);

% Calculate means for all PDSNR values for each sequence
mean_SIEMENS = mean(allPDSNR(allPDSNR(:,1)==1,2));
mean_923 = mean(allPDSNR(allPDSNR(:,1)==2,2));
mean_887 = mean(allPDSNR(allPDSNR(:,1)==3,2));
mean_923F = mean(allPDSNR(allPDSNR(:,1)==4,2));
mean_Sequences = [mean_SIEMENS, mean_923, mean_887, mean_923F];

% Calculate standard deviations for all PDSNR values for each sequence
std_SIEMENS = std(allPDSNR(allPDSNR(:,1)==1,2));
std_923 = std(allPDSNR(allPDSNR(:,1)==2,2));
std_887 = std(allPDSNR(allPDSNR(:,1)==3,2));
std_923F = std(allPDSNR(allPDSNR(:,1)==4,2));
std_Sequences = [std_SIEMENS, std_923, std_887, std_923F];

% Plot boxplots of all PDSNR values
FH2 = figure;
%boxplot(allPDSNR(:,1), allPDSNR(:,2))
scatter([1, 2, 3, 4], mean_Sequences, 50, 'bo', 'filled')
hold on
errorbar(mean_Sequences, std_Sequences, 'b', 'LineWidth', 1)
% N = length(mean_Sequences);
% errorbar(reshape([reshape(mean_Sequences',2,[]);nan(1,N)],1,[]), ...
%          reshape([reshape(std_Sequences',2,[]);nan(1,N)],1,[]))
xticklabels(sequences)
xticks([1 2 3 4])
xlim([0.5 4.5])
%ylim([-0.05 7])
xlabel('Sequence', 'FontSize', 12)
ylabel('PDSNR', 'FontSize', 12)
saveas(FH2, sprintf('%s/PDSNR-mean-std.png', parentDir))

