function plotLiverPhaseOffsets(phaseMat, fha, outDir, x, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the fitted complex harmonic from one (arbitrarily chosen)
% pixel in the liver MRE data to double check fit
%
% Renee Miller
% 12 October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of phase offsets
noff = size(phaseMat,1);

% Create x axis for phase offsets
xOff = linspace(1,noff,noff);

% Pick one pixel to check
tmpFHA = phaseMat(:,x,y);
complexFHA = fha(x,y);

% Remove dc component from tmpFHA
tmpFHA = tmpFHA - mean(tmpFHA);

% Discretise complex fit
disc = 400; %400 is arbitrary - bigger number = smoother curve
phasePlot = zeros(1,disc); 
xPlot = linspace(1,noff,disc);
for i = 1:disc
    k = xPlot(i);
    phasePlot(:,i) = abs(complexFHA).*cos(2*pi*(k-1)/noff+angle(complexFHA));
end

% Plot
FH = figure;
scatter(xOff, tmpFHA, 'ro', 'filled')
hold on
plot(xPlot, phasePlot, 'r-')
legend('Discretised Displacement', 'Fitted Harmonic', 'Location', 'Best')
xlabel('Arbitrary Time Scale', 'FontSize', 16)
ylabel('First Harmonic Amplitude', 'FontSize', 16)
saveas(FH, sprintf('%s/phase-fit.png', outDir))