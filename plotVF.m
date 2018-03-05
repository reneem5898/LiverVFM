function plotVF(uVF, U, nodes, outDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function plots the measured first harmonic amplitude for the region  
% of interest in the liver along with the the virtual displacement field for 
% the same region
%
% Inputs: 1. uVF - virtual displacement field
%         2. U - First harmonic amplitude (measured displacement)
%         3. nodes - nodal coordinates in region of interest
%         4. outDir - directory where to save fig
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Date: 1 November 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just get Y component
uVF_y = uVF(2:3:end);
u_y = U(2:3:end);

% Return measured y displacements in region of interest
u_y = u_y(nodes(:,1));

% Get limits for plots
uVF_Rlim = max(abs(real(uVF_y)))*0.8;
uVF_Ilim = max(abs(imag(uVF_y)))*0.8;
u_Rlim = max(abs(real(u_y)))*0.5;
u_Ilim = max(abs(imag(u_y)))*0.5;

% marker size
bubbleSZ = 12;

% Get indices of nodes in middle slice
indices = find(nodes(:,4)==3);

% Plot real virtual displacements (scatter), where color = value
FH = figure('Position',[50, 50, 1200, 900]);
subplot(2,2,2)
scatter3(nodes(indices,2),nodes(indices,3),nodes(indices,4), bubbleSZ, real(uVF_y(indices)), 'filled')
view([0 90])
cbh = colorbar;
caxis([-1*uVF_Rlim uVF_Rlim])
xlabel(cbh, 'FHA', 'FontSize', 12)
set(cbh, 'FontSize', 12)
axis equal tight
title('U_y^{*R}', 'FontSize', 16)
axis off

% Plot real reference displacements (scatter), where color = value
subplot(2,2,1)
scatter3(nodes(indices,2),nodes(indices,3),nodes(indices,4), bubbleSZ, real(u_y(indices)), 'filled')
view([0 90])
cbh = colorbar;
caxis([-1*u_Rlim u_Rlim])
xlabel(cbh, 'FHA', 'FontSize', 12)
set(cbh, 'FontSize', 12)
axis equal tight
title('U_y^{R}', 'FontSize', 16)
axis off

% Plot imaginary virtual displacements (scatter), where color = value
subplot(2,2,4)
scatter3(nodes(indices,2),nodes(indices,3),nodes(indices,4), bubbleSZ, imag(uVF_y(indices)), 'filled')
view([0 90])
cbh = colorbar;
caxis([-1*uVF_Ilim uVF_Ilim])
xlabel(cbh, 'FHA', 'FontSize', 12)
set(cbh, 'FontSize', 12)
axis equal tight
title('U_y^{*I}', 'FontSize', 16)
axis off

% Plot imaginary reference displacements (scatter), where color = value
subplot(2,2,3)
scatter3(nodes(indices,2),nodes(indices,3),nodes(indices,4), bubbleSZ, imag(u_y(indices)), 'filled')
view([0 90])
cbh = colorbar;
caxis([-1*u_Ilim u_Ilim])
xlabel(cbh, 'FHA', 'FontSize', 12)
set(cbh, 'FontSize', 12)
axis equal tight
title('U_y^{I}', 'FontSize', 16)
axis off

% Save image in directory
saveas(FH, sprintf('%s/uVF_Y_complex.png', outDir))

% Close figure
close(FH);