function makeWaveGifMS(phaseMat, outDir)
% ------------------------------------------------------------------------
% This script creates a gif of the liver wave iamges
%
% Written by: Renee Miller
% Date: 13 March 2017
%
% ------------------------------------------------------------------------

% Gif filename
gif_filename = sprintf('%s/liver-wave.gif', outDir);

% Colormap
cmap = awave(4096);

for s = 1:size(phaseMat,1)
    
    FH = figure(1);  
    imagesc(squeeze(phaseMat(s,1,:,:)))
    %colormap(cmap)
    cbh = colorbar;
    axis image
    axis off
    xlabel(cbh, 'FHA');
    
    % save frame to gif
    fx = getframe(FH);
    imx = frame2im(fx);
    [imind,cm] = rgb2ind(imx,256);
    if s == 1
        imwrite(imind,cm,gif_filename,'gif','Loopcount',inf) %if frame #1, write
    else
        imwrite(imind,cm,gif_filename,'gif','WriteMode','append') %if > frame #1, append
    end
end

close all;