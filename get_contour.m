function [contour_coord] = get_contour(filename, im, titleStr)

% Author: Hoi Ieng (Helen) Lam
% Date created: 20090218
% Last modified on: 

% This function gets the contour coordinates of an object in an image by
% reading the specified filename if the file exists, otherwise it will ask
% user to select the contour on the image and write the selected points
% into the file.
% =========================================================================

if (exist(filename,'file')),
    contour_coord = dlmread(filename,' ');
else
    % get contour by user input
    h = msgbox(titleStr);
    waitfor(h);
    [contour_coord] = select_contour(im,titleStr);
    dlmwrite(filename,contour_coord,' ');
    close(gcf);
end;

return;