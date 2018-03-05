function [contour] = select_contour(im, titleStr)

% Author: Hoi Ieng (Helen) Lam
% Date created: 20090218
% Last modified by: Hoi Ieng (Helen) Lam
% Last modified on: 

% This function is based on the code in 'get_dti_contour.m'
% This function obtains the contours of an object in an image by interactive user input. 

% Input:
% im: any image

% Outputs:
% contour: coordinates of the selected contour 


% Plot the image
figure('Name',titleStr);
imagesc(im); 
caxis([0 80]) % Colorscale for Siemens magnitude image of abdomen (liver)
colormap(gray);
title({'Left mouse button picks points.', 'Right mouse button picks last point of the outer contour.'});

% Ask the user whether the contour is visible for defining or not,
% default answer is 'Yes'
answer = questdlg('Is the object visible for contour selection?','Contour Selection','Yes','No','Yes');

if strcmp('No', answer) == 1,
    disp('No contour selected');
    contour = [];
    close(gcf);
    return;
end;

% Initialise string expression
answer1 = 'No';
p = []; % stores the handles for the points and spline of the outer contour

while (strcmp('No', answer1) == 1),
    
    % Display help dialogs for outer contour selection
    h1 = helpdlg('Select points from the figure to define the contour of an object.','Contour Selection');
    waitfor(h1);
    h2 = helpdlg({'Left mouse button picks points.', 'Right mouse button picks last point of the outer contour.'},'Outer Contour Selection');
    waitfor(h2);
    
    % Remove existing points and spline
    delete(p); p = []; 
    
    coord = [];
    button = 1;  % button = 1 for left button, = 3 for right button
    num_pt = 1; % Initialize the number of points 
    % Let user to select the outer contour
    while (button==1),

        % Get the coordinate (xb1, yb1) of the selected point
        [xb1,yb1,button] = ginput(1);

        % Save the coordinates of the outer contour selected by the
        % user
        coord(num_pt,:) = [xb1,yb1];
        
        if button==3, % The last point selected for outer contour
            hold on;
            p(num_pt) = plot(xb1,yb1,'rx','MarkerSize',10); 
            continue; % Terminate the selection process
        else
            % Plot the selected points in red on the same figure
            hold on;
            p(num_pt) = plot(xb1, yb1,'rx','MarkerSize',10);
        end;
        num_pt = num_pt + 1;
    end;

    % Set the last point to be the same as the first selected point, so as 
    % to complete a loop
    coord(num_pt + 1,:) = coord(1,:);

    % Generate spline for the contour and then samples the spline
    % over a finer mesh (ts)
    t = 1:num_pt+1;
    ts = 1:0.1:num_pt+1;
    xy = coord'; % convert into two rows
    xys = spline(t, xy, ts);
    hold on;
    % Plot the spline onto the image
    p(num_pt+1) = plot(xys(1,:), xys(2,:), '-r');
    hold on;

    % Ask user if the contour is fitted satisfactorily, default 
    % answer is 'yes'
    answer1 = questdlg('Is the contour satisfactory?','Contour','Yes','No','Yes');

    % Notify the user that contour selection has been completed
    if strcmp('Yes', answer1) == 1,
        h3 = helpdlg('Contour selection completed.', 'Contour Selection');
        waitfor(h3);
        contour = xys'; % the sampled points on the spine
    end;

end;
  
return;
