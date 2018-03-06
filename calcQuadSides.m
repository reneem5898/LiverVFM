function [dX, dY] = calcQuadSides(X, Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function which returns the mean element side lengths in x, y directions
%
% Written by: Renee Miller
% Date: 28 February 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


coords = [X Y];
connectivity = [1 2; 2 3; 3 4; 4 1]; % Hard coded - node connectivity for Abaqus quadrilateral elements

% Initialise variables
d = zeros(size(connectivity,1),size(coords,2));

for i = 1:length(d)
    
    coord1 = coords(connectivity(i,1),:); % Coordinate 1
    coord2 = coords(connectivity(i,2),:); % Coordinate 2
    
    % Get difference between coordinates - changed 05/03/2018 RM
    d(i,:) = abs(coord1 - coord2);
end

% Get delX and delY - changed 05/03/2018 RM
tmp = max(d);
dX = tmp(1);
dY = tmp(2);