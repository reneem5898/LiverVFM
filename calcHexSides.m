function [dX, dY, dZ] = calcHexSides(X, Y, Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function which returns the mean element side lengths in x, y and z
% directions
%
% Written by: Renee Miller
% Date: 5 August 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


coords = [X Y Z];
connectivity = [1 5; 2 6; 3 7; 4 8; 1 2; 3 4; 5 6; 7 8; 1 4; 2 3; 5 8; 6 7]; % Hard coded - node connectivity for Abaqus hexahedral elements

lenSides = zeros(size(connectivity,1),1);
for i = 1:length(lenSides)
    coord1 = coords(connectivity(i,1),:); % Coordinate 1
    coord2 = coords(connectivity(i,2),:); % Coordinate 2
    lenSides(i) = distance(coord1, coord2); % Get length between coordinates
end

% Get mean side lengths in each direction
dX = mean(lenSides(1:4));
dY = mean(lenSides(5:8));
dZ = mean(lenSides(9:12));