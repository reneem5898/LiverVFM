function x = distance(coord1, coord2)

% Calculate distance between two 3D points
x = sqrt((coord1(1) - coord2(1))^2 + (coord1(2) - coord2(2))^2 + (coord1(3) - coord2(3))^2);