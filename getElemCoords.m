function [X, Y, Z] = getElemCoords(nodes, elems, elemNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns vectors of X, Y and Z nodal coordinates for the nodes
% which make up the current element
%
% Inputs: 1) nodes - nodal coordinates (4 x # elements)
%         2) elems - element connectivity
%         3) elemNum - current element number
%
% Outputs: Vectors of X, Y, and Z coordinates (size depends on element type
% - number of nodes per element)
%
% Renee Miller
% 10 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of nodes per element
nodesPerElem = size(elems,2) - 1;

% Initialise vectors
X = zeros(nodesPerElem,1);
Y = zeros(nodesPerElem,1);
Z = zeros(nodesPerElem,1);

% Loop through nodes in current element
for j = 1:nodesPerElem
    
    % Get index of node in node list which corresponds to current node in current element
    idx = find(nodes(:,1)==elems(elemNum,j+1)); 
    
    % Save coordinates
    X(j) = nodes(idx,2);
    Y(j) = nodes(idx,3);
    Z(j) = nodes(idx,4);
end