function boundaryNodes = getBoundaryNodes(nodes, elems)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the nodes, elements and # of nodes per element, this function
% returns the nodes which are on the bounary of the model or subregion
%
% Inputs: 1) nodes - 4 x # nodes
%         2) elems - 9 x # elems
%
% Output: boundaryNodes - node list
%
% Renee Miller
% 01 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of nodes per element
nodesPerElem = size(elems,2) - 1;

%Initialise list of all nodes which comprise all elements
allNodes = sort(elems(:,2:end));
allNodes = allNodes(:);

% Variable for counting instances of each node
countNodes = zeros(size(nodes,1),1);

% Loop through nodes and count number of times that it appears in the list of element nodes
for i = 1:size(nodes,1)
    tmp = find(allNodes == nodes(i,1));
    countNodes(i) = length(tmp); % Count how many times each node is found 
end

boundaryNodes = nodes(:,1); % Start out with list of all node numbers
boundaryNodes(countNodes==nodesPerElem) = 0; % If count == 8 [meaning - if node is in 8 different elements, it is an internal node], remove from list
boundaryNodes(boundaryNodes==0) = [];
