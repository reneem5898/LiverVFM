function nodeIdcs = getNodeIdcs(nodeList, elems, elemNum, DOF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns indices of nodes in nodeList
%
% Inputs: 1) nodeList - list of node numbers (either global nodes or in subregion)
%         2) elems - element connectivity
%         3) elemNum - current element number
%         4) DOF - nodal degrees of freedom
%
% Outputs: List of nodal indices
%
% Renee Miller
% 10 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elems, 2) - 1;

% Index list
indexList = linspace(1,length(nodeList),length(nodeList));

% Initialise list of node indices
nodeIdcs = [];

% For node in current element
for n = 1:nodesPerElem
    
    % Index of node in node list
    idx = indexList(nodeList==elems(elemNum,n+1)); 
    
    for d = 1:DOF % For each DOF
        nodeIdcs = [nodeIdcs (idx*DOF)-(DOF-d)]; % Append index to list of indices
        
    end
end