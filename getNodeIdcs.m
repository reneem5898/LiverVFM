function nodeIdcs = getNodeIdcs(nodeList, elemsList, elemNum, DOF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns indices of nodes in nodeList
%
% Inputs: 1) nodeList - list of node numbers (either global nodes or in subregion)
%         2) elems - element connectivity
%         3) elemNum - current element number
%         4) DOF - nodal degrees of freedom
%         5) nodesPerElem - number of nodes per element (type of element)
%
% Outputs: List of nodal indices
%
% Renee Miller
% 10 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of nodes per element
nodesPerElem = size(elemsList,2) - 1;

% Index list
indexList = linspace(1,length(nodeList),length(nodeList));

% Initialise list of node indices
nodeIdcs = zeros(1,DOF*nodesPerElem);

count = 0;

% For node in current element
for n = 1:nodesPerElem
    
    % Index of node in node list
    idx = indexList(nodeList==elemsList(elemNum,n+1)); 
    
    for d = 1:DOF % For each DOF
        count = count + 1;
        
        nodeIdcs(count) = (idx*DOF)-(DOF-d); % Append index to list of indices
        
    end
end