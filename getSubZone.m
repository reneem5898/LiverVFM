function [nodesSubZone, elemSubZone] = getSubZone(nodes, elems, xRange, yRange, zRange, m, n, o)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function returns the nodes and elements in the given subzone region
%
% Inputs: 1) nodes: nodal coordinates
%         2) elems: element connectivity
%         3) xRange: list of values with define subzone regions in x dir
%         4) yRange: list of values with define subzone regions in y dir
%         5) zRange: list of values with define subzone regions in z dir
%         6) m: counter for subzone in x direction
%         7) n: counter for subzone in y direction
%         6) o: counter for subzone in z direction
%
% Outputs: 1) nodesSubZone: node coordinates in subzone region
%          2) elemSubZone: element connectivity for nodes in subzone region
%
% Written by: Renee Miller
% Updated: 3 January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodes per element
nodesPerElem = size(elems,2) - 1;

% Initialise local subzone node and element lists
nodesSubZone = zeros(20000, size(nodes,2)); % 20000 is arbitrary
elemSubZone = zeros(20000, size(elems,2));

countNode = 0;
% Get nodes in subzone
for i = 1:length(nodes)
    
    % Current coordinate
    x = nodes(i,2);
    y = nodes(i,3);
    z = nodes(i,4);
    
    if(x <= xRange(m,2) && x >= xRange(m,1) && y <= yRange(n,2) && y >= yRange(n,1) && z <= zRange(o,2) && z >= zRange(o,1)) % If coordinate is within current subzone range
        countNode = countNode + 1;
        nodesSubZone(countNode,:) = nodes(i,:);
    end
end

% Get just list of nodes in current subzone
nodesSubZone = nodesSubZone(any(nodesSubZone,2),:);
nodeSubNums = nodesSubZone(:,1);

countElem = 0;
% Get elements in subzone
for i = 1:size(elems,1)
    check = 0;
    
    % Count the number of nodes in each element which are
    % inside the subzone
    for j = 1:nodesPerElem
        if any(nodeSubNums==elems(i,1+j))
            check = check + 1;
        end
    end
    
    % If all 8 nodes are inside the subzone, save element to
    % subzone list
    if check == nodesPerElem
        countElem = countElem + 1;
        elemSubZone(countElem,:) = elems(i,:);
    end
    
end

% Get just list of elems in current subzone - get rid of zeros
elemSubZone = elemSubZone(any(elemSubZone,2),:);

% Get list of nodes which are in subzone elements
elemNodes = elemSubZone(:,2:end);
elemNodes = unique(elemNodes);

% Find nodes in nodeSubZone which are not in elemNodes and remove them
C = setdiff(nodeSubNums, elemNodes);

% Loop through C and replace with zeros
for i = 1:length(C)
    idx = find(nodeSubNums == C(i));
    nodesSubZone(idx,:) = zeros(1,4);
end

% Remove zeros from nodesSubZone
nodesSubZone(~any(nodesSubZone,2),:) = [];