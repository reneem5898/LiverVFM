function [Acl, rhs] = createBoundaryConstraint(boundaryNodes, nodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the boundary constraint matrices for the numeric 
% virtual displacement field from the list of boundary nodes
%
% Inputs: 1) boundaryNodes - list of node numbers on boundaries
%         2) nodes - nodal coordinates - 4 x # nodes
% 
% Outputs: 1) Acl - LHS matrix of zeros and ones - # boundary nodes x
%             (# nodes x # nodal DOF) - ones on boundary node indices
%          2) rhs - rhs matrix of complex ones
%
% Renee Miller
% 1 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOF
DOF = size(nodes,2) - 1;

% Set values of Acl at boundary nodes to one - find indices
c = 0; % Set counter

row_ind = zeros(length(boundaryNodes)*DOF,1);
col_ind = zeros(length(boundaryNodes)*DOF,1);
val = ones(length(boundaryNodes)*DOF,1);

for i = 1:length(boundaryNodes) % For node list of boundary nodes
    
    % Get local node idx of current node within the subzone
    localNodeNum = find(nodes(:,1)==boundaryNodes(i));
    
    for d = 1:DOF % For each DOF
        
        c = c + 1; % Row counter
        row_ind(c) = c;
        
        % Index of current DOF - Column counter
        col_ind(c) = (localNodeNum*DOF)-(DOF-d);
              
    end
end

% Put together sparse matrix for H
Acl = sparse(row_ind, col_ind, val);

% RHS condition to Acl
rhs = complex(zeros(size(Acl,1),1)*(1 + 1i));