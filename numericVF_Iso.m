function [uVF, eta, strain] = numericVF_Iso(U, nodesSubZone, elemSubZone, globalNodeNums, surfaceNodes, elemType, GaussPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates numeric virtual displacement field following the
% methods outlined in Connesson et al. for a harmonic displacement field
% The following conditions are specified:
% 1. fk = 0 (Ak)
% 2. fg = 1 (Ag)
% 3. sum(uVFx) = 0, sum(uVFy) = 0, sum(uVFz) = 0 (Arb)
% 4. uVF(boundaries) = 0 (Acl)
% 5. Noise minimisation (H)
%
% Inputs: 1) u - MRE displacement field
%         2) nodesSubZone - list of node numbers and nodal coordinates (numNodes x 3)
%         3) elemsSubZone - (numElems x 8) - 8 node numbers that make up the
%         element
%         4) globalNodeNums - all nodes numbers in entire model
%         5) surfaceNodes - list of surface nodes (not required)
%         6) elemType - string denoting integration type to use: 'C3D8R', 'C3D8', 'C3D8F'
%         7) GaussPoints - number of integration points (per direction) to use
%            - input to 'C3D8F' integration type
%
% Outputs: 1) uVF - special and optimised virtual displacement field
%          2) eta - sensitivity value
%          3) strain - strain calculated from given displacements (used to 
%             compare with strains output from Abaqus - verification)
%
% Note: This version of the numericVF function calculates a virtual field 
% for a subzone in the model.
%
% Written by: Renee Miller
% Date: 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Construct Acl matrix - boundary nodes uVF = 0
% Find boundary nodes and set displacement condition to zero

disp('Constructing boundary node constraint matrix...');

% Get boundary nodes, unless given
if isempty(surfaceNodes)
    boundaryNodes = getBoundaryNodes(nodesSubZone, elemSubZone);
else
    boundaryNodes = surfaceNodes;
end

% Create boundary constraint matrices
[Acl, rhs_Acl] = createBoundaryConstraint(boundaryNodes, nodesSubZone);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Construct Ak, Ag and H matrices

disp('Constructing Ak, Ag and Hg constraint matrices...');

% Uniform strain elements - C3D8R
if strcmp(elemType, 'C3D8R')
    [Ak, Ag, H, strain] = isoC3D8R(U, nodesSubZone, elemSubZone, globalNodeNums);

% Selectively reduced integration type element - C3D8
elseif strcmp(elemType, 'C3D8')
    [Ak, Ag, H, strain] = isoC3D8(U, nodesSubZone, elemSubZone, globalNodeNums);

% Fully integrated element 
else
    [Ak, Ag, H, strain] = isoC3D8F(U, nodesSubZone, elemSubZone, globalNodeNums, GaussPoints);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Construct Arb - condition where the sum of all displacements in each
% direction: x, y and z, respectively are 0

disp('Constructing the constraint on the sum of virtual displacements in x, y and z...');

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Construct repeated identity matrix: I I I I ... N times (N = number of
% nodes)
Arb = zeros(DOF,length(nodesSubZone)*DOF);
for i = 1:size(Arb,1)
    Arb(i,i:DOF:end) = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Compile constraint matrices and solve for virtual displacement field

disp('Solving for the virtual displacement field...')

A = [Acl; Arb; Ak; Ag];
A = sparse(A);

% Compile LHS
%z1 = spalloc(size(A,1),size(A,1),0);
z1 = sparse(size(A,1),size(A,1));
LHS = [H A.'; A z1]; % Need .' rather than ' since A is complex - want just transpose, not conjugate transpose

% Compile RHS
%     Acl       %%%%%%%%% Arb %%%%%%%%%%%      Ak         Ag
Zg = [rhs_Acl; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 0 + 0*1i; 1.0 + 0.0*1i];
z2 = complex(zeros(size(H,1),1)); % Minimise H
RHS = [z2; Zg];

% Solve for virtual displacements and lagrangian multipliers which satisfy conditions
x = LHS\RHS;

% Return uVF (virtual displacement field)
uVF = x(1:size(H,1));

% Return eta (sensitivity)
eta = uVF' * H * uVF;
etaReal = real(uVF)' * real(H) * real(uVF);
etaImag = imag(uVF)' * imag(H) * imag(uVF);

eta = [eta, etaReal, etaImag];