function [FK, FG, b, strain] = solveIsoVFM_C3D8R(u, uVF, rho, omega, nodesSubZone, elemSubZone, globalNodeNums)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the RHS (b) and LHS (fk and fg) of equation 7 from
% Pierron et al. 2013
%
% Inputs: 1) u - MRE displacement field
%         2) uVF - virtual displacement field
%		  3) rho - density
%         4) omega - angular frequency
%         5) nodeCoords - list of nodal coordinates (numNodes x 3)
%         6) elems - (numElems x 8) - 8 node numbers that make up the
%         element
%         7) globalNodeNums - node numbers of all nodes in entire model
%
% Outputs: 1) fk & fg (LHS) and 2) b (RHS)
%
% Written by: Renee Miller
% Date: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Weighting for "single gauss point" - just used for constraints
w = 2;

% Flanagan 1981 - Uniform Strain Paper
% Components of shape functions
Sigma = [1, 1, 1, 1, 1, 1, 1, 1];
Lambda1 = [-1, 1, 1, -1, -1, 1, 1, -1];
Lambda2 = [-1, -1, 1, 1, -1, -1, 1, 1];
Lambda3 = [-1, -1, -1, -1, 1, 1, 1, 1];
Gamma1 = [1, 1, -1, -1, -1, -1, 1, 1];
Gamma2 = [1, -1 -1, 1, -1, 1, 1, -1];
Gamma3 = [1, -1, 1, -1 1, -1, 1, -1];
Gamma4 = [-1, 1, -1, 1, 1, -1, 1, -1];

% Initialise matrix to save material parameters
FK = 0;
FG = 0;
b = 0;

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1),7); % 7 = 1 element label + 6 strain components (e11, e22, e33, e12, e23, e13)

count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an iclar allndication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Get vectors of nodal coordinates
    [X, Y, Z] = getElemCoords(nodesSubZone, elemSubZone, i);
    
    %% Indexing - there are two types of node indexing - local and global. 
    % Local refers to the node index in the subzone and global refers to 
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
    
    % Get displacements for particular element using indices
    Ue = u(globalNodeIdcs);
    UeVF = uVF(localNodeIdcs);
    
    %% Uniform strain method - Flanagan 1981
    % This element type was implemented since it is the same as C3D8R element
    % type in Abaqus
    
    % Make B matrix - Flanagan 1981 - Appendix I
    B = makeBmatrix(X, Y, Z);
    
    % Calculate element volume using uniform strain formulation - Flanagan 1981
    detJ = calcVolumeUniformStrain(X, Y, Z);
    
    % Divide B matrix by element volume to get B matrix for element
    B = B/detJ;
        
    % Calculate strain: B*U
    eV = B*Ue; % Strain of measured displacements
    eVF_V = B*UeVF; % Strain of virtual displacement field
    
    % Save strain from measured displacement field -
    % compare to Abaqus (CHECK)
    count = count + 1;
    strain(count,:) = [i eV.'];
    
    % Convert strains to square matrices
    e = [eV(1) 0.5*eV(4) 0.5*eV(5); ...
        0.5*eV(4) eV(2) 0.5*eV(6);...
        0.5*eV(5) 0.5*eV(6) eV(3)];
    eVF =[eVF_V(1) 0.5*eVF_V(4) 0.5*eVF_V(5); ...
        0.5*eVF_V(4) eVF_V(2) 0.5*eVF_V(6);...
        0.5*eVF_V(5) 0.5*eVF_V(6) eVF_V(3)];
    
    % Calculate portion of constitutive law governed by bulk
    % modulus
    fk = trace(e)*trace(eVF)*detJ;
    
    % Calculate portion of constitutive law which contributes to shear
    % motion [2nd part should be zero since trace(eVF) = 0]
    % Note: trace(A*B') = double dot product of two matrices
    fg = 2*(trace(e*eVF.') - (1/3)*trace(e)*trace(eVF))*detJ;
    
    %%%%%%% RHS %%%%%%%
    
    % Single integration point - only for determining the volume of the element
    m = 0; n = 0; o = 0;
    
    % Shape functions - Flanagan 1981
    ShapeFuns = 1/8 * Sigma + 1/4 * m * Lambda1 + 1/4 * n * Lambda2 + 1/4 * o * Lambda3 + ...
        1/2 * n * o * Gamma1 + 1/2 * m * o * Gamma2 + 1/2 * m * n * Gamma3 + ...
        1/2 * m * n * o * Gamma4;
    
    % Compile N
    %  _                                                       _
    % | N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 0 |
    % | 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 0 |
    % | 0 0 N1 0 0 N2 0 0 N3 0 0 N4 0 0 N5 0 0 N6 0 0 N7 0 0 N8 |
    %  -                                                       -
    N = []; % Initialise the list of shape functions with a column of zeros
    for sf = 1:length(ShapeFuns) % Loop through each individual shape function
        N = [N eye(DOF).*ShapeFuns(sf)]; % Append identity shape functions (3 x 3)
    end    
    
    % RHS = Integral(rho * omega^2* u * uVF * dV)
    sh = N'*N;
    sh = diag(sum(sh,1)); % Lumped mass matrix
    acc = UeVF.'*rho*(omega^2)*sh*Ue*detJ;
    
    % Sum RHS and LHS for all elements
    FK = FK + fk*w; %LHS - bulk component
    FG = FG + fg*w; %LHS - shear component
    b = b + acc*w; %RHS - acceleration/inertial term
    
end

% Close wait bar
close(WH);
