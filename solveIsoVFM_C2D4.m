function [FK, FG, b, strain] = solveIsoVFM_C2D4(U, uVF, rho, omega, nodesSubZone, elemSubZone, globalNodeNums, GP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the RHS (b) and LHS (fk and fg) of equation 7 from
% Pierron et al. 2013
%
% Inputs: 1) u - MRE displacement field
%         2) uVF - virtual displacement field
%		  3) rho - density
%         4) omega - angular frequency
%         5) nodesSubZone - list of nodal coordinates (numNodes x 3) in ROI
%         6) elemSubZone - (numElems x 8) - node numbers that make up the
%         elements in the ROI
%         7) globalNodeNums - node numbers of all nodes in entire model
%         8) GaussPoints - number of Gauss points to use (per direction)
%
% Outputs: 1) fk & fg (LHS)
%          2) b (RHS)
%          3) strain - to validate against output from Abaqus
%
% Written by: Renee Miller
% Date: 6 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Zeta = integration points
% w = weights
switch GP
    case 1
        zeta = 0;
        w = 2;
    case 2
        zeta = [-sqrt(1/3); sqrt(1/3)];
        w = [1; 1];
    case 3
        zeta = [-sqrt(3/5); 0; sqrt(3/5)];
        w = [5/9; 8/9; 5/9];
end

% Initialise matrix to save material parameters
FK = 0;
FG = 0;
b = 0;

% Create waitbar
WH = waitbar(0, 'Assembling global matrices...');

% Initialise variable to create list of elemenets with negative jacobians
elemNegJac = [];

% Number of nodes per element
nodesPerElem = size(elemSubZone,2) - 1;

% Variable to save strain at Gauss points (to comopare to Abaqus output)
strain = zeros(size(elemSubZone,1)*size(zeta,1)*size(zeta,1)*size(zeta,1),5); % 7 = 1 element label + 6 strain components (e11, e22, e33, e12, e23, e13)

count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an iclar allndication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))
    
    % Initialise variables and matrices
    A1_1 = 0; A2_1 = 0; B_1 = 0;
    A1_2 = 0; A2_2 = 0; B_2 = 0;
    
    % Get vectors of nodal coordinates
    [X, Y] = getElemCoords2D(nodesSubZone, elemSubZone, i);
    
    %% Indexing - there are two types of node indexing - local and global.
    % Local refers to the node index in the subzone and global refers to
    % the node index in the whole model.
    
    % Get indics of nodal DOF in global list
    globalNodeIdcs = getNodeIdcs(globalNodeNums, elemSubZone, i, DOF);
    
    % Get indics of nodal DOF in global list
    localNodeIdcs = getNodeIdcs(nodesSubZone(:,1), elemSubZone, i, DOF);
    
    % Get displacements for particular element using indices
    Ue = U(globalNodeIdcs);
    UeVF = uVF(localNodeIdcs);

    % Loop through gauss points
    for j = 1:length(zeta)
        n = zeta(j);
        
        for k = 1:length(zeta)
            m = zeta(k);
            
            % Gauss integration coordinate = (m,n) %
            
            % Shape functions
            N1 = 1/4*(1-m)*(1-n);
            N2 = 1/4*(1+m)*(1-n);
            N3 = 1/4*(1+m)*(1+n);
            N4 = 1/4*(1-m)*(1+n);
            
            
            ShapeFuns = [N1 N2 N3 N4];
            
            % Compile N
            %  _                           _
            % | N1 0 0 N2 0 0 N3 0 0 N4 0 0 |
            % | 0 N1 0 0 N2 0 0 N3 0 0 N4 0 |
            % | 0 0 N1 0 0 N2 0 0 N3 0 0 N4 |
            %  -                           -
            N = []; % Initialise the list of shape functions with a column of zeros
            for sf = 1:length(ShapeFuns) % Loop through each individual shape function
                N = [N eye(DOF).*ShapeFuns(sf)]; % Append identity shape functions (3 x 3)
            end
            
            % Evaluate the derivative of the shape functions at m, n
            dNdm = 0.25 * [-1*(1-n) (1-n) (1+n) -1*(1+n)];
            dNdn = 0.25 * [-1*(1-m) -1*(1+m) (1+m) (1-m)];
                        
            DN = [dNdm; dNdn];
            
            % Calculate Jacobian matrix
            jac = DN*[X Y];
            
            % Determinant of Jacobian
            detJ = det(jac);
            
            % Check for negative jacobians (distorted elements)
            if inv(jac) <= 0
                %disp('Error: negative Jacobian in element');
                elemNegJac = [elemNegJac i];
            end
            
            %%%%%%%% LHS matrix %%%%%%%%
            
            % Multiply inverse of Jacobian times the derivative of shape functions
            dNdXY = jac\DN;
            
            % Calculate B matrix (strain matrix)
            B = [];
            for c = 1:nodesPerElem %Loop through number of nodes per element
                Bi = [dNdXY(1,c)     0        ; ...
                      0            dNdXY(2,c) ; ...
                      dNdXY(2,c)   dNdXY(1,c)];
                B = [B Bi];
            end
            
            % Calculate strain: B*U
            eV = B*Ue; % Strain of measured displacements
            eVF_V = B*UeVF; % Strain of virtual displacement field
            
            % Save strain from measured displacement field -
            % compare to Abaqus (CHECK)
            count = count + 1;
            strain(count,:) = [i count eV.'];
            
            % Convert strains to square matrices
            e = [eV(1) 0.5*eV(3); ...
                0.5*eV(3) eV(2) ];
            eVF = [eVF_V(1) 0.5*eVF_V(3); ...
                   0.5*eVF_V(3) eVF_V(2) ];
            
            % Calculate portion of constitutive law governed by bulk
            % modulus
            fk = trace(e)*trace(eVF)*detJ;
            
            % Calculate portion of constitutive law which contributes to shear
            % motion [2nd part should be zero since trace(eVF) = 0]
            % Note: trace(A*B') = double dot product of two matrices
            fg = 2*(trace(e*eVF.') - (1/3)*trace(e)*trace(eVF))*detJ;
            
            %%%%%%% RHS %%%%%%%
            
            % RHS = Integral(rho * omega^2* u * uVF * dV)
            sh = N'*N;
            sh = diag(sum(sh,1)); % Lumped mass matrix
            acc = UeVF.'*rho*(omega^2)*sh*Ue*detJ;
            
            % Sum the weighted functions
            A1_1 = A1_1 + w(k).*fk;
            A2_1 = A2_1 + w(k).*fg;
            B_1 = B_1 + w(k).*acc;
            
        end
        
        % Sum the weighted functions
        A1_2 = A1_2 + w(j).*A1_1;
        A2_2 = A2_2 + w(j).*A2_1;
        B_2 = B_2 + w(j).*B_1;
        
        % Clear the inner sums
        A1_1 = A1_1.*0;
        A2_1 = A2_1.*0;
        B_1 = B_1.*0;
        
    end
    
    % Sum RHS and LHS for all elements
    FK = FK + A1_2; %LHS - bulk component
    FG = FG + A2_2; %LHS - shear component
    b = b + B_2; %RHS - acceleration/inertial term
    
end

% Get number of elements with negative jacobians
countNegJac = length(unique(elemNegJac));
disp(sprintf('There were %d elements with negative Jacobians.', countNegJac));

% Close wait bar
close(WH);
