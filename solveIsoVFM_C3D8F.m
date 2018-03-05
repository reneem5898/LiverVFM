function [FK, FG, b, strain] = solveIsoVFM_C3D8F(U, uVF, rho, omega, nodesSubZone, elemSubZone, globalNodeNums, GaussPoints)

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
%         8) GaussPoints - number of Gauss points to use (per direction)
%
% Outputs: 1) fk & fg (LHS) 
%          2) b (RHS)
%          3) strain - to validate against out put from Abaqus
%
% Written by: Renee Miller
% Date: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Zeta = integration points
% w = weights
switch GaussPoints
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
strain = zeros(size(elemSubZone,1)*size(zeta,1)*size(zeta,1)*size(zeta,1),7); % 7 = 1 element label + 6 strain components (e11, e22, e33, e12, e23, e13)

count = 0; % make counter to assign strain values

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    % Update waitbar to give user an iclar allndication of time
    percentComplete = i/size(elemSubZone,1);
    waitbar(percentComplete, WH, sprintf('%.2f%% complete...', 100*percentComplete))

    % Initialise variables and matrices
    A1_1 = 0; A2_1 = 0; B_1 = 0;
    A1_2 = 0; A2_2 = 0; B_2 = 0;
    A1_3 = 0; A2_3 = 0; B_3 = 0; 
    
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
    Ue = U(globalNodeIdcs);
    UeVF = uVF(localNodeIdcs);
    
    for j = 1:length(zeta)
        o = zeta(j);
        
        for k = 1:length(zeta)
            n = zeta(k);
            
            for l = 1:length(zeta)
                m = zeta(l);
                
                % Gauss integration coordinate = (m,n,o) %
                
                % Shape functions
                N1 = 1/8*(1-m)*(1-n)*(1-o);
                N2 = 1/8*(1+m)*(1-n)*(1-o);
                N3 = 1/8*(1+m)*(1+n)*(1-o);
                N4 = 1/8*(1-m)*(1+n)*(1-o);
                N5 = 1/8*(1-m)*(1-n)*(1+o);
                N6 = 1/8*(1+m)*(1-n)*(1+o);
                N7 = 1/8*(1+m)*(1+n)*(1+o);
                N8 = 1/8*(1-m)*(1+n)*(1+o);
                
                ShapeFuns = [N1 N2 N3 N4 N5 N6 N7 N8];
                
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
                
                % Evaluate the derivative of the shape functions at m, n, o
                DN = 0.125*[-1*(1-n)*(1-o) (1-n)*(1-o) (1+n)*(1-o) -1*(1+n)*(1-o) -1*(1-n)*(1+o) (1-n)*(1+o) (1+n)*(1+o) -1*(1+n)*(1+o);...
                    -1*(1-m)*(1-o) -1*(1+m)*(1-o) (1+m)*(1-o) (1-m)*(1-o) -1*(1-m)*(1+o) -1*(1+m)*(1+o) (1+m)*(1+o) (1-m)*(1+o);
                    -1*(1-m)*(1-n) -1*(1+m)*(1-n) -1*(1+m)*(1+n) -1*(1-m)*(1+n) (1-m)*(1-n) (1+m)*(1-n) (1+m)*(1+n) (1-m)*(1+n)];
                
                % Calculate Jacobian matrix
                jac = DN*[X Y Z];
                
                % Determinant of Jacobian
                detJ = det(jac);
                
                % Check for negative jacobians (distorted elements)
                if inv(jac) <= 0
                    %disp('Error: negative Jacobian in element');
                    elemNegJac = [elemNegJac i];
                end
                
                %%%%%%%% LHS matrix %%%%%%%%
                
                % Multiply inverse of Jacobian times the derivative of shape functions
                dNdXYZ = jac\DN;
                
                % Calculate B matrix (strain matrix)
                B = [];
                for c = 1:nodesPerElem %Loop through number of nodes per element
                    Bi = [dNdXYZ(1,c)     0              0; ...
                          0            dNdXYZ(2,c)       0; ...
                          0               0         dNdXYZ(3,c); ...
                          dNdXYZ(2,c)  dNdXYZ(1,c)       0; ...
                          dNdXYZ(3,c)     0         dNdXYZ(1,c); ...
                          0            dNdXYZ(3,c)  dNdXYZ(2,c)];
                    B = [B Bi];
                end
                
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
                
                % RHS = Integral(rho * omega^2* u * uVF * dV)
                sh = N'*N;
                sh = diag(sum(sh,1)); % Lumped mass matrix
                acc = UeVF.'*rho*(omega^2)*sh*Ue*detJ;
                
                % Sum the weighted functions
                A1_1 = A1_1 + w(l).*fk;
                A2_1 = A2_1 + w(l).*fg;
                B_1 = B_1 + w(l).*acc;
                
            end
            
            % Sum the weighted functions
            A1_2 = A1_2 + w(k).*A1_1;
            A2_2 = A2_2 + w(k).*A2_1;
            B_2 = B_2 + w(k).*B_1;
            
            % Clear the inner sums
            A1_1 = A1_1.*0;
            A2_1 = A2_1.*0;
            B_1 = B_1.*0;
            
        end
        
        % Sum the weighted functions
        A1_3 = A1_3 + w(j).*A1_2;
        A2_3 = A2_3 + w(j).*A2_2;
        B_3 = B_3 + w(j).*B_2;
        
        % Clear the inner sums
        A1_2 = A1_2.*0;
        A2_2 = A2_2.*0;
        B_2 = B_2.*0;
        
    end
    
    % Final element LHS and RHS
    A1 = A1_3;
    A2 = A2_3;
    B = B_3;
    
    % Sum RHS and LHS for all elements
    FK = FK + A1; %LHS - bulk component
    FG = FG + A2; %LHS - shear component
    b = b + B; %RHS - acceleration/inertial term
    
end

% Get number of elements with negative jacobians
countNegJac = length(unique(elemNegJac));
disp(sprintf('There were %d elements with negative Jacobians.', countNegJac));

% Close wait bar
close(WH);
