function [FK, FG, b, strain] = solveIsoVFM_C3D8(u, uVF, rho, omega, nodesSubZone, elemSubZone, globalNodeNums)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the RHS (b) and LHS (fk and fg) of equation 7 from
% Pierron et al. 2013
%
% Inputs: 1) u - MRE displacement field
%         2) uVF - virtual displacement field
%		  3) rho - density
%         4) omega - angular frequency
%         5) nodeCoords - list of nodal coordinates (numNodes x 3)
%         6) elems - (numElems x 9) - index + 8 node numbers that make up the
%         element
%         7) globalNodeNums - node numbers of all nodes in entire model
%         8) elemType - element type: 'C3D8R', 'C3D8', 'C3D8F'
%         9) GaussPoints - number of Gauss points per direction (1, 2, 3)
%
% Outputs: 1) fk & fg (LHS) 
%          2) b (RHS)
%          3) strain - to validate with Abaqus
%
% Written by: Renee Miller
% Date: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nodal DOFs
DOF = size(nodesSubZone,2) - 1;

% Gauss point coordinates
zeta = [-sqrt(1/3); sqrt(1/3)];
% Integration point weights
w = [1; 1];
% Number of Gauss points per element
GP = 8;

% Initialise variables
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
strain = zeros(size(elemSubZone,1)*size(zeta,1)*size(zeta,1)*size(zeta,1),8); % 7 = 1 element label + gauss point + 6 strain components (e11, e22, e33, e12, e23, e13)

% Loop through each element
for i = 1:size(elemSubZone,1)
    
    count = 0; % make counter to assign strain values
    
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
    Ue = u(globalNodeIdcs);
    UeVF = uVF(localNodeIdcs);
    
    % Make B matrix - Flanagan 1981 - Appendix I
    Br = makeBmatrix(X, Y, Z);
       
    % Calculate element volume using uniform strain formulation - Flanagan 1981
    elementVolume = calcVolumeUniformStrain(X, Y, Z);
    
    % Divide B matrix by element volume to get B matrix for element
    % Reduced integration used for volumetric part of strain 
    Br = Br/elementVolume;
    
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
                
                % Calculate full integration B matrix (strain matrix)
                Bf = [];
                for c = 1:nodesPerElem %Loop through number of nodes per element
                    Bi = [dNdXYZ(1,c)     0              0; ...
                          0            dNdXYZ(2,c)       0; ...
                          0               0         dNdXYZ(3,c); ...
                          dNdXYZ(2,c)  dNdXYZ(1,c)       0; ...
                          dNdXYZ(3,c)     0         dNdXYZ(1,c); ...
                          0            dNdXYZ(3,c)  dNdXYZ(2,c)];
                    Bf = [Bf Bi];
                end
                
                % Combine rows 1-3 of Br with rows 4-6 of Bf to get 
                % "selectively reduced" B matrix
                %B = [Br(1:3,:); Bf(4:6,:)];
                
                % Calculate strain: B*U
                eV_r = Br*Ue; % Strain of measured displacements
                eV_f = Bf*Ue; % Fully integrated strain
                eVF_V_r = Br*UeVF; % Strain of virtual displacement field
                eVF_V_f = Bf*UeVF; % Fully integrated virtual strain
                
                % Save strain from measured displacement field - compare to Abaqus (CHECK)
                count = count + 1;
                idx = (i-1)*GP + count;
                strain(idx,:) = [i count eV_r(1:3).' eV_f(4:6).'];
                
                % Convert strains to square strain tensor
                eF = [eV_f(1) 0.5*eV_f(4) 0.5*eV_f(5); ...
                    0.5*eV_f(4) eV_f(2) 0.5*eV_f(6);...
                    0.5*eV_f(5) 0.5*eV_f(6) eV_f(3)];
                
                eR = [eV_r(1) 0.5*eV_r(4) 0.5*eV_r(5); ...
                    0.5*eV_r(4) eV_r(2) 0.5*eV_r(6);...
                    0.5*eV_r(5) 0.5*eV_r(6) eV_r(3)];
                
                eVF_f =[eVF_V_f(1) 0.5*eVF_V_f(4) 0.5*eVF_V_f(5); ...
                    0.5*eVF_V_f(4) eVF_V_f(2) 0.5*eVF_V_f(6);...
                    0.5*eVF_V_f(5) 0.5*eVF_V_f(6) eVF_V_f(3)];
                
                eVF_r =[eVF_V_r(1) 0.5*eVF_V_r(4) 0.5*eVF_V_r(5); ...
                    0.5*eVF_V_r(4) eVF_V_r(2) 0.5*eVF_V_r(6);...
                    0.5*eVF_V_r(5) 0.5*eVF_V_r(6) eVF_V_r(3)];
                                                            
                % Calculate portion of constitutive law governed by bulk
                % modulus
                fk = trace(eR)*trace(eVF_r)*detJ;
                
                % Calculate portion of constitutive law which contributes to shear
                % motion [2nd part should be zero since trace(eVF) = 0]
                % Note: trace(A*B') = double dot product of two matrices
                fg = 2*(trace(eF*eVF_f.') - (1/3)*trace(eF)*trace(eVF_f))*detJ; 
                
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
