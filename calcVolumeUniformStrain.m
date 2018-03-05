function detJ = calcVolumeUniformStrain(X, Y, Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the element volume using the uniform strain
% formulation - Flanagan 1981
%
% Inputs: X, Y and Z nodal coordinate vectors
%
% Output: detJ - volume of element
%
% Renee Miller
% 28 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IP% Integration point
IP = [0, 0, 0];
m = IP(1); n = IP(2); o = IP(3);

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

% Dertivative of shape functions - Flanagan 1981
% dN1 = 1/4 * Lambda1 + 1/2 * o * Gamma2 + 1/2 * n * Gamma3 + 1/2 * n * o * Gamma4; %% Abaqus documentation has factor of 1/2
% dN2 = 1/4 * Lambda2 + 1/2 * o * Gamma1 + 1/2 * m * Gamma3 + 1/2 * m * o * Gamma4;
% dN3 = 1/4 * Lambda3 + 1/2 * n * Gamma1 + 1/2 * m * Gamma2 + 1/2 * m * n * Gamma4;
dN1 = 1/4 * Lambda1 + 1/2 * o * Gamma2 + 1/2 * n * Gamma3 + n * o * Gamma4; %% Original paper has no factor of 1/2
dN2 = 1/4 * Lambda2 + 1/2 * o * Gamma1 + 1/2 * m * Gamma3 + m * o * Gamma4;
dN3 = 1/4 * Lambda3 + 1/2 * n * Gamma1 + 1/2 * m * Gamma2 + m * n * Gamma4;

dN = [dN1; dN2; dN3];

% Calculate Jacobian for current element
jac = dN*[X Y Z];

% Element volume
detJ = det(jac);