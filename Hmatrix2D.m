function h = Hmatrix2D(B, detJ, delX, delY)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the element h matrix (Connesson et al. 2015)
%
% Inputs: 1) B matrix (strain matrix)
%         2) detJ - element volume (determinant of Jacobian matrix)
%		  3) delX, delY - lengths of the sides of the element
%
% Output: h matrix
%
% Renee Miller (rmil520@auckland.ac.nz)
% 2 March 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate H matrix for current element
% V(Shear) ~ E(fg^2(e,eVF)) = (1/2)*stdU^2*UeVF'*H*UeVF

% Calculate each component of h
h1 = (2/(delX^2))*B(1,:)'*B(1,:);
h2 = (2/(delY^2))*B(2,:)'*B(2,:);
h12 = (8/delX^2 + 8/delY^2)*(0.5*B(3,:))'*(0.5*B(3,:));

h = (detJ^2)*(h1 + h2 + h12); %Hg