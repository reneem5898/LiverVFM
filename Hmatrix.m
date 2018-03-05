function h = Hmatrix(B, detJ, delX, delY, delZ)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the element h matrix (Connesson et al. 2015)
%
% Inputs: 1) B matrix (strain matrix)
%         2) detJ - element volume (determinant of Jacobian matrix)
%		  3) delX, delY and delZ - lengths of the sides of the element
%
% Output: h matrix
%
% Renee Miller
% 1 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate H matrix for current element
% V(Shear) ~ E(fg^2(e,eVF)) = (1/2)*stdU^2*UeVF'*H*UeVF

% t = 1/3 * Tr(e*)
t = 1/3 * sum(B(1:3,:),1);

% Calculate each component of h
h1 = (2/(delX^2))*(B(1,:) - t)'*(B(1,:) - t);
h2 = (2/(delY^2))*(B(2,:) - t)'*(B(2,:) - t);
h3 = (2/(delZ^2))*(B(3,:) - t)'*(B(3,:) - t);
h12 = (8/delX^2 + 8/delY^2)*(0.5*B(4,:))'*(0.5*B(4,:));
h13 = (8/delX^2 + 8/delZ^2)*(0.5*B(5,:))'*(0.5*B(5,:));
h23 = (8/delY^2 + 8/delZ^2)*(0.5*B(6,:))'*(0.5*B(6,:));

h = (detJ^2)*(h1 + h2 + h3 + h12 + h13 + h23); %Hg