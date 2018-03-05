function [nG] = calcNormSensitivity(eta, shear)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate normalised sensitivity value for isotropic shear modulus
% (relative sensitivity to noise)
%
% Renee Miller (rmil520@aucklanduni.ac.nz)
% 1 November 2017
%
% See Avril 2004 or my thesis for documentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sensitivity: n = sqrt ( Q' H Q )
n = sqrt(shear.^2 .* eta(1));

% Normalise sensitivity value
nG = n ./ shear; % Note - nG is real only but has trailing 0.0000i
nG = real(nG); % Remove trailing 0.0000i