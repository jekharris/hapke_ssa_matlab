% author: Jennifer Harris 

% A function to calculate reflectance using the Hapke function and a
% nonisotropic scattering function p(g) /= 1, p(g) == legrende polynomial

% INPUT

% w = the SSA
% X = the wavelength vector (only in here to make the hapkeRvariables minimisation work, should just cancel out)
% inc = the incidence angle (in degrees)
% emi = the emission angle (in degrees)
% g = the phase angle (in degrees)

% CALCULATIONS

% need to ensure that I feed the vectors (w,X) in as row vectors or at
% least as all the same orientation vectors, row or column

function rc = hapke_reflectanceLegrendeP(w,X)

% define angle parameters
inc = 30;
emi = 0;
g = 30;

mu = cosd(emi);
mu0 = cosd(inc);
mug = cosd(g);

% define phase function parameters
b = -0.4;
c = 0.25;

r0 = (1-sqrt(1-w))./(1+sqrt(1-w));

H_mu = 1./ ( 1 - w.*mu .*(r0 + ( (1-2*r0.*mu)/2) * log((1+mu)/mu) ) );
H_mu0 = 1./ ( 1 - w.*mu0 .*(r0 + ( (1-2*r0.*mu0)/2) * log((1+mu0)/mu0) ) );

% The Legendre polynomial p(g) as used in Sklute et al 2015 https://doi.org/10.2138/am-2015-4824

p_g1 = 1 + b.*mug + c .* (1.5 * (mug.^2) - 0.5);

rc =  ((w./(4*(mu+mu0))) .* (p_g1 -1 + (H_mu0.*H_mu)))+ ((X./X)-1);

return;