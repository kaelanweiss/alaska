function [m,Tb,Sb] = solve3EqnSchulz22(u,T,S,z,u_thresh)
% Function to solve the 3 equation melt parameterization of Jenkins et al.
% (2010) with updated transfer velocities suggested by Schulz et al. 2022.
%
% [m,Tb,Sb] = solve3Eqn(u,Tw,Sw,z)
%
% Input
%   u: ocean velocity (m/s)
%   T: ocean temperature (C)
%   S: ocean salinity (psu)
%   z: depth (m) (positive downward, z=0 at surface)
%   u_thresh: threshold velocity (m/s)
%
% Output
%   m: melt rate (m/s)
%   Tb: boundary (freezing) temperature
%   Sb: boundary salinity
%
% KJW
% 7 Apr 2024

% physical parameters
rho_i = 916; % ice density (kg/m3)
Lf = 334000; % latent heat of fusion (water) (J/kg)
rho_w = 1030; % water density (kg/m3)
cw = 3974; % water specific heat (J/C kg)
lambda1 = -0.0573; % liquidus slope (C/psu)
lambda2 = 0.0832; % liquidus intercept (C)
lambda3 = -7.53e-8; % liquidus pressure coefficient (C/Pa)

% model parameters
Cd = 0.0097;
GammaT = 0.011;
GammaS = 3.1e-4;
gammaT_star = 1e-3; % m/s

% size of inputs
sz = size(u);

% turn depth into 
if numel(z)==1
    z = z*ones(sz);
end

% calculate absolute pressure
P = 101300 + rho_w*9.81*z;

% calculate gammaT, gammaS
idx = u>u_thresh;
gammaT = gammaT_star*ones(sz);
gammaT(idx) = gammaT(idx) + GammaT*sqrt(Cd)*(u(idx)-u_thresh);
gammaS = .07*gammaT;

% define coefficients
a1 = rho_i*Lf;
a2 = rho_w*cw*gammaT;
a3 = a2.*T;
b1 = lambda2 + lambda3*P;
c1 = rho_i;
c3 = rho_w*gammaS;
c2 = c3.*S;


% solve boundary salinity
A = -a2*lambda1;
B = a3 - a2.*b1 + a1.*c3./c1;
C = -a1.*c2./c1;
 
Sb = (-B + sqrt(B.^2 - 4.*A.*C))./(2*A);

% boundary temperature
Tb = b1 + lambda1*Sb;

% melt rate
m = (a3 - a2.*Tb)/a1;






