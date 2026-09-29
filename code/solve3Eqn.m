function [m,Tb,Sb] = solve3Eqn(u,T,S,z,varargin)
% Function to solve the 3 equation melt parameterization of Jenkins et al.
% (2010). The 3 equations form a nonlinear system due to the m*Sb term.
%
% [m,Tb,Sb] = solve3Eqn(u,Tw,Sw,z)
% [m,Tb,Sb] = solve3Eqn(u,Tw,Sw,z,params)
%
% Input
%   u: ocean velocity (m/s)
%   T: ocean temperature (C)
%   S: ocean salinity (psu)
%   z: depth (m) (positive downward, z=0 at surface)
%   params: (optional) vector or structure storing model parameters
%       vector: [Cd, GammaT, GammaS]
%       struct: params.fld with flds Cd, GammaT, and GammaS
%
% Output
%   m: melt rate (m/s)
%   Tb: boundary (freezing) temperature
%   Sb: boundary salinity
%
% KJW
% 7 Feb 2024

% default model parameters
Cd = 0.0097;
GammaT = 0.011;
GammaS = 3.1e-4;
defaultParams = [Cd, GammaT, GammaS];

% parse inputs
p = inputParser;
validParams = @(x) (isvector(x) && length(x)==3) || isstruct(x);
addOptional(p,'params',defaultParams,validParams);
parse(p,varargin{:});
if isstruct(p.Results.params)
    Cd = p.Results.params.Cd;
    GammaT = p.Results.params.GammaT;
    GammaS = p.Results.params.GammaS;
else
    Cd = p.Results.params(1);
    GammaT = p.Results.params(2);
    GammaS = p.Results.params(3);
end

% physical parameters
rho_i = 916; % ice density (kg/m3)
Lf = 334000; % latent heat of fusion (water) (J/kg)
rho_w = 1030; % water density (kg/m3)
cw = 3974; % water specific heat (J/C kg)
lambda1 = -0.0573; % liquidus slope (C/psu)
lambda2 = 0.0832; % liquidus intercept (C)
lambda3 = -7.53e-8; % liquidus pressure coefficient (C/Pa)

% size of inputs
sz = size(u);

% turn depth into a vector
if numel(z)==1
    z = z*ones(sz);
end

% calculate absolute pressure
P = 101300 + rho_w*9.81*z;

% define coefficients
a1 = rho_i*Lf;
a2 = rho_w*cw*sqrt(Cd)*GammaT*u;
a3 = a2.*T;
b1 = lambda2 + lambda3*P;
c1 = rho_i;
c3 = rho_w*sqrt(Cd)*GammaS*u;
c2 = c3.*S;


% solve boundary salinity
A = -a2*lambda1;
B = a3 - a2.*b1 + a1.*c3./c1;
C = -a1.*c2./c1;
 
Sb = (-B + sqrt(B.^2 - 4.*A.*C))./(2*A);

% boundary temperature
Tb = b1 + lambda1*Sb;

% melt rate
m = a2.*(T - Tb)/a1;






