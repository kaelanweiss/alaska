function m = solveFM17(U,u,w,Tw,Ti,L)
% Function to solve the iceberg side melt parameterization of FitzMaurice
% et al. 2017.
%
% m = solve3Eqn(U,u,w,Tw,Ti,L)
%
% Input
%   U: ocean velocity magnitude (m/s)
%   u: horizontal ocean velocity (m/s)
%   w: vertical ocean velocity (m/s)
%   Tw: ocean temperature (C)
%   Ti: iceberg temperature (C)
%   L: iceberg scale
%
% Output
%   m: melt rate (m/day)
%
% KJW
% 7 Apr 2024

% constants
K = 0.58;

% melt rate (detached plume)
m = K*(U.^0.8).*(Tw-Ti)./(L.^0.2);

% melt rate (attached plume)
alpha = U./(sqrt(2)*w);
Tp = alpha.*Tw;
m_a = K*(w.^0.8).*(Tp-Ti)./(L.^0.2);

% combine according to u vs w
idx_a = w>abs(u);
m(idx_a) = m_a(idx_a);