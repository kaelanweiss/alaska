function beta = adcpIceSlope(r1,r2)
% Function to calculate ice slope relative to the Sig1000 ADCP using two
% opposing beam distances. Geometry is calculated for the Sig1000, but
% could be generalized to other ADCPs by changing the beam angle and
% distance between ADCP coordinate origin and location of beam convergence.
% Positive slope results from r1>r2. Swapping inputs will result in the
% same, but opposite, angle.
%
% beta = adcpIceSlope(r1,r2)
%
% Input
%   r1: location of surface in beam 1 as measured along ADCP z-axis
%   r2: location of surface in beam 2 as measured along ADCP z-axis
%
% Output
%   beta: angle (rad) of ice relative to the ADCP
%
% KJW
% 20 Oct 2025

% ADCP geometry
beam_angle = 25; % (deg)
conv_dist = 0.163; % (m) distance between ADCP coordinate origin and location of beam convergence

% calculate slope
beta = atan((r1-r2)/((r1+r2+2*conv_dist)*tand(beam_angle)));