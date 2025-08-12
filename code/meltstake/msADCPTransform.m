function adcp = msADCPTransform(adcp,varargin)
% Function to perform ADCP transform for meltstake Sig1000. QC is
% performed, and velocities are transformed into instrument and ice
% coordinates. Instrument tilt is taken into account. Requires adcp
% structure output from "parseNortekDeployment" function.
%
% adcp = msADCPTransform(adcp)
% adcp = msADCPTransform(adcp,cor_min,amp_min)
%
% Input
%   adcp: adcp structure with fields "burst" and "attitude" (as output by
%       parseNortekDeployment function)
%   cor_min: (optional) correlation minimum for QC, default 50
%   amp_min: (optional) amplitude minimum for QC, default 0
%
% Output
%   adcp: adcp structure with new subfields burst.vel_xyz (instrument
%       coordinates), burst.vel_ice (pitch and tilt correction,
%       transformation described by INST2ICE), and burst.vel_err (error
%       velocities)
%   
% version control:
% 18 Mar 2025 - coordinate system updated
% 31 Jul 2025 - added "corrected" pitch and roll fields
%
% KJW
% 30 Jan 2024

% parse varargin
cor_min = 50;
amp_min = 0;
switch nargin
    case 3
        cor_min = varargin{1};
        amp_min = varargin{2};
end

% instrument-to-ice transformation matrix
INST2ICE = [ -1   0   0   0   0;... % inst x      ice x
              0  -1   0   0   0;... % inst y      ice z 
              0   0   0   0  -1;... % inst z1 --> ice y (b5)
              0   0 -.5 -.5   0;... % inst z2     ice y (b1-4)
              0   0   1  -1   0];   % inst b5     error 
% vel component labels
ice_lbls = {'u (left)','w (up)','v b5 (ice-out)','v b1-4','err (1,3-2,4)'};
err_lbls = {'v1,3-v5','v2,4-v5','v1,3-v2,4'};

% helpful values
[nt,nc,nb] = size(adcp.burst.vel);

% velocity field (use vel_unw if available)
if isfield(adcp.burst,'vel_unw')
    vel = adcp.burst.vel_unw;
    vel_unw_flag = true;
else
    vel = adcp.burst.vel;
    vel_unw_flag = false;
end

% quality control
qa_cor = adcp.burst.cor<cor_min;
qa_amp = adcp.burst.amp<amp_min;
qa = qa_cor | qa_amp;
vel(qa) = nan;

% embed 4x4 beam-to-instrument matrix in a 5x5 matrix (if 5 beam data)
B2I = eye(nb);
B2I(1:size(adcp.burst.beam2xyz,1),1:size(adcp.burst.beam2xyz,1)) = adcp.burst.beam2xyz;

% time series of pitch and roll (attitude and burst are on offset time axes)
adcp.attitude.pitch_corrected = adcp.attitude.roll-270; % "roll" is positive rotation about +x_inst, at 270 when +y_inst points down
adcp.attitude.roll_corrected = -adcp.attitude.pitch; % when +y_inst points down, "pitch" is zero and increases with negative rotation about the +z_inst axis
roll = interp1angle(adcp.attitude.time,adcp.attitude.roll_corrected*pi/180,adcp.burst.time);
pitch = interp1angle(adcp.attitude.time,adcp.attitude.pitch_corrected*pi/180,adcp.burst.time);

% preallocate
vel_xyz = nan(nb,nc,nt);
vel_ice = nan(nb,nc,nt);

% loop through time steps
for k = 1:nt
    % matrix of beam velocities at every cell for one timestamp
    velk = squeeze(vel(k,:,:))';

    % beam to instrument transformation
    vel_xyz(:,:,k) = B2I*velk;
    
    % instrument to ice coordinates (with vertical w)
    % this is a little confusing since pitch and roll are applied to
    % opposite axes than what would usually be done due to the way the ADCP
    % reports pitch and roll in the horizontal configuration, this is
    % because the ADCP reports pitch and roll extrinsically
    CR = cos(roll(k));
    SR = sin(roll(k));
    CP = cos(pitch(k));
    SP = sin(pitch(k));

           %   (u   w  v5 v1-4  err)
    Rpitch = [  1   0   0   0   0;... % u
                0  CR  SR   0   0;... % w'
                0 -SR  CR   0   0;... % v5'
                0 -SR   0  CR   0;... % v1-4'
                0   0   0   0   1];   % err
    
            %   (u   w  v5 v1-4  err)
    Rroll = [  CP  SP   0   0   0;... % u'
              -SP  CP   0   0   0;... % w'
                0   0   1   0   0;... % v5
                0   0   0   1   0;... % v1-4
                0   0   0   0   1];   % err
    
    R = Rpitch*Rroll;
    vel_ice(:,:,k) = R*INST2ICE*vel_xyz(:,:,k);
end

% put dimensions back in order
vel_xyz = permute(vel_xyz,[3 2 1]);
vel_ice = permute(vel_ice,[3 2 1]);

% add new velocity fields and useful info to the adcp.burst structure
t_proc = datetime('now','TimeZone','UTC');
adcp.burst.vel_xyz = vel_xyz;
adcp.burst.vel_ice = vel_ice;
adcp.burst.vel_err = cat(3,vel_xyz(:,:,3)-vel_xyz(:,:,5),vel_xyz(:,:,4)-vel_xyz(:,:,5),vel_xyz(:,:,3)-vel_xyz(:,:,4));
adcp.burst.processing = struct('vel_unwrap',vel_unw_flag,'cor_min',cor_min,'amp_min',amp_min,'beam2xyz',B2I,'inst2ice',INST2ICE,'ice_lbls',{ice_lbls},'err_lbls',{err_lbls},'proc_timestamp',sprintf('%sZ',t_proc));

