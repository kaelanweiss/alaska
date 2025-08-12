function adv = msADVTransform(adv,varargin)
% Function to perform ADV transform for meltstake Vector. QC is performed,
% and velocities are transformed into instrument and ice coordinates. Tilt 
% is taken into account through the adcp structure output from the 
% "parseNortekDeployment" function. The coordinate system is consistent
% with those reported in the 2024 and 2025 GRL papers (i.e. v points away
% from the ice).
%
% adv = msADVTransform(adv)
% adv = msADVTransform(adv,attitude)
% adv = msADVTransform(adv,cor_min,snr_min)
% adv = msADVTransform(adv,attitude,cor_min,snr_min)
%
% Input
%   adv: adv structure as output by "parseADVDeployment" function
%   attitude: (optional) structure with fields time, pitch, and roll in 
%       degrees; pitch and roll need to actually correspond to those 
%       measurements, this means R = -pitch and P = roll-270 if using ADCP 
%       output
%   cor_min: (optional) correlation minimum for QC, default 90
%   snr_min: (optional) amplitude minimum for QC, default 0
%
% Output
%   adv: adv structure with new fields adv.vel_xyz (instrument
%       coordinates), adv.vel_ice (attitude corrected), adv.processing
%
% version control:
% 31 Jul 2025 - made script aware of "corrected" ADCP tilt
%
% KJW
% 6 Mar 2024

% parse varargin
cor_min = 90;
snr_min = 0;
attitude = struct('time',extrema(adv.time),'pitch',[0 0],'roll',[0 0]);
tilt_corrected = 0;
switch nargin
    case 2
        attitude = varargin{1};
        tilt_corrected = 1;
    case 3
        cor_min = varargin{1};
        snr_min = varargin{2};
    case 4
        attitude = varargin{1};
        cor_min = varargin{2};
        snr_min = varargin{3};
        tilt_corrected = 1;
end

% QA
qa_vel = adv.snr<snr_min | adv.cor<cor_min;
vel = adv.vel;
vel(qa_vel) = nan;

% beam 2 instrument (consistent with positive velocity toward receivers and
% coordinate system on pg 39 of the comprehensive manual NOT PG 21!!)
adv.vel_xyz = (adv.beam2xyz*vel')';

% beam to ice (pre-attitude correction)
xyz2ice = [ 0  1  0;...
            0  0 -1;...
           -1  0  0];
adv.vel_ice = (xyz2ice*adv.vel_xyz')';
ice_lbls = {'u (left)','v (away from ice)','w (up)'};

% attitude correction
if tilt_corrected
    % get pitch and roll, correct if necessary, interpolate
    if isfield(attitude,'pitch_corrected')
        roll = interp1angle(attitude.time,attitude.roll_corrected,adv.time);
        pitch = interp1angle(attitude.time,attitude.pitch_corrected,adv.time);
    else
        roll = -interp1angle(attitude.time,attitude.pitch*pi/180,adv.time); % when +y_inst points down, "pitch" is zero and increases with negative rotation about the +z_inst axis
        pitch = interp1angle(attitude.time,(attitude.roll-270)*pi/180,adv.time); % "roll" is positive rotation about +x_inst, at 270 when +y_inst points down
    end

    % recreate attitude struct for output
    attitude = struct('time',adv.time,'pitch',pitch,'roll',roll);
    
    % rotate at each time step
    for i = 1:length(adv.time)
        veli = adv.vel_ice(i,:)';
        
        % matrix elements
        cr = cos(roll(i));
        sr = sin(roll(i));
        cp = cos(pitch(i));
        sp = sin(pitch(i));
    
        % rotation matrices
        % roll is negative rotation about +y_ice
        R = [ cr   0 -sr;...
               0   1   0;...
              sr   0  cr];
    
        % pitch is negative rotation about +x_ice
        P = [  1   0   0;...
               0  cp  sp;...
               0 -sp  cp];
        
        % apply rotation
        veli = R*P*veli;
        adv.vel_ice(i,:) = veli';
    end
end

% add processing information
adv.processing = struct('cor_min',cor_min,'snr_min',snr_min,...
    'xyz2ice',xyz2ice,'ice_lbls',{ice_lbls},'tilt_corrected',tilt_corrected);
adv.processing.attitude = attitude;





