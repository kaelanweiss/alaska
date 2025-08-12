% Script to parse and process the downward-looking ADCP data on the
% BlueBoat.
%
% KJW
% 5 Aug 2024

clear

%%% USER INPUT %%%
% specify data location
% data_path = '../data/';
data_path = 'F:\bear_202409\data\raw\adcp';
dep_name = '102620_bear_20240908';
save_output = 1;

%%% USER INPUT %%%
% specify time window (if desired)
t1 = datetime(2024,9,8,21,25,0);
t2 = datetime(2024,9,8,22,12,0);

% parse ADCP MAT files output from MIDAS
adcp_raw = parseNortekDeployment(fullfile(data_path,dep_name),t1,t2);

% plot raw data
ax1 = adcpQuickPlot(figure(1),adcp_raw,'vel',0.75*[-1 1],NaT,[0 100],1);
text(ax1(1),0.01,1,'beam velocity, raw','units','normalized','verticalalignment','bottom')
ax2 = adcpQuickPlot(figure(2),adcp_raw,'amp',[30 85],NaT,[0 100],1);
text(ax2(1),0.01,1,'beam amplitude','units','normalized','verticalalignment','bottom')
ax3 = adcpQuickPlot(figure(3),adcp_raw,'cor',[40 100],NaT,[0 100],1);
text(ax3(1),0.01,1,'beam correlation','units','normalized','verticalalignment','bottom')

%% QC
adcp = adcp_raw;

%%% USER INPUT %%%
% correlation limit
cor_min = 80;
adcp.burst.vel(adcp.burst.cor<cor_min) = nan;

%%% USER INPUT %%%
% bt fom limit
fom_max = 1e3;
adcp.bt.vel(adcp.bt.fom>fom_max) = nan;

% plot QC data
ax4 = adcpQuickPlot(figure(4),adcp,'vel',0.75*[-1 1],NaT,[0 100],1);
text(ax4(1),0.01,1,'beam velocity, post QC','units','normalized','verticalalignment','bottom')

% BT
figure(5); clf
flds = {'vel','distance','fom'};
for i = 1:3
    ax5(i) = subplot(3,1,i);
    plot(adcp.bt.time,adcp.bt.(flds{i}))
    ylabel(['BT ' flds{i}])
end
linkaxes(ax5,'x')

%% Remove Bottom Track Velocities
% make sure bottom track aligns in time with velocity
if length(adcp.bt.time) ~= length(adcp.burst.time)
    bt_vel_interp = nan(length(adcp.burst.time),adcp.bt.nbeams);
    for i = 1:adcp.bt.nbeams
        bt_vel_interp(:,i) = interp1(adcp.bt.time,adcp.bt.vel(:,i),adcp.burst.time,'nearest');
    end
    adcp.bt.vel = bt_vel_interp;
end

% clean up BT velocities using a hanning filter
%%% USER INPUT %%%
% specify hanning filter width
bt_hann_width = 0.5; % sec, you'll want to keep this pretty low compared to the sample rate

adcp.bt.vel_filt = adcp.bt.vel;
for i = 1:size(adcp.bt.vel,2)
    adcp.bt.vel_filt(:,i) = hannFilter(adcp.bt.vel(:,i),bt_hann_width*adcp.burst.samplerate,'omitnan');
end

% calculate 5th beam BT vel
bt5 = (1/cosd(25))*mean([mean(adcp.bt.vel_filt(:,[1 3]),2) mean(adcp.bt.vel_filt(:,[2 4]),2)],2,'omitnan');
adcp.bt.vel_filt = cat(2,adcp.bt.vel_filt,bt5);

% remove BT velocities
adcp.burst.vel_bt = adcp.burst.vel - permute(repmat(adcp.bt.vel_filt,[1 1 adcp.burst.ncells]),[1 3 2]);

% plot velocity with BT removed
ax6 = adcpQuickPlot(figure(6),adcp,'vel_bt',0.5*[-1 1],NaT,[0 100],1);
text(ax6(1),0.01,1,'beam velocity, BT removed','units','normalized','verticalalignment','bottom')

%% Coordinate Transformations
% beam to XYZ (instrument)
R_xyz = zeros(5);
% R_xyz(1:4,1:4) = adcp.burst.beam2xyz;

% rearrange so that output is [x y z(b5) z(b1234) error]'
R_xyz(1:2,1:4) = adcp.burst.beam2xyz(1:2,1:4); % x, y
R_xyz(3,5) = 1; % z from beam 5
R_xyz(4,1:4) = sum(adcp.burst.beam2xyz(3:4,:),1)/2; % z from beams 1-4
R_xyz(5,1:4) = adcp.burst.beam2xyz(3,:) - adcp.burst.beam2xyz(4,:); % error vel
adcp.burst.vel_xyz_lbls = {'x','y','z (b5)','z (b1234)','error'};

% XYZ to ENU
%%% USER INPUT %%%
% magnetic declination correction (set to zero if set in MIDAS)
% +14.7 deg at Bear Glacier
mag_dec = 14.7;

R_mag = eye(5);
R_mag(1:2,1:2) = [cosd(mag_dec) sind(mag_dec); -sind(mag_dec) cosd(mag_dec)];

% use AHRS rotation matrix for u, v, w(b5), w(b1234), and leave 5 (error) alone
R_enu = reshape(adcp.attitude.ahrsrotationmatrix,[size(adcp.attitude.ahrsrotationmatrix,1) 3 3]);

% preallocate transformed velocity and rearrange dimensions to make linear 
% algebra easier (beam, time, cell)
adcp.burst.vel_xyz = nan*permute(adcp.burst.vel,[3 1 2]);
adcp.burst.vel_enu = adcp.burst.vel_xyz;

% calculate transforms (loop through cells)
for i = 1:size(adcp.burst.vel,2)
    % beam to XYZ
    adcp.burst.vel_xyz(:,:,i) = R_xyz*squeeze(adcp.burst.vel_bt(:,i,:))';

    % XYZ to ENU
    R1 = eye(5);
    R1(1:3,1:3) = R_enu(i,:,:); % AHRS rotation
    R1(4,:) = [R_enu(i,3,1) R_enu(i,3,2) 0 R_enu(i,3,3) 0];
    R2 = R_mag*R1; % magnetic declination
    adcp.burst.vel_enu(:,:,i) = R2*adcp.burst.vel_xyz(:,:,i);
end
adcp.burst.vel_enu_lbls = {'u','v','w (b5)','w (b1234)','error'};

% permute velocities back
adcp.burst.vel_xyz = permute(adcp.burst.vel_xyz,[2 3 1]);
adcp.burst.vel_enu = permute(adcp.burst.vel_enu,[2 3 1]);

% plot transformed velocity
ax7 = adcpQuickPlot(figure(7),adcp,'vel_enu',1*[-1 1],NaT,[0 100],1);
text(ax7(1),0.01,1,'ENU velocity','units','normalized','verticalalignment','bottom')
for i = 1:5
    text(ax7(i),0.01,1,adcp.burst.vel_enu_lbls{i},'units','normalized','verticalalignment','top')
end

% save output
if save_output
    save(fullfile(data_path,dep_name,'adcp.mat'),'adcp','-v7.3')
end