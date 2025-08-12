% Script to turn spider deployment into a single usable format.
%
% KJW
% 26 Sep 2022

clear;

addpath('../../../code')

% set berg name
bergs = {'20220824_teaparty',...
         '20220824_singingflower',...
         '20220825_singingflower'};
berg = bergs{2};

fprintf('%s\n',berg);
save_processed_output = false;

% set paths
raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw';
mat_path = 'F:/Alaska2022/data/iceberg_surveys/mat';
proc_path = 'F:/Alaska2022/data/iceberg_surveys/proc';

% full paths
depname_raw = fullfile(raw_path,berg,'spider');
depname_mat = fullfile(mat_path,berg,'spider');
depname_proc = fullfile(proc_path,berg,'spider');

if ~exist(depname_raw,'dir')
    warning('raw directory: %s not found, exiting...',depname_raw);
    return
end

% make mat and proc directories if they don't exist
if ~exist(depname_mat,'dir')
    mkdir(depname_mat);
end
if ~exist(depname_proc,'dir')
    mkdir(depname_proc);
    mkdir(fullfile(depname_proc,'figures'));
end
if ~exist(fullfile(depname_proc,'figures'),'dir')
    mkdir(fullfile(depname_proc,'figures'));
end

%% turn raw data into mat files, or load
% force raw data to be reparsed
force_parse = true;

% load deployment times file
tlim = parseDeploymentTimes(depname_raw);

% RBR
fprintf('\nRBR...\n');
if exist(fullfile(depname_mat,'rbr.mat'),'file') && ~force_parse
    % load mat files
    fprintf('loading\n');
    rbr = load(fullfile(depname_mat,'rbr.mat'));
    rbr = rbr.rbr;
else
    % make mat file
    fprintf('parsing\n');
    rbr = parseRBRDeployment(depname_raw,tlim(1),tlim(2));
    fprintf('saving %s\n',fullfile(depname_mat,'rbr.mat'));
    save(fullfile(depname_mat,'rbr.mat'),'rbr');
end

% ADCP
fprintf('\nADCP...\n');
if exist(fullfile(depname_mat,'adcp.mat'),'file') && ~force_parse
    % load mat file
    fprintf('loading\n');
    adcp = load(fullfile(depname_mat,'adcp.mat'));
    adcp = adcp.adcp;
else
    % make mat file
    fprintf('parsing\n');
    adcp = parseNortekDeployment(depname_raw,tlim(1),tlim(2));
    echo = adcp.echo;
    adcp = rmfield(adcp,'echo');
    fprintf('saving %s\n',fullfile(depname_mat,'adcp.mat'));
    save(fullfile(depname_mat,'adcp.mat'),'adcp');
    fprintf('saving %s\n',fullfile(depname_mat,'echo.mat'));
    save(fullfile(depname_mat,'echo.mat'),'echo');
    clear echo
end

% load beam distance file
load(fullfile(depname_raw,'beam_distance.mat'))

% load Dragon ROV if available
use_dragon = false;
rov_file = fullfile(mat_path,berg,'dragon','rov.mat');
if exist(rov_file,'file')
    rov = load(rov_file);
    rov = rov.rov;
    use_dragon = true;
end

%% convert to distance-from-ice
fprintf('\nconverting to distance-from-ice...\n');
adcp.burst = beam2ice(adcp.burst,beam_distance);

%% QA velocity
qa = struct();
qa.corr_thresh = 50;
qa.ev_thresh = 0.25;

% correlation
for j = 1:adcp.burst.nbeams
    corr_mask = adcp.burst.cor(:,:,j)<qa.corr_thresh;
    tmp = adcp.burst.vel(:,:,j);
    tmp(corr_mask) = nan;
    adcp.burst.vel(:,:,j) = tmp;
end

%% beam to instrument coordinates
% reshape to make more efficient (loop through cells not time)
adcp.burst.vel = permute(adcp.burst.vel,[3 1 2]);
% preallocate
adcp.burst.vel_xyz = nan*adcp.burst.vel;

% loop through cells
fprintf('\ntransforming beam to instrument...\n');
for i = 1:size(adcp.burst.vel,3)
    % x
    adcp.burst.vel_xyz(1,:,i) = adcp.burst.beam2xyz(1,[1 3])*adcp.burst.vel([1 3],:,i);
    % y
    adcp.burst.vel_xyz(2,:,i) = -adcp.burst.beam2xyz(2,[2 4])*adcp.burst.vel([2 4],:,i);
    % z1
    adcp.burst.vel_xyz(3,:,i) = adcp.burst.beam2xyz(3,[1 3])*adcp.burst.vel([1 3],:,i);
    % error
    adcp.burst.vel_xyz(4,:,i) = adcp.burst.beam2xyz(4,[2 4])*adcp.burst.vel([2 4],:,i);
    adcp.burst.vel_xyz(4,:,i) = adcp.burst.vel_xyz(4,:,i) - adcp.burst.vel_xyz(3,:,i);
end

%% time bin
fprintf('\ntime binning...\n');
% reshape again to make time coordinate last
while size(adcp.burst.vel,3) ~= length(adcp.burst.time)
    adcp.burst.vel = permute(adcp.burst.vel,[3 1 2]);
    adcp.burst.vel_xyz = permute(adcp.burst.vel_xyz,[3 1 2]);
end

% calculate time bins
bin_width = 10; % seconds
adcp.burst.bin_width = bin_width;
% beam
[vel_binned,tbins] = time_bin(adcp.burst.time,adcp.burst.vel,bin_width/86400,'omitnan');
adcp.burst.tbins = tbins';
adcp.burst.vel_binned = vel_binned;
% xyz
[vel_binned,tbins] = time_bin(adcp.burst.time,adcp.burst.vel_xyz,bin_width/86400,'omitnan');
adcp.burst.vel_xyz_binned = vel_binned;

% reshape again to original (nt x nc x nb)
for tmp = {'vel','vel_xyz','vel_binned','vel_xyz_binned'}
    tmp = tmp{1};
    adcp.burst.(tmp) = permute(adcp.burst.(tmp),[3 1 2]);
end


%% %%%%%%%%% scratch %%%%%%%%%
% make attitude plot
fprintf('\nplotting...\n')
% change roll limits
roll = adcp.attitude.roll;
idx = roll<0;
roll(idx) = roll(idx)+360;
adcp.attitude.roll = roll;

% clean up heading wrap if needed
if max(abs(diff(adcp.attitude.heading))) > 350
    adcp.attitude.heading = cleanAngleWrap(adcp.attitude.heading);
end

% berg title
berg_splt = strsplit(berg,'_');
berg_ttl = sprintf('%s %s (local)',berg_splt{2},datestr(datenum(berg_splt{1},'yyyymmdd'),'dd-mmm'));

% plot
save_figs = false;
flds = {'pitch','roll','heading','pressure'};
fs = 14;
lw = 1;
figure(1); clf;
ax = [];
for i = 1:length(flds)
    subplot(length(flds),1,i)
    ax(end+1) = gca;
    plot(adcp.attitude.time,adcp.attitude.(flds{i}),'linewidth',lw)
    grid;
    xlim(tlim)
    ylim(extrema(adcp.attitude.(flds{i})))
    datetick('x','mmmdd HH:MM','keeplimits')
    ylabel(flds{i},'fontsize',fs)
    if strcmp(flds{i},'pressure')
        axis ij
        if use_dragon
            hold on
            plot(rov.time,-rov.altitudeAMSL,'k')
            ylim(extrema(-rov.altitudeAMSL))
            legend({'spider','dragon'},'location','southeast','fontsize',fs-2);
        end
    end
    if i == 1
        title(sprintf('spider attitude - %s',berg_ttl),'fontsize',fs)
    end
end

if save_figs
    print(gcf,fullfile(depname_proc,'figures',sprintf('%s_attitude.png',berg)),'-dpng','-r300')
end

%% plot velocity stuff
save_figs = false;
dti = 2;
maxD = 1;
idxt = 1:dti:length(adcp.burst.time);
idxr = adcp.burst.range <= maxD;

flds = {'vel','cor','amp'};
fld_units = {'m/s','m/s','%','dB'};
CMAP = {'ba2','ba2','matter','amp'};
beam_ttls = {'right','up','toward1','error'};%{'right','up','left','down'};
ylbl = 'dist from ice (m)';
CLIM = {0.08*[-1 1],[0 100],extrema(adcp.burst.amp)};
ax = [];
for i = 1:3
    figure(i+1); clf;
    for j = 1:adcp.burst.nbeams
        subplot(adcp.burst.nbeams,1,j);
        ax(end+1) = gca;
        pcolor(adcp.burst.time(idxt),adcp.burst.range(idxr),adcp.burst.(flds{i})(idxt,idxr,j)')
        shading flat
        cmocean(CMAP{i})
        clim(CLIM{i})
        xlim(extrema(adcp.burst.time(idxt)))
        ylim(extrema(adcp.burst.range(idxr)))
        datetick('x','mmmdd HH:MM','keeplimits');
        if j == 1
            axpos1 = get(gca,'position');
            text(min(xlim),max(ylim)+0.1*diff(ylim),sprintf('spider - %s',berg_ttl),'fontsize',fs+2,'fontweight','bold');
            set(gca,'position',axpos1);
        elseif j == adcp.burst.nbeams
            axpos2 = get(gca,'position');
            cbar = colorbar;
            cbar.Position = [sum(axpos2([1 3]))+.01 axpos2(2) .015 sum(axpos1([2 4]))-axpos2(2)];
            cbar.Label.String = fld_units{i};
            cbar.Label.FontSize = fs;
        end
        ylabel(ylbl,'fontsize',fs)
        title(sprintf('%s: %s',strrep(flds{i},'_','\_'),beam_ttls{j}),'fontsize',fs)
    end
end
linkaxes(ax,'x')

if save_figs
    for i = 1:length(flds)
        fig_name = sprintf('%s_%s_preprocessed.png',berg,flds{i});
        print(figure(i+1),fullfile(depname_proc,'figures',fig_name),'-dpng','-r300')
    end
end

%% plot binned velocity
flds = {'vel_binned','vel_xyz_binned'};
fld_ttls = {'vel\_beam','vel\_xyz'};
beam_ttls = {{'right','up','left','down'},{'v_x (right)','v_y (up)','v_{z1} (toward ice)','v_{z1}-v_{z2} (error)'}};
CLIM = {0.05*[-1 1],0.065*[-1 1]};
maxD = [1 1];

% plot
for i = 1:2
    figure(3+i); clf;
    ax = [];
    idxr = adcp.burst.range <= maxD(i);
    for j = 1:4
        subplot(4,1,j);
        ax(end+1) = gca;
        pcolor(adcp.burst.tbins,adcp.burst.range(idxr),adcp.burst.(flds{i})(:,idxr,j)')
        shading flat
        cmocean('ba2')
        clim(CLIM{i})
        xlim(extrema(adcp.burst.tbins))
        ylim(extrema(adcp.burst.range(idxr)))
        datetick('x','mmmdd HH:MMZ','keeplimits');
        if j == 1
            axpos1 = get(gca,'position');
            text(min(xlim),max(ylim)+0.1*diff(ylim),sprintf('spider - %s',berg_ttl),'fontsize',fs,'fontweight','bold');
            set(gca,'position',axpos1);
        elseif j == adcp.burst.nbeams
            axpos2 = get(gca,'position');
            cbar = colorbar;
            cbar.Position = [sum(axpos2([1 3]))+.01 axpos2(2) .015 sum(axpos1([2 4]))-axpos2(2)];
            cbar.Label.String = 'm/s';
            cbar.Label.FontSize = fs;
        end
        ylabel(ylbl,'fontsize',fs)
        if j == 1
            title(sprintf('%s (%ds average): %s',fld_ttls{i},adcp.burst.bin_width,beam_ttls{i}{j}),'fontsize',fs)
        else
            title(beam_ttls{i}{j},'fontsize',fs)
        end
    end
    linkaxes(ax,'x')
end
