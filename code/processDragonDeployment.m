% Script to take raw dragon deployment data and turn it into useful
% profiles. I suspect that this will be a work-in-progress for quite some
% time as more data streams are incorporated, and we figure out the best
% way to present the data.
%
% TODO:
%
% KJW
% 13 Sep 2022

clear;

% set berg name
bergs = {'20220820_marchhare',...
         '20220821_marchhare',...
         '20220822_oyster',...
         '20220823_queenofhearts',...
         '20220824_singingflower',...
         '20220825_singingflower'};
berg = bergs{2};

fprintf('%s\n',berg);
save_processed_output = true;

% set paths
raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw';
mat_path = 'F:/Alaska2022/data/iceberg_surveys/mat';
proc_path = 'F:/Alaska2022/data/iceberg_surveys/proc';

% full paths
depname_raw = fullfile(raw_path,berg,'dragon');
depname_mat = fullfile(mat_path,berg,'dragon');
depname_proc = fullfile(proc_path,berg,'dragon');

if ~exist(depname_raw,'dir')
    warning('raw directory: %s not found, exiting...',depname_raw);
    return
end

% get rhib name
platforms = extractfield(dir(fullfile(raw_path,berg)),'name');
rhib_name = platforms{strcmpi(platforms,'aries') | strcmpi(platforms,'polly')};

% make directories if they don't exist
if ~exist(depname_mat,'dir')
    mkdir(depname_mat);
end
if ~exist(depname_proc,'dir')
    mkdir(depname_proc);
end

%% turn raw data into mat files, or load
% force raw data to be reparsed
force_parse = true;

% load profile file
times = load(fullfile(depname_raw,'profile_times.mat'));
n = length(times.t1);
t_start = times.t1(1)-60/86400;
t_end = times.t2(end)+60/86400;

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
    rbr = parseRBRDeployment(depname_raw,t_start,t_end);
    fprintf('saving %s\n',fullfile(depname_mat,'rbr.mat'));
    save(fullfile(depname_mat,'rbr.mat'),'rbr');
end

% ROV
fprintf('\nROV...\n');
if exist(fullfile(depname_mat,'rov.mat'),'file') && ~force_parse
    % load mat file
    fprintf('loading\n');
    rov = load(fullfile(depname_mat,'rov.mat'));
    rov = rov.rov;
else
    % make mat file
    fprintf('parsing\n');
    rov = parseROVTelemetry(depname_raw,t_start,t_end);
    fprintf('saving %s\n',fullfile(depname_mat,'rov.mat'));
    save(fullfile(depname_mat,'rov.mat'),'rov');
end
nt_rov = length(rov.time);

% ADCP
fprintf('\nADCP...\n');
if exist(fullfile(depname_mat,'adcp.mat'),'file') && ~force_parse
    % load mat file
    fprintf('loading\n');
    adcp = load(fullfile(depname_mat,'adcp.mat'));
    adcp = adcp.adcp;
    use_adcp = true;
else
    % make mat file
    fprintf('parsing\n');
    adcp = parseNortekDeployment(depname_raw,{},t_start,t_end);
    if ~isempty(fieldnames(adcp))
        use_adcp = true;
        fprintf('saving %s\n',fullfile(depname_mat,'adcp.mat'));
        save(fullfile(depname_mat,'adcp.mat'),'adcp');
    else
        use_adcp = false;
        fprintf('no adcp record found\n');
    end
end
% trim off echosounder data, too bulky
if use_adcp && isfield(adcp,'echo') && 0
    adcp = rmfield(adcp,'echo');
end
if use_adcp
    adcp_flds = fieldnames(adcp);
end

% RHIB GPS
fprintf('\nRHIB GPS...\n');
if exist(fullfile(depname_mat,'rhib_gps.mat'),'file') && ~force_parse
    % load mat file
    fprintf('loading\n');
    rhib_gps = load(fullfile(depname_mat,'rhib_gps.mat'));
    rhib_gps = rhib_gps.rhib_gps;
else
    % make mat file
    fprintf('parsing\n');
    rhib_gps = parseRHIBGPS(fullfile(raw_path,berg),rhib_name,t_start,t_end);
    fprintf('saving %s\n',fullfile(depname_mat,'rhib_gps.mat'));
    save(fullfile(depname_mat,'rhib_gps.mat'),'rhib_gps');
end
nt_rhib = length(rhib_gps.time);

%% process data streams into profiles
clear profs
profs(n) = struct();

% % rake/ctd sensor positions
rake_idx = ~strcmp('concerto',extractfield(rbr,'model'));
inst_idx = find(rake_idx);
ctd_idx = find(~rake_idx,1);
pos = extractfield(rbr(rake_idx),'pos');

% this helps indexing a slightly noisy time vector
% the fastest time vector is the ctd at 16Hz (.0625 sec)
t_pad = .015/86400;

% loop through profiles
% each profs(i) has fields
%   .rake
%   .ctd
%   .rov
%   .gps
%   (.adcp)
%
% QC for vehicle motion and attitude
%   1) pitch within a certain range
%   2) roll within a certain range
%   3) depth increasing
%   4) yaw rate below a certain value
qc = struct();
qc.pitch_range = 12; % deg
qc.roll_range = 12; % deg
qc.yawrate_range = 12; % deg/s
qc.fallrate_min = .03; % dbar/s
qc.fallrate_max = .15; % dbar/s

verbose = true; % debugging
for i = 1:n
    %%% build QC vector %%%
    
    % index into data streams
    rov_slice = rov.time >= times.t1(i)-0.5/86400-t_pad & rov.time < times.t2(i)-0.5/86400; % rov timestamps are not locked to full seconds, so center window around the full second to go for closest point
    ctd_qc_slice = rbr(ctd_idx).time >= times.t1(i)-t_pad-2/86400 & rbr(ctd_idx).time < times.t2(i)+2/86400;
    ctd_slice = rbr(ctd_idx).time >= times.t1(i)-t_pad & rbr(ctd_idx).time < times.t2(i);
    rake_slice = rbr(inst_idx(1)).time>=times.t1(i)-t_pad & rbr(inst_idx(1)).time < times.t2(i)+t_pad;
    rhib_slice = rhib_gps.time >= times.t1(i)-t_pad & rhib_gps.time < times.t2(i);
    if use_adcp % adcp (by channel)
        adcp_slice = struct();
        for j = 1:length(adcp_flds)
            fldj = adcp_flds{j};
            if isfield(adcp.(fldj),'time')
                adcp_slice.(fldj) = adcp.(fldj).time>=times.t1(i)-t_pad & adcp.(fldj).time<=times.t2(i);
            end
        end
    end
    
    % fastest vector to interpolate onto
    qc.ctd_time = rbr(ctd_idx).time(ctd_qc_slice);
    qc.ctd_depth = rbr(ctd_idx).values(ctd_qc_slice,strcmp(rbr(ctd_idx).channels,'Depth'));
    ctd_sr = mode(round(1./(86400*diff(qc.ctd_time))));
    
    % get fields (no need to use angle interpolation since 'nearest'
    % algorithm is used
    qc.pitch = interp1(rov.time,rov.pitch,qc.ctd_time,'nearest','extrap');
    qc.roll = interp1(rov.time,rov.roll,qc.ctd_time,'nearest','extrap');
    qc.yawrate = interp1(rov.time,rov.yawRate,qc.ctd_time,'nearest','extrap');
    qc.pitchrate = interp1(rov.time,rov.pitchRate,qc.ctd_time,'nearest','extrap');
    qc.rollrate = interp1(rov.time,rov.rollRate,qc.ctd_time,'nearest','extrap');
    
    % remove mode from pitch/roll
    qc.pitch = qc.pitch - mode(round(qc.pitch,1));
    qc.roll = qc.roll - mode(round(qc.roll,1));
      
    % calculate fall rate from averaged ctd depth
    qc.fallrate = diff(meanFilter(qc.ctd_depth,ctd_sr+1))./(86400*diff(qc.ctd_time)); qc.fallrate = [qc.fallrate; qc.fallrate(end)];
    qc.fallrate = meanFilter(qc.fallrate,ctd_sr+1);
    
    % evaluate QC criteria
    qc.bad = false(size(qc.ctd_time));
    qc.bad = qc.bad | abs(qc.pitch)>qc.pitch_range;
    qc.bad = qc.bad | abs(qc.roll)>qc.roll_range;
    qc.bad = qc.bad | abs(qc.yawrate)>qc.yawrate_range;
    qc.bad = qc.bad | qc.fallrate<qc.fallrate_min;
    qc.bad = qc.bad | qc.fallrate>qc.fallrate_max;
    qc.bad_frac = length(find(qc.bad))/length(qc.bad);
    
    % subsample to each data stream (might not need rov)
    qc.bad_ctd = logical(interp1(qc.ctd_time,qc.bad*1,rbr(ctd_idx).time(ctd_slice),'nearest'));
    qc.bad_rake = logical(interp1(qc.ctd_time,qc.bad*1,rbr(inst_idx(1)).time(rake_slice),'nearest'));
    qc.bad_rov = logical(interp1(qc.ctd_time,qc.bad*1,rov.time(rov_slice),'nearest'));
    if use_adcp
        qc.bad_adcp = struct();
        for j = 1:length(adcp_flds)
            fldj = adcp_flds{j};
            if isfield(adcp.(fldj),'time')
                qc.bad_adcp.(fldj) = logical(interp1(qc.ctd_time,qc.bad*1,adcp.(fldj).time(adcp_slice.(fldj)),'nearest'));
            end
        end
    end
    
    % save QC information in each profile
    profs(i).qc = qc;
    
    %%% build profile %%%
    
    % use depth from ctd
    
    % ctd
    ctd_flds = {'depth','temperature','salinity','conductivity'};
    profs(i).ctd = struct();
    profs(i).ctd.pos = rbr(ctd_idx).pos;
    profs(i).ctd.time = rbr(ctd_idx).time(ctd_slice);
    for tmp = ctd_flds
        fld = tmp{1};
        profs(i).ctd.(lower(fld)) = rbr(ctd_idx).values(ctd_slice,strcmpi(rbr(ctd_idx).channels,fld));
    end
    profs(i).ctd.fallrate = interp1(qc.ctd_time,qc.fallrate,profs(i).ctd.time,'nearest');
    
    % rake
    profs(i).rake = struct();
    profs(i).rake.pos = pos;
    profs(i).rake.time = rbr(inst_idx(1)).time(rake_slice); % use first instrument time vector (they're all the same)
    profs(i).rake.depth = interp1(profs(i).ctd.time,profs(i).ctd.depth,profs(i).rake.time,'previous');
    m = length(profs(i).rake.time);
    profs(i).rake.temperature = NaN(m,length(pos)); % preallocate temperature matrix
    profs(i).rake.temperature(:,1) = rbr(inst_idx(1)).values(rake_slice,strcmp(rbr(inst_idx(1)).channels,'Temperature'));
    for j = 2:length(inst_idx) % loop through the rest
        rake_slice = rbr(inst_idx(j)).time>=times.t1(i)-t_pad & rbr(inst_idx(j)).time<=times.t2(i);
        profs(i).rake.temperature(:,j) = rbr(inst_idx(j)).values(rake_slice,strcmp(rbr(inst_idx(j)).channels,'Temperature'));
    end

    % rov
    profs(i).rov = struct();
    rov_flds = fieldnames(rov)';
    for tmp = rov_flds
        fld = tmp{1};
        if length(rov.(fld))==nt_rov
            profs(i).rov.(fld) = rov.(fld)(rov_slice);
        end
    end
    profs(i).rov.depth = interp1(profs(i).ctd.time,profs(i).ctd.depth,profs(i).rov.time,'previous');
    
    % clean up heading a little bit to minimize wrapping
    if max(abs(diff(profs(i).rov.heading)))>350
        cuts = min(extrema(profs(i).rov.heading)):max(extrema(profs(i).rov.heading));
        performance = nan*(cuts);
        for j = 1:length(cuts)
            new_hdg = profs(i).rov.heading;
            if cuts(j) < 0
                new_hdg(new_hdg<cuts(j)) = new_hdg(new_hdg<cuts(j))+360;
            else
                new_hdg(new_hdg>=cuts(j)) = new_hdg(new_hdg>=cuts(j))-360;
            end
            performance(j) = max(abs(diff(new_hdg)));
        end
        [~,cut_idx] = min(performance);
        cut_best = cuts(cut_idx);
        if cut_best < 0 
            profs(i).rov.heading(profs(i).rov.heading<cut_best) = profs(i).rov.heading(profs(i).rov.heading < cut_best)+360;
        else
            profs(i).rov.heading(profs(i).rov.heading>=cut_best) = profs(i).rov.heading(profs(i).rov.heading>=cut_best)-360;
        end
    end

    % adcp
    if use_adcp
        profs(i).adcp = struct();
        for j = 1:length(adcp_flds)
            fldj = adcp_flds{j};
            if isfield(adcp.(fldj),'time') % channel has a time vector, split it into profiles
                subflds = fieldnames(adcp.(fldj));
                for k = 1:length(subflds)
                    fldk = subflds{k};
                    nd1 = size(adcp.(fldj).(fldk),1);
                    if nd1 == length(adcp_slice.(fldj)) % time slice
                        profs(i).adcp.(fldj).(fldk) = adcp.(fldj).(fldk)(adcp_slice.(fldj),:,:);
                    else % copy full subfield
                        profs(i).adcp.(fldj).(fldk) = adcp.(fldj).(fldk);
                    end
                end
                if ~isfield(adcp.(fldj),'pressure') % if no pressure axis, make one by interpolating
                    profs(i).adcp.(fldj).pressure = interp1(adcp.attitude.time,adcp.attitude.pressure,profs(i).adcp.(fldj).time);
                end
            else % else just copy the whole thing
                profs(i).adcp.(fldj) = adcp.(fldj);
            end
        end
        % change roll from [-90, 270) to [0, 360)
        if isfield(profs(i).adcp,'attitude') && isfield(profs(i).adcp.attitude,'roll')
            adcp_roll = profs(i).adcp.attitude.roll;
            adcp_roll(adcp_roll<0) = adcp_roll(adcp_roll<0)+360;
            profs(i).adcp.attitude.roll = adcp_roll;
        end
    end
    
    % rhib gps
    profs(i).rhib_gps = struct();
    rhib_flds = fieldnames(rhib_gps)';
    for tmp = rhib_flds
        fld = tmp{1};
        if length(rhib_gps.(fld))==nt_rhib
            profs(i).rhib_gps.(fld) = rhib_gps.(fld)(rhib_slice);
        else
            profs(i).rhib_gps.(fld) = rhib_gps.(fld);
        end
    end    
    
    %%% apply qc %%%
    platforms = {'rake'};
    for j = 1:length(platforms)
        platform = platforms{j};
        m = length(profs(i).(platform).time);
        
        bad_blocks = findBlocks(profs(i).qc.(['bad_' platform]));
        if ~isempty(bad_blocks)
            rm_vec = [];
            nan_vec = bad_blocks(:,1);
            for k = 1:size(bad_blocks,1) % build vector of line indices to remove (leave one for nan-ing)
                rm_vec = [rm_vec (bad_blocks(k,1)+1):bad_blocks(k,2)];
            end

            flds_all = fieldnames(profs(i).(platform))';
            for tmp = flds_all % nan out one line and chop the rest
                fld = tmp{1};
                if size(profs(i).(platform).(fld),1)==m
                    profs(i).(platform).(fld)(nan_vec,:) = nan;
                    profs(i).(platform).(fld)(rm_vec,:) = [];
                end
            end
        end
    end
    
    % adcp
    if use_adcp
        for j = 1:length(adcp_flds)
            fldj = adcp_flds{j};
            if ~isfield(profs(i).adcp.(fldj),'time') || strcmp(fldj,'attitude')
                continue
            end
            m = length(profs(i).adcp.(fldj).time);
            bad_blocks = findBlocks(profs(i).qc.bad_adcp.(fldj));
            
            if ~isempty(bad_blocks)
                rm_vec = [];
                nan_vec = bad_blocks(:,1);
                for k = 1:size(bad_blocks,1) % build vector of line indices to remove (leave one for nan-ing)
                    rm_vec = [rm_vec (bad_blocks(k,1)+1):bad_blocks(k,2)];
                end
                
                subflds = fieldnames(profs(i).adcp.(fldj))';
                for tmp = subflds
                    fldk = tmp{1};
                    if size(profs(i).adcp.(fldj).(fldk),1)==m
                        profs(i).adcp.(fldj).(fldk)(nan_vec,:,:) = nan;
                        profs(i).adcp.(fldj).(fldk)(rm_vec,:,:) = [];
                    end
                end
            end
        end
    end
end

%% save processed output
if save_processed_output
    proc_fname = fullfile(depname_proc,'profiles.mat');
    fprintf('saving file %s\n',proc_fname)
    save(proc_fname,'profs');
end

% debugging
if verbose
    fprintf('Fraction removed from each profile (time)\n');
    for i = 1:n
        fprintf('%d) %.2f\n',i,profs(i).qc.bad_frac)
    end
end

