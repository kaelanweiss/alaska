function [adv,raw] = parseADVDeployment(depname,varargin)
% Function to parse Vector ADV deployment.
%
% [adv,raw] = parseADVDeployment(depname)
% [adv,raw] = parseADVDeployment(depname,t1,t2)
%
% Input
%   depname: path to folder containing files output from Vector software
%       after parsing .VEC file (i.e. .dat, .hdr, .pck, .vhd, .sen, .ssl)
%   t1: (optional) start of time window to parse (datetime)
%   t2: (optional) end of time window to parse (datetime)
% 
% KJW
% 16 Sep 2023

%%% Check that directory exists %%%
if ~exist(depname,'dir')
    warning('directory does not exist: %s\n',depname);
    return
end

%%% Look for instrument.txt file %%%
if exist(fullfile(depname,'','instruments.txt'),'file')
    info = parseInstrumentFile(depname,'');
    fprintf('instrument file found, using time offset dt=%.1f\n',info.t_offset);
    if ~isnan(info.t_offset)
        dt = info.t_offset;
    else
        dt = seconds(0);
    end
    use_inst_file = true;
else
    fprintf('instrument file not found\n');
    dt = seconds(0);
    use_inst_file = false;
end

%%% Parse varargin %%%
% Optional time window
if nargin == 1
    t1 = datetime(0,1,1);
    t2 = datetime(3000,1,1);
else
    t1 = varargin{1};
    t2 = varargin{2};
end

% Grab header file
d = dir(fullfile(depname,'*.hdr'));
dname = d.name(1:end-4);

% Set parsing formats for file types
exts = {'vhd','sen','pck'};
fmts = {'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f',...
        '%d %d %d %d %d %d %b %b %f %f %f %f %f %f %d %d',...
        '%d %f %d %d %d'};
dat_fmt = '%d %d %f %f %f %d %d %d %f %f %f %d %d %d %f %f %f %d';

% Preallocate raw output
raw = struct();

%%% parse specific lines from header file %%%
% header file
raw.hdr = struct();
fid = fopen(fullfile(depname,[dname '.hdr']),'r');
hdr_lines = {};
while ~feof(fid)
    hdr_lines{end+1} = fgetl(fid);
end
fclose(fid);

% sample rate
i = 1;
while isempty(regexp(hdr_lines{i},'Sampling rate','once'))
    i = i + 1;
end
sr_str = strsplit(hdr_lines{i},' ');
sr = str2double(strip(sr_str{end-1}));

% number of measurements
i = 1;
while isempty(regexp(hdr_lines{i},'Number of measurements','once'))
    i = i + 1;
end
nmeas_str = strsplit(hdr_lines{i},' ');
nmeas = str2double(strip(nmeas_str{end}));

% velocity range
i = 1;
while isempty(regexp(hdr_lines{i},'Nominal velocity range','once'))
    i = i + 1;
end
vr_str = strsplit(hdr_lines{i},' ');
vel_range = str2double(strip(vr_str{end-1}));

% burst interval
i = 1;
while isempty(regexp(hdr_lines{i},'Burst interval','once'))
    i = i + 1;
end
bi_str = strsplit(hdr_lines{i},' ');
burst_interval = str2double(strip(bi_str{end-1}));

% samples per burst
i = 1;
while isempty(regexp(hdr_lines{i},'Samples per burst','once'))
    i = i + 1;
end
spb_str = strsplit(hdr_lines{i},' ');
samples_per_burst = str2double(strip(spb_str{end}));

% sampling volume
i = 1;
while isempty(regexp(hdr_lines{i},'Sampling volume','once'))
    i = i + 1;
end
sv_str = strsplit(hdr_lines{i},' ');
sampling_volume = str2double(strip(sv_str{end-1}));

% transmit length
i = 1;
while isempty(regexp(hdr_lines{i},'Transmit length','once'))
    i = i + 1;
end
txl_str = strsplit(hdr_lines{i},' ');
transmit_length = str2double(strip(txl_str{end-1}));

% beam to instrument transformation matrix
i = 1;
while isempty(regexp(hdr_lines{i},'Transformation matrix','once'))
    i = i + 1;
end
b2i_lines = hdr_lines(i:i+2);
beam2xyz = nan(3);
for i = 1:3
    C = strsplit(b2i_lines{i},' ');
    beam2xyz(i,:) = str2double(C(end-2:end));
end

% preallocate dat
Cdatfmt = strsplit(dat_fmt,' ');
raw.dat = cell(1,length(Cdatfmt));
for i = 1:length(Cdatfmt)
    if strcmp(Cdatfmt{i},'%f')
        raw.dat{i} = -99*ones(1,nmeas,'double')';
    elseif strcmp(Cdatfmt{i},'%d')
        raw.dat{i} = -99*ones(1,nmeas,'int32')';
    end
end

% parse dat files
dat_files = dir(fullfile(depname,[dname '*.dat']));
pos = 1;
for i = 1:length(dat_files)
    fid = fopen(fullfile(depname,dat_files(i).name),'r');
    dati = textscan(fid,dat_fmt);
    ni = length(dati{1});
    for j = 1:length(dati)
        raw.dat{j}(pos:pos+ni-1) = dati{j};
    end
    pos = pos + ni;
    fclose(fid);
    fprintf('%d/%d\n',i,length(dat_files))
end

% other files
for i = 1:length(exts)
    if ~strcmp(exts{i},'dat')
        fid = fopen(fullfile(depname,[dname '.' exts{i}]),'r');
        raw.(exts{i}) = textscan(fid,fmts{i});
        fclose(fid);
    end
end

% sort and name
adv = struct();
if use_inst_file
    adv.timeCoordinationOffset = dt;
else
    adv.timeCoordinationOffset = NaT;
end
adv.beam2xyz = beam2xyz;
adv.sample_rate = sr;
adv.n_measurements = nmeas;
adv.vel_range = vel_range;
adv.burst_interval = burst_interval;
adv.samples_per_burst = samples_per_burst;
adv.sampling_volume = sampling_volume;
adv.transmit_length = transmit_length;

% dat
adv.burst_num = double(raw.dat{1});
adv.ens_num = double(raw.dat{2});
adv.vel = [raw.dat{3} raw.dat{4} raw.dat{5}];
adv.amp = [raw.dat{6} raw.dat{7} raw.dat{8}];
adv.snr = [raw.dat{9} raw.dat{10} raw.dat{11}];
adv.cor = [raw.dat{12} raw.dat{13} raw.dat{14}];
adv.pressure = raw.dat{15};
adv.chk = raw.dat{18};

% vhd
adv.hdr_time = datetime([raw.vhd{3} raw.vhd{1} raw.vhd{2} raw.vhd{4} raw.vhd{5} raw.vhd{6}]) + dt;
adv.dp = cat(3,[raw.vhd{15} raw.vhd{16} raw.vhd{17}],[raw.vhd{20} raw.vhd{21} raw.vhd{22}]);
adv.dprobe_ave = [raw.vhd{18} raw.vhd{23}];
adv.dsvol_ave = [raw.vhd{19} raw.vhd{24}];

% sen
adv.sen_time = datetime([raw.sen{3} raw.sen{1} raw.sen{2} raw.sen{4} raw.sen{5} raw.sen{6}]) + dt;
adv.batteryvoltage = raw.sen{9};
adv.soundspeed = raw.sen{10};
adv.hdg = raw.sen{11};
adv.ptc = raw.sen{12};
adv.rol = raw.sen{13};
adv.temp = raw.sen{14};

% pck
n_pck = length(raw.pck{1});
idx_prof = [1; find(diff(raw.pck{1})<0)+1; n_pck+1];
n_prof = length(idx_prof)-1;
n_cell = max(diff(idx_prof));
adv.pck_amp = nan(n_prof,n_cell,3);
adv.pck_dist = nan(n_prof,n_cell);
adv.pck_time = NaT(n_prof,1);
burst_length = seconds(adv.samples_per_burst/adv.sample_rate);
for i = 1:n_prof
    slice = idx_prof(i):idx_prof(i+1)-1;
    ni = length(slice);
    adv.pck_dist(i,1:ni) = raw.pck{2}(slice);
    adv.pck_amp(i,1:ni,:) = [raw.pck{3}(slice) raw.pck{4}(slice) raw.pck{5}(slice)];
    adv.pck_time(i) = adv.hdr_time(floor((i+1)/2)) + mod(i+1,2)*burst_length;
end

% time vector
adv.time = adv.hdr_time(adv.burst_num) + seconds((1/sr)*(adv.ens_num-1));
adv.burst_time = adv.time(logical([1;diff(adv.burst_num)]));

% trim down to chosen time frame (this is not the best way to do this since
% you end up parsing way more than you need, but oh well)
flds = fields(adv);
time_flds = {'hdr_time','sen_time','pck_time','time','burst_time'};
for i = 1:length(time_flds) % loop through time fields
    nt = length(adv.(time_flds{i}));
    idx = adv.(time_flds{i}) >= t1 & adv.(time_flds{i}) <= t2;
    for j = 1:length(flds) % oh god this is a terrible way to do this
        if size(adv.(flds{j}),1) == nt
            adv.(flds{j}) = adv.(flds{j})(idx,:,:);
        end
    end
end

end