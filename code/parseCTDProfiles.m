function ctd = parseCTDProfiles(rsk_files,dz,varargin)
% Function to read RSK CTD files, find downcast profiles, perform some
% light processing, and bin vertically.
%
% ctd = parseCTDProfiles(rsk_files,dz)
% ctd = parseCTDProfiles(rsk_files,dz,tlim)
%
% Input
%
%
% Output
%
%
% KJW
% 1 Feb 2024

% parse varargin
t1 = 0;
t2 = 1e6;
switch nargin
    case 3
        tlim = varargin{1};
        t1 = datenum(tlim(1));
        t2 = datenum(tlim(2));
end

% deal with single vs. cell array of files
if isstr(rsk_files)
    rsk_files = {rsk_files};
end
nfiles = length(rsk_files);

% preallocate ctd output
data = cell(1,nfiles);
time = cell(1,nfiles);
z = cell(1,nfiles);
z_all = 0;
np_all = 0;

% channels to bin (depth last as binning axis)
flds = {'Temperature','Conductivity','Salinity','Density Anomaly','Speed of Sound','Depth'};

for i = 1:nfiles
    % open RSK file
    rsk = RSKopen(rsk_files{i});
    rsk = RSKreaddata(rsk,'t1',t1,'t2',t2); % this should detect profiles

    % continue if no data within time interval
    if ~isfield(rsk,'data')
        continue
    end
    
    % check that profiles were found, try again, then continue
    if ~isfield(rsk,'profiles') || length(rsk.profiles.downcast) <= 1 % weird edge case
        rsk = RSKfindprofiles(rsk);
    end
    if ~isfield(rsk,'profiles')
        continue
    end

    % adjust z=0 level, recalculate depth
    if ~any(strcmp('Sea Pressure',{rsk.channels.longName}))
        rsk = RSKderiveseapressure(rsk);
    end
    sp_chn = find(strcmp('Sea Pressure',{rsk.channels.longName}));
    sp = rsk.data.values(:,sp_chn);
    sp_surf = sp(sp<1);
    sp0 = median(round(sp_surf,1));
    rsk.data.values(:,sp_chn) = sp - sp0;
    rsk = RSKderivedepth(rsk);

    % derive missing channels
    if ~any(strcmp('Salinity',{rsk.channels.longName}))
        rsk = RSKderivesalinity(rsk);
    end
    if ~any(strcmp('Density Anomaly',{rsk.channels.longName}))
        rsk = RSKderivesigma(rsk);
    end
    if ~any(strcmp('Speed of Sound',{rsk.channels.longName}))
        rsk = RSKderivesoundspeed(rsk);
    end
    
    % useful values
    np = length(rsk.profiles.downcast.tstart);
    sp_chn = find(strcmp('Sea Pressure',{rsk.channels.longName}));
    sp_max = max(rsk.data.values(:,sp_chn));
    t_downcast = [rsk.profiles.downcast.tstart rsk.profiles.downcast.tend];

    % despike/de-lag

    % preallocate (T, C, S, sigma, soundspeed)
    zi = (0:dz:(sp_max+10))';
    nz = length(zi);
    datai = nan(nz,np,5);
    timei = nan(nz,np);

    % bin vertically
    ch_num = nan(size(flds));
    for j = 1:length(flds)
        ch_num(j) = find(strcmpi(flds{j},{rsk.channels.longName}));
    end
    % loop through profiles
    for j = 1:np
        idxj = rsk.data.tstamp>=t_downcast(j,1) & rsk.data.tstamp<=t_downcast(j,2);
        dataj = rsk.data.values(idxj,ch_num(1:end-1));
        depthj = rsk.data.values(idxj,ch_num(end));
        timej = rsk.data.tstamp(idxj);
        % vert bins
        binnedj = vertBin(dataj,depthj,zi);
        datai(:,j,:) = permute(binnedj,[1 3 2]);
        timei(:,j) = vertBin(timej,depthj,zi);
    end

    data{i} = datai;
    time{i} = timei;
    z{i} = zi;
    
    % keep track of total profiles and longest depth axis
    np_all = np_all + np;
    if max(z_all) < max(zi)
        z_all = zi;
    end
end

% join all data
ctd.z = z_all;
nz = length(z_all);

ctd.time = dn2dt(joinData(time,nz));
ctd.t0 = min(ctd.time,[],1,'omitnan');

data_joined = joinData(data,nz);
ctd.T = data_joined(:,:,1);
ctd.C = data_joined(:,:,2);
ctd.S = data_joined(:,:,3);
ctd.sigma = data_joined(:,:,4);
ctd.soundspeed = data_joined(:,:,5);

% get rid of duplicate profiles from overlapping files (this is
% inefficient)
[ctd.t0,uidx] = unique(ctd.t0);
ctd.time = ctd.time(:,uidx);

flds_final = {'T','C','S','sigma','soundspeed'};
for i = 1:length(flds_final)
    fldi = flds_final{i};
    ctd.(fldi) = ctd.(fldi)(:,uidx);
end

%%% subfunctions %%%

function out = vertBin(in,depth,z)
% Function to vertically bin data. Cells are centered around z-values.
% Input has dimensions nt x nchannels

% Useful values
nz = length(z);
nc = size(in,2);
dz = z(2)-z(1);

% preallocate
out = nan(nz,nc);

% loop over depth bins
for i = 1:nz
    % get indices in depth bin
    idxz = depth>=(z(i)-0.5*dz) & depth<(z(i)+0.5*dz);
    % skip bin if there's no data
    if ~sum(idxz)
        continue
    end

    % loop through channels
    for j = 1:nc
        % skip if there are a bunch of nans
        if sum(isnan(in(idxz,j)))/sum(idxz) > 0.75
            continue
        end
        % calculate mean
        out(i,j) = mean(in(idxz,j),'omitnan');
    end
end

function joined = joinData(data_cells,nz)
% Function to join data saved in cells into a single matrix.
    joined = [];
    n_cells = length(data_cells);
    for i = 1:n_cells
        if ~isempty(data_cells{i})
            % pad to match z length
            nzi = size(data_cells{i},1);
            pad = nan(nz-nzi,size(data_cells{i},2),size(data_cells{i},3));
            data_pad = cat(1,data_cells{i},pad);

            % cat along 2nd (profile) dimension
            joined = cat(2,joined,data_pad);
        end
    end






