function [edges,tdrill] = findADVIceEdge(adv,ndrill,x0,L,L0,hann_width,varargin)
% Function to find ice location in ADV probe check echosounder profiles.
% All inputs after "adv" can be NaN to use default values.
%
% [edges,tdrill] = findADVIceEdge(adv,ndrills,x0,L,L0,hann_width)
% [edges,tdrill] = findADVIceEdge(adv,ndrills,x0,L,L0,hann_width,range_limit)
%
% Input
%   adv: structure containing pck data (output from parseADVDeployment)
%   x0: initial boundary location guess (scalar or 3-vector, NaN sets to
%       default value 200)
%   L: width of taper window (scalar or 3-vector, NaN sets to
%       default value 300)
%   L0: width of constant-value portion of taper window (scalar or 3-vector, NaN sets to
%       default value 100)
%   ndrill: number of drill events to attempt to find (NaN sets to default
%       value 0)
%   hann_width: width of Hanning filter applied along distance axis of
%       backscatter profiles (scalar or 3-vector, NaN sets to
%       default value 5)
%   range_limit: (optional) if scalar, minimum distance to look for edge,
%       if 2-element vector, second element is maximum distance
%
% Output
%   edges: structure with ice edge positions
%   tdrill: timestamps of when drilling likely occurred
%
% KJW
% 5 December 2023

% varargin
range_lim = [25 500];
if nargin == 7
    if numel(varargin{1}) == 1
        range_lim(1) = varargin{1};
    elseif numel(varargin{1}) == 2
        range_lim = varargin{1};
    else
        warning('Invalid range limits, reverting to default [25 500]')
    end
end

% set default parameters if nans are input
if isnan(ndrill)
    ndrill = 0;
end

if isnan(x0)
    x0 = 200;
end

if isnan(L)
    L = 300;
end

if isnan(L0)
    L0 = 100;
end

if isnan(hann_width)
    hann_width = 5;
end

% turn scalars into vectors for each beam
if isscalar(x0)
    x0 = x0*[1 1 1];
end
if isscalar(L)
    L = L*[1 1 1];
end
if isscalar(L0)
    L0 = L0*[1 1 1];
end
if isscalar(hann_width)
    hann_width = hann_width*[1 1 1];
end

% preallocate and redefine variables
edges = struct;
amp = adv.pck_amp;
dist = adv.pck_dist;
[nt,nc,nb] = size(amp);

%%% find drilling operations %%%
tdrill = NaT(ndrill,1);
if ndrill
    damp = abs(diff(amp,1,1));
    damp2_sum = sum(sum(damp.^2,3,'omitnan'),2);
    damp_rms = sqrt(damp2_sum);
    
    for i = 1:ndrill
        [~,imax] = max(damp_rms);
        tdrill(i) = adv.pck_time(imax);
        damp_rms(imax-3:imax+3) = 0;
    end
    tdrill = sort(tdrill);
end

%%% combine profiles separated by 5 sec %%%
dt = diff(adv.pck_time);
% check if amp needs to be front-padded
if dt(1)==seconds(5)
    burst1_idx = 2;
    amp_comb = amp;
    time_comb = adv.pck_time;
else
    burst1_idx = 1;
    amp_comb = cat(1,amp(1,:,:),amp);
    time_comb = cat(1,adv.pck_time(1),adv.pck_time);
end

% check if amp needs to be rear-padded
a = burst1_idx==1;
b = ~logical(mod(nt,2)); % is even
if (a && b) || (~a && ~b)
    amp_comb = cat(1,amp_comb,amp_comb(end,:,:));
    time_comb = cat(1,time_comb,time_comb(end));
end

% take mean
damp = diff(amp_comb,1,1);
amp = amp_comb(1:2:end,:,:)  + 0.5*damp(1:2:end,:,:);
time = time_comb(1:2:end);
nt = length(time);

% only use the burst-start distance axis
dist = dist(burst1_idx,:);
if isnan(dist(end))
    dist(end) = dist(end-1) + diff(dist(end-2:end-1));
end

% realign tdrill with new time axis
for i = 1:ndrill
    tdrill(i) = time(find((time-tdrill(i))<0,1,'last'));
end

%%% edge-finding %%%
% create smoothed amplitude matrix
for i = 1:nt
    for j = 1:nb
        amp(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width(j));
    end
end

% max jump flag setup
dz = min(abs(diff(dist)));
max_jump = 5*dz;

% find edges
x0_prev = x0;
epos = nan(nt,nb);
eidx = 0;
for i = 1:nt % loop through time
    for j = 1:nb % loop through beams
        % grab a single profile
        ampij = amp(i,:,j)';

        % taper profile around expected position, don't taper if meltstake
        % just drilled
        Lj = L(j);
        if x0(j)-x0_prev(j) <= -max_jump
            Lj = 1e6;
            fprintf('MAX JUMP!!\n')
        end
        if ~any(tdrill==time(i)) && eidx ~= 1
            ampij = taper(x0(j),Lj,L0(j),dist,ampij);
        end
        
        % find correct point in profile
        E_fld = centeredDifference(ampij);
        E_fld(dist<range_lim(1) | dist>range_lim(2)) = 0; % apply range limit
        [~,eidx] = max(E_fld);
        epos(i,j) = dist(eidx);

        % keep track of previous point(s)
        x0_prev(j) = x0(j);
        x0(j) = dist(eidx);
    end
end

% put output in structure
edges.pos = epos;
edges.desc = 'max CD (3pt) with hanning filter';
edges.hann_width = hann_width;
edges.time = time;

%%% subfunctions %%%
function prof = taper(x0,L,L0,dist,amp)
% Applies trapezoid weighting function to an amplitude profile.
    w = L/(L-L0) - 2/(L-L0)*abs(dist'-x0);
    w(w<0) = 0;
    w(w>1) = 1;
    prof = w.*amp;

function d = centeredDifference(v)
    d = [0; v(3:end)-v(1:end-2); 0];


