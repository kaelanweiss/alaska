function [edges,tdrill] = findADVIceEdge(adv,ndrill,x0,L,L0,hann_width,varargin)
% Function to find ice location in ADV probe check echosounder profiles.
%
% [edges,tdrill] = findADVIceEdge(adv,ndrills,x0,L,L0,hann_width)
% [edges,tdrill] = findADVIceEdge(adv,ndrills,x0,L,L0,hann_width,min_dist)
%
% Input
%   adv: structure containing pck data (output from parseADVDeployment)
%   x0: 
%   L: width of taper window
%   L0: width of constant-value portion of taper window
%   ndrill: 
%   hann_width: 
%   range_limit: (optional) if scalar, minimum distance to look for edge,
%       if 2-element vector, second element is maximum distance
%
% Output
%   edges: 
%   tdrill: 
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
if length(x0)==1
    x0 = x0*[1 1 1];
end
if length(L)==1
    L = L*[1 1 1];
end
if length(L0)==1
    L0 = L0*[1 1 1];
end

% preallocate and redefine variables
edges = struct;
amp = adv.pck_amp;
dist = adv.pck_dist;
[nt,nc,nb] = size(amp);
% amp(repmat(dist,[1 1 nb])<min_dist) = nan;

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

%%% combine time-adjacent profiles %%%
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
% centered-difference d/dz operator

% derivative operator (3pt stencil, unscaled)
D3 = zeros(nc) - diag(ones(nc-1,1),-1) + diag(ones(nc-1,1),1);
D3([1 end],:) = 0;

% create smoothed amplitude matrix
for i = 1:nt
    for j = 1:nb
        amp(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
    end
end

% look-behind setup
% m = 8; % number of previous points
% k = 1; % steepness of weighting (0: flat, ~10: previous point)
% 
% x = -(m-1):0;
% A = sum(exp(k*x))^-1;
% wlb = A*exp(k*x);
% x0 = repmat(x0,[m 1]);
% drill_flag = 0;
dz = min(abs(diff(dist)));
x0_prev = x0;
max_jump = 5*dz;

% find edges
epos = nan(nt,nb);
for i = 1:nt
    for j = 1:nb
        % grab a single profile and fill nans
        ampij = amp(i,:,j)';
%         ampij(isnan(ampij)) = ampij(find(~isnan(ampij),1));
        
%         if j == 2
%             figure(99); clf; hold on
%             plot(dist,ampij,'k.-') % debug
%         end

        % taper profile around expected position
        % x0j = wlb*x0(:,j);
        Lj = L(j);
        if x0(j)-x0_prev(j) <= -max_jump
            Lj = 1e6;
            fprintf('MAX JUMP!!\n')
        end

        if ~any(tdrill==time(i))
            ampij = taper(x0(j),Lj,L0(j),dist,ampij);
%         else
%             drill_flag = 1;
%             idxt = (1:ndrill)*(tdrill==time(i));
%             fprintf('drill %d: %s == %s\n',idxt,tdrill(idxt),time(i))
        end
        
        % find correct point
        cent_diff = D3*ampij;
        cent_diff(dist<range_lim(1) | dist>range_lim(2)) = 0;
        [~,eidx] = max(cent_diff);
        epos(i,j) = dist(eidx);

        % keep track of previous point(s)
        x0_prev(j) = x0(j);
        x0(j) = dist(eidx);
%         x0(:,j) = [dist(eidx); x0(1:m-1,j)];
%         if drill_flag
%             x0(:,j) = x(1,j);
%             drill_flag = 0;
%         end

%         if j == 2
%             plot(dist,ampij,'r.-') % debug
%             plot(x0(j)*[1 1],[0 max(ampij)],'b--') % debug
%             title(sprintf('%s',time(i)))
%         end
    end
    a=1;
end
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


