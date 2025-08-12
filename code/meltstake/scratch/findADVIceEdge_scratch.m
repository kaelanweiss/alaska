function [edges,tdrill] = findADVIceEdge(adv,varargin)
% Function to find ice location in ADV probe check echosounder profiles.
%
% format
%
% Input
%   adv: structure containing pck data (output from parseADVDeployment)
%   x0: (optional)
%   L: (optional)
%   ndrill: (optional)
%
% Output
%   edges: 
%   tdrill: 
%
% KJW
% 5 December 2023

% parse varargin
switch nargin
    case 1 % default values
        x0 = 150*[1 1 1];
        L = 1e6*[1 1 1];
        ndrill = 0;
    case 2 % only ndrills specified
        x0 = 150*[1 1 1];
        L = 1e6*[1 1 1];
        ndrill = varargin{1};
    case 3 % only window parameters specified
        x0 = varargin{1};
        L = varargin{2};
        ndrill = 0;
        if length(x0)==1
            x0 = x0*[1 1 1];
        end
        if length(L)==1
            L = L*[1 1 1];
        end
    case 4 % window parameters and ndrills specified
        x0 = varargin{1};
        L = varargin{2};
        if length(x0)==1
            x0 = x0*[1 1 1];
        end
        if length(L)==1
            L = L*[1 1 1];
        end
        ndrill = varargin{3};
end

% preallocate and redefine variables
edges = struct;
amp = adv.pck_amp;
dist = adv.pck_dist;
[nt,nc,nb] = size(amp);
amp(repmat(dist,[1 1 nb])<25) = nan;

%%% find drilling operations %%%
tdrill = NaT(ndrill,1);
if ndrill
    damp = abs(diff(amp,1,1));
    damp2_sum = sum(sum(damp.^2,3,'omitnan'),2);
    damp_rms = sqrt(damp2_sum);
    
    for i = 1:ndrill
        [~,imax] = max(damp_rms);
        tdrill(i) = adv.pck_time(imax);
        damp_rms([imax-3:imax+3]) = 0;
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

%%% edge-finding with varying techniques %%%
% % 1. max of profile
% amp1 = amp;
% 
% [~,eidx1] = max(amp1,[],2);
% eidx1 = squeeze(eidx1);
% 
% epos1 = nan*eidx1;
% for i = 1:nt
%     for j = 1:nb
%         epos1(i,j) = dist(eidx1(i,j));
%     end
% end
% out(1).pos = epos1;
% out(1).desc = 'max of profile';
% 
% figure(201); clf
% plot(time,epos1)
% 
% % 2. max positive difference in profile
% amp2 = amp;
% 
% [~,eidx2] = max(diff(amp2,1,2),[],2);
% eidx2 = squeeze(eidx2);
% 
% epos2 = nan*eidx2;
% for i = 1:nt
%     for j = 1:nb
%         epos2(i,j) = dist(eidx2(i,j));
%     end
% end
% out(2).pos = epos2;
% out(2).desc = 'max positive difference in profile';
% 
% figure(202); clf
% plot(time,epos2)
% 
% % 3. max of profile with hanning filter
% hann_width = 5;
% amp3 = amp;
% 
% for i = 1:nt
%     for j = 1:nb
%         amp3(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
%     end
% end
% 
% [~,eidx3] = max(amp3,[],2);
% eidx3 = squeeze(eidx3);
% 
% epos3 = nan*eidx3;
% for i = 1:nt
%     for j = 1:nb
%         epos3(i,j) = dist(eidx3(i,j));
%     end
% end
% out(3).pos = epos3;
% out(3).desc = 'max of profile with hanning filter';
% out(3).hann_width = hann_width;
% 
% figure(203); clf
% plot(time,epos3)
% 
% % 4. max diff of profile with hanning filter
% hann_width = 5;
% amp4 = amp;
% 
% for i = 1:nt
%     for j = 1:nb
%         amp4(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
%     end
% end
% 
% [~,eidx4] = max(diff(amp4,1,2),[],2);
% eidx4 = squeeze(eidx4);
% 
% epos4 = nan*eidx4;
% for i = 1:nt
%     for j = 1:nb
%         epos4(i,j) = dist(eidx4(i,j));
%     end
% end
% out(4).pos = epos4;
% out(4).desc = 'max positive diff of profile with hanning filter';
% out(4).hann_width = hann_width;
% 
% figure(204); clf
% plot(time,epos4)

% 5. centered-difference d/dz operator
% derivative operator (3pt stencil, unscaled)
D3 = zeros(nc) - diag(ones(nc-1,1),-1) + diag(ones(nc-1,1),1);
D3([1 end],:) = 0;

% Hanning filtering width
hann_width = 5;

% create smoothed amplitude matrix
amp5 = amp;
for i = 1:nt
    for j = 1:nb
        amp5(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
    end
end

% find edges
L0 = 100;
epos5 = nan(nt,nb);
for i = 1:nt
    for j = 1:nb
        % grab a single profile and fill nans
        ampij = amp5(i,:,j)';
        ampij(isnan(ampij)) = ampij(find(~isnan(ampij),1));
        
        if j == 2
            figure(99); clf; hold on
            plot(dist,ampij,'k.-') % debug
        end

        % taper profile around expected position
        if ~any(tdrill==time(i))
            ampij = taper(x0(j),L(j),L0,dist,ampij,'tri'); % runs into trouble with nans in dist?
        else
            idxt = (1:ndrill)*(tdrill==time(i));
            fprintf('drill %d: %s == %s\n',idxt,tdrill(idxt),time(i))
        end
        
        [~,eidx5] = max(D3*ampij);
        epos5(i,j) = dist(eidx5);
        x0(j) = dist(eidx5);
        if j == 2
            plot(dist,ampij,'r.-') % debug
            plot(x0(j)*[1 1],[0 max(ampij)],'b--') % debug
            title(sprintf('%s',time(i)))
        end
    end
    a = 1;
end
edges.pos = epos5;
edges.desc = 'max CD (3pt) with hanning filter';
edges.hann_width = hann_width;

figure(205); clf
plot(time,epos5,'.-')

% % 6. centered-difference d/dz operator
% % derivative operator (5pt stencil, unscaled)
% %D5 = zeros(nc) - 8*diag(ones(nc-1,1),-1) + 8*diag(ones(nc-1,1),1) + diag(ones(nc-2,1),-2) - diag(ones(nc-2,1),2);
% D5 = zeros(nc) - diag(ones(nc-1,1),-1) + diag(ones(nc-1,1),1) - diag(ones(nc-2,1),-2) + diag(ones(nc-2,1),2);
% D5([1 2 end-1 end],:) = 0;
% 
% hann_width = 1;
% amp6 = amp;
% 
% for i = 1:nt
%     for j = 1:nb
%         amp6(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
%     end
% end
% 
% epos6 = nan(nt,nb);
% for i = 1:nt
%     for j = 1:nb
%         ampij = amp6(i,:,j)';
%         ampij(isnan(ampij)) = ampij(find(~isnan(ampij),1));
% 
%         [~,eidx6] = max(D5*ampij);
%         epos6(i,j) = dist(eidx6);
%     end
% end
% out(6).pos = epos6;
% out(6).desc = 'max CD (5pt) with hanning filter';
% out(6).hann_width = hann_width;
% 
% figure(206); clf
% plot(time,epos6)
% 
% % 7. weighted mean sum
% hann_width = 1;
% wms_power = 20;
% amp7 = amp;
% 
% for i = 1:nt
%     for j = 1:nb
%         amp7(i,:,j) = hannFilter(squeeze(amp(i,:,j)),hann_width);
%     end
% end
% 
% epos7 = nan(nt,nb);
% for i = 1:nt
%     for j = 1:nb
%         ampij = amp7(i,:,j)';
%         ampij(isnan(ampij)) = ampij(find(~isnan(ampij),1));
%         
%         maxij = max(ampij);
%         ampij = ampij/maxij;
%         
%         sumij = sum(ampij.^wms_power);
%         wampij = (ampij.^wms_power).*(1:length(ampij))';
%         
%         eidx7 = round(sum(wampij)/sumij);
%         epos7(i,j) = dist(eidx7);
%     end
% end
% out(7).pos = epos7;
% out(7).desc = sprintf('WMS (p=%d) with hanning filter',wms_power);
% out(7).hann_width = hann_width;
% 
% figure(207); clf
% plot(time,epos7)

% for i = 1:length(edges)
%     edges(i).time = time;
% end
edges.time = time;

% plot over .pck profiles
figure(299); clf
clear ax
for i = 1:3
    ax(i) = axes(figure(299),'position',axgridpos(1,3,i,.05,.05,.05,.08));
    hold on
    pcolor(time,dist,amp(:,:,i)')
    shading interp
    clim(extrema(amp(:,dist>50,:)))
    datetick('x','keeplimits')
    
    plot(time,edges.pos(:,i),'k.-')
%     for j = 5%1:length(out)
%         plot(time,edges(j).pos(:,i),'k.-')%,color=colors(j))
%     end

    linkaxes(ax)
    xlim(ax(1),extrema(time))
    ylim(ax(1),extrema(dist))
end

%%% subfunctions %%%
function prof = taper(x0,L,L0,dist,amp,shape)
    if strcmp(shape,'tri')
        w = L/(L-L0) - 2/(L-L0)*abs(dist'-x0);
        w(w<0) = 0;
        w(w>1) = 1;
        prof = w.*amp;
    elseif strcmp(shape,'hann')
        a = 1;
    end


