% Script to play sandbox edge-finding routine in Spider ADCP data (and
% possibly beyond)
%
% KJW
% 4/24/23
%
% This is dialed for SF24 just using velocity data.
clear
try
    adcp;
catch
%     load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
%     load F:/alaska2022/data/iceberg_surveys/mat/20220824_singingflower/spider/adcp_exclude.mat
    load F:/AK_202305/adcp/103041_AK1/adcp.mat
end

% QA
% manual exlusions
% sz_vel = size(adcp.burst.vel);
% sz_vel(1) = 1;
% exclude_idx = repmat(exclude_idx,sz_vel);
% adcp.burst.vel(exclude_idx) = nan;
% correlation filter
qa_cor = adcp.burst.cor <= 50;
adcp.burst.vel(qa_cor) = nan;

% plot data
tlim = [0 Inf];
rlim = [0 1];
flds = {'vel','amp','cor'};
clims = {0.02*[-1 1],[70 85],[50 100]};
dt = 2;

ax = struct();
for i = 1:3
    fnum = i;
    figure(fnum); clf;
    ax.(flds{i}) = adcpQuickPlot(figure(fnum),adcp,flds{i},clims{i},tlim,rlim,dt);
end

%% qa (check how the routine does without qa first)
khann = 2*8+1;
for k = 1:3
    for j = 1:size(adcp.burst.(flds{k}),3)
        for i = 1:size(adcp.burst.(flds{k}),2)
            adcp.burst.(flds{k})(:,i,j) = hannFilter(adcp.burst.(flds{k})(:,i,j),khann);
        end
    end
end

%% findBlocks

vmax = 0.001;
cmin = 98;
amin = 84;
mask = true([size(adcp.burst.vel) 3]);

mask(:,:,:,1) = abs(adcp.burst.vel)<=vmax;
mask(:,:,:,2) = adcp.burst.cor>=cmin;
mask(:,:,:,3) = adcp.burst.amp>=amin;

clear ax_bw
bn = 2;
figure(5); clf
for i = 1:3
    ax_bw(i) = axes(figure(5),'position',axgridpos(3,1,i,0.05,0.05,0.1,0.08));
    pcolor(mask(:,:,bn,i)')
    shading flat
    colormap bone
end

pxlim = [5 20];
mask = mask(:,pxlim(1):pxlim(2),:,:);

%%
edge_idx = nan(size(squeeze(mask(:,1,:,:))));
for k = 1:size(mask,4) % fields 
    for j = 1:size(mask,3) % beams
        for i = 1:size(mask,1) % time
            % find largest block of true values
            [blcks,wdths] = findBlocks(mask(i,:,j,k));
            if ~isempty(blcks)
                [~,imax] = max(wdths);
                edge_idx(i,j,k) = pxlim(1)+blcks(imax,1)-1;
            end
        end
    end
end
edge_idx(:,:,1) = edge_idx(:,:,1);

kmed = 8*30+1;
edge_idx_med = medfilt1(edge_idx,kmed,'omitnan','truncate');
edge_idx_all = ceil(medfilt1(mean(edge_idx_med(:,:,1),3),2*kmed,'omitnan','truncate'));

edge_idx_all(isnan(edge_idx_all)) = 1;%length(adcp.burst.range);
ice_pos = adcp.burst.range(edge_idx_all);

%%

for i = 1:3
    hold(ax_bw(i),'on')
    %plot(ax_bw(i),edge_idx(:,bn,i),'r')
    plot(ax_bw(i),edge_idx_all(:,bn),'r','linewidth',1)

    for j = 1:size(adcp.burst.vel,3)
        hold(ax.(flds{i})(j),'on')
        plot(ax.(flds{i})(j),adcp.burst.time,ice_pos(:,j),'r')
    end
end


