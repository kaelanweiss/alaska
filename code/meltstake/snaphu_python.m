% script to figure out 1) how to pass data between matlab and python in WSL
% and 2) how to best use the SNAPHU package

clear

raw_dir = 'F:/meltstake/data/raw';
unw_dir = 'F:/meltstake/data/unwrap';

dep = 'ms02_20240712_2020';
% dep = 'ms02_20240715_2001';
% dep = 'ms01_20240717_0050';

if ~exist(fullfile(unw_dir,dep),'dir')
    mkdir(fullfile(unw_dir,dep))
end

load(fullfile(raw_dir,dep,'adcp','adcp.mat'))

ax_bm = adcpQuickPlot(figure(1),adcp,'vel',extrema(adcp.burst.vel),NaT,[0 1],1);
ax_unw = adcpQuickPlot(figure(2),adcp,'vel_unw',extrema(adcp.burst.vel),NaT,[0 1],1);

%% maximum size of tile (machine limitation)
L_MAX = 32000;
% L_MAX = 10000;

% calculate width of tiles
ovrlp = 128;
nt = size(adcp.burst.vel,1);
n_tiles = ceil(nt/(L_MAX-ovrlp));
t_incr = ceil(nt/n_tiles);
% tile indices
tile_idx = [1:t_incr:n_tiles*t_incr;...
            (t_incr:t_incr:n_tiles*t_incr)+ovrlp];
tile_idx(end,end) = nt;
% save tiles
for i = 1:n_tiles
    k1 = tile_idx(1,i);
    k2 = tile_idx(2,i);
    vel = adcp.burst.vel(k1:k2,:,:);
    cor = adcp.burst.cor(k1:k2,:,:);
    indices = (k1:k2)';
    save(fullfile(unw_dir,dep,sprintf('velcor%02d.mat',i)),'vel','cor','indices')
end

% SNAPHU happens in python

%% load tiles
clear snaphu
snaphu(n_tiles) = struct();
for i = 1:n_tiles
    load_struct = load(fullfile(unw_dir,dep,sprintf('vel_unw%02d.mat',i)));
    snaphu(i).phi = load_struct.phi_unw;
    snaphu(i).v_max = load_struct.v_max;
    snaphu(i).cc = load_struct.conn_comp;
    snaphu(i).indices = cast((load_struct.indices(1):load_struct.indices(2))','uint32');
end

%% stitch together phase tiles
phi_unw = nan(size(adcp.burst.vel));
phi_unw(snaphu(1).indices,:,:) = snaphu(1).phi;
snaphu(1).used_indices = extrema(snaphu(1).indices);

% find overlap phase offsets
phi_offsets = zeros(5,n_tiles);
for i = 1:n_tiles-1 % loop through overlap sections
    ovrlp_range = extrema(intersect(snaphu(i).indices,snaphu(i+1).indices));
    idx1 = snaphu(i).indices >= ovrlp_range(1) & snaphu(i).indices <= ovrlp_range(2);
    idx2 = snaphu(i+1).indices >= ovrlp_range(1) & snaphu(i+1).indices <= ovrlp_range(2);
    for j = 1:5 % loop through beams
        phi1 = snaphu(i).phi(idx1,:,j);
        phi2 = snaphu(i+1).phi(idx2,:,j);
        phi_offsets(j,i+1) = median(round(phi2-phi1,3),'all');
        snaphu(i+1).used_indices = snaphu(i+1).indices([find(~idx2,1) end]);
        phi_unw(snaphu(i).indices(end)+1:snaphu(i+1).indices(end),:,j) = snaphu(i+1).phi(find(~idx2,1):end,:,j) -  sum(phi_offsets(j,1:i+1));
    end
end

phi_med = squeeze(median(phi_unw,2,'omitnan'));
figure(2); clf
plot(phi_med,'.-')

%% plot some panels
% phi_offsets = zeros(5,n_tiles);
% for i = 1:5
%     for j = 1:n_tiles
%         phi_offsets(i,j) = 0; %round(mode(snaphu(j).phi(:,:,i),'all'));
%     end
% end

clear ax
for i = 1:5
    figure(10+i); clf
    CLIM = [0 1];
    for j = 1:n_tiles
        ax(i,j) = subplot(n_tiles,1,j);
        pcolor(ax(i,j),snaphu(j).indices,adcp.burst.range,snaphu(j).phi(:,:,i)' - sum(phi_offsets(i,1:j)))
        shading flat
        colorbar
        cmocean('curl')

        CLIM = extrema([extrema(snaphu(j).phi(:,:,i));CLIM]);
    end
%     set(ax(i,:),'CLim',CLIM)
end
linkaxes(ax)
set(ax,'CLim',1.5*[-1 1])

%%
clear ax_all
figure(21); clf
for i = 1:5
    ax_all(i) = subplot(5,1,i);
    pcolor(ax_all(i),1:nt,adcp.burst.range,phi_unw(:,:,i)')
    shading flat
    colorbar
    cmocean('curl')
end
linkaxes(ax_all)
set(ax_all,'CLim',1.5*[-1 1])

% perform manual steps to shift phase offsets within beams

%% rescale velocities
vel_unw = nan(size(phi_unw));
for i = 1:n_tiles
    k1 = snaphu(i).used_indices(1);
    k2 = snaphu(i).used_indices(2);
    for j = 1:5
        vel_unw(k1:k2,:,j) = 2*snaphu(i).v_max(j)*phi_unw(k1:k2,:,j);
    end
end

adcp.burst.vel_unw = vel_unw;
ax_unw = adcpQuickPlot(figure(22),adcp,'vel_unw',1*extrema(adcp.burst.vel),NaT,[0 1],1);

%% save unwrapped velocity when ready
save(fullfile(unw_dir,dep,'vel_unw.mat'),'vel_unw','phi_unw')