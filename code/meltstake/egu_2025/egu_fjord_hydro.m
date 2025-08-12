% Script to plot open-fjord context during meltstake deployments
%
% KJW
% 9 Apr 2025

clear

addpath('..')

proc_dir = 'F:/meltstake/data/proc';

% ms deployment info
ms_tbl = loadMSInfo(26:28,'manualwindows');
[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);
dep_depth = ms_tbl.depth(uidx);

% load ctd data
lgd_lbls = cell(ndeps,1);
ctd = cell(ndeps,1);
for i = 1:ndeps
    load_data = load(fullfile(proc_dir,dep_names{i},'ctd.mat'));
    ctd{i} = load_data.ctd;
    % use polly data if available
    npp = size(ctd{i}(strcmp('polly',{ctd{i}.vessel})).T,2);
    if npp
        ctd{i} = ctd{i}(strcmp('polly',{ctd{i}.vessel}));
    else
        ctd{i} = ctd{i}(~strcmp('polly',{ctd{i}.vessel}));
    end
    t_start = ms_tbl.Start(uidx(i));
    t_start.Format = 'MMMdd yyyy';
    lgd_lbls{i} = sprintf('dep %d (%s)',i,t_start); 
end 
ctd = cell2mat(ctd);

%% data processing
for i = 1:ndeps
    [nz,np] = size(ctd(i).T);
    % clean up ground-strikes in salinity
    iz_min = 0;
    for j = 1:np
        iz_bottom = nz - find(~isnan(flip(ctd(i).S(:,j))),1) + 1;
        [~,idx_filt] = sigmaFilter(ctd(i).S(iz_bottom-10:iz_bottom,j),1,1,2);
        if any(idx_filt(end-1:end))
            % if any outliers detected, clear bottom meter
            ctd(i).S(iz_bottom-3:iz_bottom,j) = ctd(i).S(iz_bottom-4,j);
            ctd(i).PD(iz_bottom-3:iz_bottom,j) = ctd(i).PD(iz_bottom-4,j);
        end
        iz_min = max([iz_min iz_bottom]);
    end

    % calculate mean profiles
    ctd(i).depth_mean = ctd(i).depth(1:iz_min);
    ctd(i).T_mean = mean(ctd(i).T(1:iz_min,:),2,'omitnan');
    ctd(i).S_mean = mean(ctd(i).S(1:iz_min,:),2,'omitnan');
    ctd(i).rho_mean = mean(ctd(i).PD(1:iz_min,:),2,'omitnan')+1000;

    ctd(i).T_std = std(ctd(i).T(1:iz_min,:),0,2,'omitnan');
    ctd(i).T_min = min(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    ctd(i).T_max = max(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    ctd(i).S_std = std(ctd(i).S(1:iz_min,:),0,2,'omitnan');
    ctd(i).S_min = min(ctd(i).S(1:iz_min,:),[],2,'omitnan');
    ctd(i).S_max = max(ctd(i).S(1:iz_min,:),[],2,'omitnan');
    ctd(i).rho_std = std(ctd(i).PD(1:iz_min,:),0,2,'omitnan');
    ctd(i).rho_min = min(ctd(i).PD(1:iz_min,:),[],2,'omitnan')+1000;
    ctd(i).rho_max = max(ctd(i).PD(1:iz_min,:),[],2,'omitnan')+1000;

    % calculate N2
    dz = diff(ctd(i).depth(1:2));
    rho_mean_smth = hannFilter(ctd(i).rho_mean,round(3/dz)+1);
    drho_dz = gradient(rho_mean_smth)/dz;
    ctd(i).N2 = 9.81*drho_dz./rho_mean_smth/(2*pi)^2;
    ctd(i).N = ctd(i).N2;
    ctd(i).N(ctd(i).N<0) = nan;
    ctd(i).N = sqrt(ctd(i).N);

    % find melstake depth index
    [~,ctd(i).zidx] = min(abs(ctd(i).depth_mean-dep_depth(i)));

    % find range of N values
    z_range = 5; % m
    dzi = round(z_range/dz);
    ctd(i).N_ms = extrema(ctd(i).N((ctd(i).zidx-dzi):(ctd(i).zidx+dzi)));

end

%% plot
fs = 11;
ms = 10;
pad = [.03 .05 .1 .1];
shift = [0 -.05];
figsize = [22 12];
fig = setFigureSize(figure(1),figsize);
clf
plt_lbls = {'T','S','\rho','N'};
plt_units = {'\circC','psu','kg/m^3','cps'};
N_scale = 1;
load clrs_egu.mat
clrs(2,:) = 0.9*clrs(2,:);

clear ax
for i = 4:-1:1
    ax(i) = axes(fig,'position',axgridpos(1,4,i,pad,shift),'fontsize',fs);
    hold on
    box on
    axis ij
    title(ax(i),[plt_lbls{i} ' [' plt_units{i} ']'],'fontsize',fs+1)
    ax(i).XAxisLocation = 'top';
%     xlabel(ax(i),['[' plt_units{i} ']'],'fontsize',fs)
end

% profiles
f_alpha = 0.4;
for i = 1:3
    % shaded confidence intervals
    plotProfConf(ax(1),ctd(i).depth_mean,ctd(i).T_min,ctd(i).T_max,clrs(i,:),f_alpha)
    plotProfConf(ax(2),ctd(i).depth_mean,ctd(i).S_min,ctd(i).S_max,clrs(i,:),f_alpha)
    plotProfConf(ax(3),ctd(i).depth_mean,ctd(i).rho_min,ctd(i).rho_max,clrs(i,:),f_alpha)
%     plotProfConf(ax(2),ctd(i).depth_mean,ctd(i).S_mean,ctd(i).S_std,clrs(i,:),f_alpha)
%     plotProfConf(ax(3),ctd(i).depth_mean,ctd(i).rho_mean,ctd(i).rho_std,clrs(i,:),f_alpha)
end
for i = 1:3
    % T
    plot(ax(1),ctd(i).T_mean,ctd(i).depth_mean,'color',clrs(i,:),'linewidth',2)
    % S
    plot(ax(2),ctd(i).S_mean,ctd(i).depth_mean,'color',clrs(i,:),'linewidth',2)
    % rho
    plot(ax(3),ctd(i).rho_mean,ctd(i).depth_mean,'color',clrs(i,:),'linewidth',2)
    % N2
    plot(ax(4),ctd(i).N*N_scale,ctd(i).depth_mean,'color',clrs(i,:),'linewidth',2)
end

% ms depths
clear h
for i = 1:3
    % T
    h(i) = plot(ax(1),ctd(i).T_mean(ctd(i).zidx),ctd(i).depth(ctd(i).zidx),['k' smbs{i} '-'],'color',clrs(i,:),'markerfacecolor',clrs(i,:),'markeredgecolor',0*clrs(i,:),'linewidth',1.5,'markersize',ms);
    h(4) = plot(ax(1),msTable2Vector(ms_tbl.T(ms_tbl.Number==dep_nums(i))),ctd(i).depth(ctd(i).zidx)*ones(sum(ms_tbl.Number==dep_nums(i)),1),['k' smbs{i}],'markerfacecolor',0.5*[1 1 1],'linewidth',1,'markersize',ms-6);
    % S
    plot(ax(2),ctd(i).S_mean(ctd(i).zidx),ctd(i).depth(ctd(i).zidx),['k' smbs{i}],'color',clrs(i,:),'markerfacecolor',clrs(i,:),'markeredgecolor',0*clrs(i,:),'linewidth',1.5,'markersize',ms)
    % rho
    plot(ax(3),ctd(i).rho_mean(ctd(i).zidx),ctd(i).depth(ctd(i).zidx),['k' smbs{i}],'color',clrs(i,:),'markerfacecolor',clrs(i,:),'markeredgecolor',0*clrs(i,:),'linewidth',1.5,'markersize',ms)
    % N2
    plot(ax(4),ctd(i).N(ctd(i).zidx)*N_scale,ctd(i).depth(ctd(i).zidx),['k' smbs{i}],'color',clrs(i,:),'markerfacecolor',clrs(i,:),'markeredgecolor',0*clrs(i,:),'linewidth',1.5,'markersize',ms)
end

% x limits
xlim(ax(1),[4.5 8])
xlim(ax(2),[24 31.5])
xlim(ax(3),[1018 1025])
xlim(ax(4),[0 0.006]*N_scale)
linkaxes(ax,'y')

% y labels
for i = 2:4
    set(ax(i),'yticklabels',[])
end
ylabel(ax(1),'depth [m]','fontsize',fs)

% legend
lgd = legend(ax(1),h,cat(1,lgd_lbls,{'onboard T'}),'location','southeast','fontsize',fs);
% lgd.Title.String = 'Meltstakes';
lgd.Title.FontSize = fs;
lgd.Position = lgd.Position + [.15 0 0 0];


%%% subfunctions %%%
function hp = plotProfConf(ax,x,y_l,y_h,clr,alpha)
    % function to create shaded confidence region, inputs should be column
    % vectors
    x_patch = [x; flip(x)];
    y_patch = [y_l; flip(y_h)];
    hp = patch(ax,y_patch,x_patch,clr,'facealpha',alpha,'edgecolor','none');
end