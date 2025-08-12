% Script to plot spectra of boundary layer forcing
%
% KJW
% 11 Apr 2025

clear

addpath('..')

raw_dir = 'F:/meltstake/data/raw';

% ms deployment info
ms_tbl = loadMSInfo(26:28,'manualwindows');
[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);
dep_depth = ms_tbl.depth(uidx);

% load data
adcp = cell(ndeps,1);
T = cell(ndeps,1);
tlims = cell(ndeps,1);
rlims = nan(ndeps,1);
for i = 1:ndeps
    % adcp
    load_data = load(fullfile(raw_dir,dep_names{i},'adcp','adcp.mat'));
    adcp{i} = load_data.adcp;
    tlims{i} = extrema([ms_tbl.Start(ms_tbl.Number==dep_nums(i)) ms_tbl.End(ms_tbl.Number==dep_nums(i))]);
    rlims(i) = min(ms_tbl.rmax(ms_tbl.Number==dep_nums(i)));
    % T
    load_data = load(fullfile(raw_dir,dep_names{i},'rbr','T.mat'));
    T{i} = load_data.T;
end
adcp = cell2mat(adcp);
T = cell2mat(T);

% drilling times
dt_blank = seconds(10);
tdrill = cell(ndeps,1);
tdrill{1} = datetime(2024,7,12,20,41,25) + minutes(20)*(0:3)';
tdrill{2} = datetime(2024,7,15,20,21,6) + minutes(20)*(0:7)';
tdrill{3} = datetime(2024,7,17,1,11,53) + minutes(20)*(0:6)';

%% calculate mean component timeseries
clear psd T_psd
psd(ndeps) = struct;
T_psd(ndeps) = struct;
clear ax4
figure(4); clf
ax4 = axes; hold on
for i = 1:ndeps
    % collect relevant data from adcp structure
    idxt = adcp(i).burst.time >= tlims{i}(1) & adcp(i).burst.time <= tlims{i}(2);
    idxr = adcp(i).burst.range <= rlims(i);
    psd(i).fs = adcp(i).burst.samplerate;
    psd(i).time = adcp(i).burst.time(idxt);
    % mean components
%     idxr = round(find(idxr,1,'last')/2); %%% CHANGED FOR SINLGE POINT %%%
    psd(i).uvw = squeeze(mean(adcp(i).burst.vel_ice(idxt,idxr,[1 3 2]),2,'omitnan'));
%     psd(i).uvw = squeeze(mean(adcp(i).burst.vel(idxt,idxr,[3 4 5]),2,'omitnan'));
%     plot(ax4,psd(i).time,psd(i).uvw,'k-')

    % blank out drill times
    idx_drill = false(size(psd(i).time));
    for j = 1:length(tdrill{i})
        idxj = psd(i).time >= tdrill{i}(j) & psd(i).time <= tdrill{i}(j)+dt_blank;
        idx_drill = idx_drill | idxj;
    end
    psd(i).uvw(idx_drill,:) = nan;
%     plot(ax4,psd(i).time,psd(i).uvw)
    % find nans
    idx_nan = isnan(psd(i).uvw(:,1));
    % fill with linear interpolation
    uvw_fill = psd(i).uvw;
    idx_fill = idx_drill | idx_nan;
    for j = 1:3
        uvw_fill(idx_fill,j) = interp1(psd(i).time(~idx_fill),psd(i).uvw(~idx_fill,j),psd(i).time(idx_fill),'linear');
        uvw_fill(:,j) = hannFilter(uvw_fill(:,j),5);
%         plot(ax4,psd(i).time(idx_drill),uvw_fill(idx_drill,j),'r.-')
    end
    psd(i).uvw(idx_fill,:) = uvw_fill(idx_fill,:);

    % T
    T_psd(i).time = T(i,1).time;
    T_psd(i).fs = seconds(diff(T(i,1).time([1 2])))^-1;
    T_psd(i).T = [T(i,1).values T(i,2).values T(i,3).values];
end

%% calculate velocity PSDs
% % confidence interval
% r = 1; %1.96; % Hanning window reduction of DoFs
% alpha = 0.05;
% M = 2*nf_bin/r;
% P_ci = M./chi2inv([1-alpha/2 alpha/2],M);

% binning
nf_u = [4 15 40 60]';
nf_u_stops = [0 80 401 1001]'/3600;
% confidence intervals
alpha = 0.05;

for i = 1:ndeps
    % raw PSD
    [psd(i).Praw,psd(i).fraw] = psd_matlab(psd(i).uvw,psd(i).fs,'notaper');
    psd(i).dfraw = diff(psd(i).fraw([1 2]));
    % binned PSD
    [Pbin,fbin,Mbin] = progressiveBin(psd(i).fraw(2:end),psd(i).Praw(2:end,:),nf_u,nf_u_stops);
    psd(i).P = Pbin;
    psd(i).df = diff(fbin([1 2]));
    psd(i).f = fbin;
    psd(i).M = Mbin;
    psd(i).s_ci = [Mbin./chi2inv(1-alpha/2,Mbin) Mbin./chi2inv(alpha/2,Mbin)];

%     % plot
%     figure(20+i); clf
%     hold on
%     for j = 1:3
% %         plot(psd(i).fraw*3600,psd(i).Praw(:,j),'color',clrs(j,:))
%         plot(psd(i).f*3600,psd(i).P(:,j),'linewidth',1)
%     end
% %     plot(psd(i).f*3600,psd(i).P,'k-','linewidth',1.5)
%     set(gca,'xscale','log','yscale','log')
end

%% temperature PSDs
% binning
nf_T = [4 10 20 60]';
nf_T_stops = [0 80 400 1001]'/3600;
% confidence intervals
alpha = 0.05;

for i = 1:ndeps
    % raw PSD
    T_in = T_psd(i).T - repmat(mean(T_psd(i).T),[size(T_psd(i).T,1) 1]);
    [T_psd(i).Praw,T_psd(i).fraw] = psd_matlab(T_in,T_psd(i).fs);
    T_psd(i).dfraw = diff(T_psd(i).fraw([1 2]));
    % binned PSD
    [Pbin,fbin,Mbin] = progressiveBin(T_psd(i).fraw(2:end),T_psd(i).Praw(2:end,:),nf_T,nf_T_stops);
    T_psd(i).P = Pbin;
    T_psd(i).df = diff(fbin([1 2]));
    T_psd(i).f = fbin;
    T_psd(i).M = Mbin;
    T_psd(i).s_ci = [Mbin./chi2inv(1-alpha/2,Mbin) Mbin./chi2inv(alpha/2,Mbin)];
    
%     % plot
%     figure(30+i); clf
%     hold on
%     for j = 1:3
% %         plot(psd(i).fraw*3600,psd(i).Praw(:,j),'color',clrs(j,:))
%         plot(T_psd(i).f*3600,T_psd(i).P(:,j),'linewidth',1)
%     end
% %     plot(psd(i).f*3600,psd(i).P,'k-','linewidth',1.5)
%     set(gca,'xscale','log','yscale','log')
end

%% plot
fs = 11;
lw = 1.5;
pad = [.03 .05 .05 .12];
shift = [0.03 0.04];
figsize = [20 10];
fig1 = setFigureSize(figure(1),figsize);
clf(fig1)
fig2 = setFigureSize(figure(2),figsize);
clf(fig2)

clrs = [237 248 177; 127 205 187; 44 127 184]/255;
clrs = [colors(5); colors(1); colors(7)];
clrs = hexTriplet({'0c2c84','1d91c0','7fcdbb'});
clrs(3,:) = clrs(3,:)/1.5;
% clrs = flipud(clrs);

T_clr = [0.2 1 0.5]/1.1;
T_clrs = (T_clr'*[0.95 0.5 0.2])';

N_ms = [.002 .0036;...
        .0014 .004;...
        .002 .0032]; % cycles/sec

f_scale = 1;%60*60;
vp = 0;

clear ax1 ax2 hu hT
for i = 1:ndeps
    % velocity
    ax1(i) = axes(fig1,'position',axgridpos(1,ndeps,i,pad,shift),'fontsize',fs);
    hold(ax1(i),'on')
    box(ax1(i),'on')
    set(ax1(i),'xscale','log','yscale','log')

    % local N
%     plot(ax1(i),N_ms(i,1)*[1 1]*f_scale,[1e-7 1e3],'color',0.5*[1 1 1],'linewidth',1.5)
    patch(ax1(i),[N_ms(i,:) flip(N_ms(i,:))],[1e-7 1e-7 1e3 1e3],0.8*[1 1 1],'facealpha',1,'edgecolor','none');
    text(ax1(i),1.1*N_ms(i,2)*f_scale,1e-4/(10^vp),'N_{local}','Rotation',0,'fontsize',fs,'horizontalalignment','left')

    % shaded regions
    for j = 1:3
        plotPConf(ax1(i),psd(i).f*f_scale,psd(i).P(:,j).*psd(i).s_ci(:,1).*psd(i).f.^vp,psd(i).P(:,j).*psd(i).s_ci(:,2).*psd(i).f.^vp,clrs(j,:),0.3);
    end

    % psds
    for j = [3 1 2]
        hu(j) = plot(ax1(i),psd(i).f*f_scale,psd(i).P(:,j).*psd(i).f.^vp,'-','color',clrs(j,:),'linewidth',lw);
%         plot(ax1(i),psd(i).fraw(nf_u(1)-1)*f_scale,psd(i).Praw(1,j).*psd(i).fraw(nf_u(1)-1).^vp,'o','color',clrs(j,:),'markerfacecolor',clrs(j,:),'markersize',6)
    end
%     plot(ax1(i),N_ms(i)*[1 1]*f_scale,[1e-7 1e3],'color',0.5*[1 1 1],'linewidth',1.5)
%     text(ax1(i),1.1*N_ms(i)*f_scale,1e-4/(10^vp),'N_{local}','Rotation',0,'fontsize',fs,'horizontalalignment','left')
    plot(ax1(i),10.^[-1 0.8],(6e-4)*(10.^([-1 0.8]*(-5/3))),'k--','linewidth',lw)
    text(ax1(i),1,3e-3,'f^{-5/3}','fontsize',fs)

    % temperature
    ax2(i) = axes(fig2,'position',axgridpos(1,ndeps,i,pad,shift),'fontsize',fs);
    hold(ax2(i),'on')
    box(ax2(i),'on')
    set(ax2(i),'xscale','log','yscale','log')

    patch(ax2(i),[N_ms(i,:) flip(N_ms(i,:))],[1e-7 1e-7 1e3 1e3],0.8*[1 1 1],'facealpha',1,'edgecolor','none');

    % shaded regions
    for j = 1:3
        plotPConf(ax2(i),T_psd(i).f*f_scale,T_psd(i).P(:,j).*T_psd(i).s_ci(:,1).*T_psd(i).f.^vp,T_psd(i).P(:,j).*T_psd(i).s_ci(:,2).*T_psd(i).f.^vp,T_clrs(j,:),0.3);
    end

    % psds
    for j = [3 2 1]
        hT(j) = plot(ax2(i),T_psd(i).f*f_scale,T_psd(i).P(:,j).*T_psd(i).f.^vp,'-','color',T_clrs(j,:),'linewidth',lw);
    end
%     patch(ax2(i),[N_ms(i,:) flip(N_ms(i,:))],[1e-7 1e-7 1e3 1e3],0.8*[1 1 1],'facealpha',1,'edgecolor','none');
%     plot(ax2(i),N_ms(i)*[1 1]*f_scale,[1e-7 1e3],'color',0.5*[1 1 1],'linewidth',1.5)
    text(ax2(i),1.1*N_ms(i,2)*f_scale,5e-2/(100^vp),'N_{local}','Rotation',0,'fontsize',fs,'horizontalalignment','left')
    plot(ax2(i),10.^[-2 -0.],5e-3*(10.^([-2 -0.]*(-5/3))),'k--','linewidth',lw)
    text(ax2(i),1e-1,5e-1,'f^{-5/3}','fontsize',fs)

    a = 5;

end

% axis limits
linkaxes(ax1)
xlim(ax1(1),[8e-1 2e4]*f_scale/3600)
ylim(ax1(1),[8e-7*(4^vp) 2e1/(100^vp)])

linkaxes(ax2)
xlim(ax2(1),[8e-1 5e3]*f_scale/3600)
ylim(ax2(1),[1e-4*(4^vp) 2e2/(100.^vp)])

% axis labels
for i = 1:3
    set(ax1(i),'xtick',(10.^(-10:5)));
    set(ax2(i),'xtick',(10.^(-10:5)));
    set(ax1(i),'ytick',(10.^(-10:2:5)));
    set(ax2(i),'ytick',(10.^(-10:5)));
%     set(ax1(i),'xtick',(10.^(-1:5))*f_scale/3600);
%     set(ax2(i),'xtick',(10.^(-1:5))*f_scale/3600);
    xlabel(ax1(i),'f [cps]','fontsize',fs)
    xlabel(ax2(i),'f [cps]','fontsize',fs)
    title(ax1(i),sprintf('deployment %d',i),'fontsize',fs)
    title(ax2(i),sprintf('deployment %d',i),'fontsize',fs)
end
for i = 2:3
    set(ax1(i),'yticklabels',[])
    set(ax2(i),'yticklabels',[])
end
ylabel(ax1(1),'PSD [m^2/s^2 cps^{-1}]','fontsize',fs)
ylabel(ax2(1),'PSD [K^2 cps^{-1}]','fontsize',fs)

% legends
lgdu = legend(ax1(1),hu,{'u','v','w'},'fontsize',fs,'location','northeast');
legend(ax2(1),hT,{'0.1 m','0.2 m','0.5 m'},'fontsize',fs,'location','northeast');

%%% subfunctions %%%
function hp = plotPConf(ax,x,y_l,y_h,clr,alpha)
% function to create shaded confidence region, inputs should be column
% vectors
    x_patch = [x; flip(x)];
    y_patch = [y_l; flip(y_h)];
    hp = patch(ax,x_patch,y_patch,clr,'facealpha',alpha,'edgecolor','none');
end

function rgb = hexTriplet(hex_clrs)
% function to convert hex codes to scaled rgb triplets
% hex_clrs should be a cell array of 6-digit hex codes
    nr = length(hex_clrs);
    rgb = nan(nr,3);
    for i = 1:nr
        for j = 1:3
            rgb(i,j) = hex2dec(hex_clrs{i}([1 2]+2*(j-1)))/256;
        end
    end
end