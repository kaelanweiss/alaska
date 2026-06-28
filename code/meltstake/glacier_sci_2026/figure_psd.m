% Script to plot PSD of ADCP velocity.
%
% KJW
% 27 Jun 2026
clear


load ..\glacier_grl_2025\platform_comparison_data\ms_psd.mat
load glacier_clrs.mat
load F:/meltstake/data/proc/glacier_2024_ctd.mat

dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

ms_depth = nan(3,1);
for i = 1:3
    ms_depth(i) = mean(seg_tbl.depth(seg_tbl.Number==dep_tbl.Number(i)));
end

%% prep adcp data
% where to sample for PSDs and velocities
y_s = [0.18] + 0.05*[-1 1]; % distance from wall 
for i = 1:length(adcp)
    adcp(i).range = round(double(adcp(i).range),4);
    y_max = adcp(i).range(end);
%     [~,idx_s] = min(abs(y_max-adcp(i).range-y_s));
    idx_s = (y_max-adcp(i).range) >= min(y_s) & (y_max-adcp(i).range) <= max(y_s);
    adcp(i).idx_s = idx_s;
    adcp(i).uvw = squeeze(mean(adcp(i).vel_ice(:,idx_s,:),2));
    adcp(i).Puvw = squeeze(mean(adcp(i).P_ice(:,idx_s,:),2));
end
if ~(sum(adcp(1).idx_s) == sum(adcp(2).idx_s) && sum(adcp(1).idx_s) == sum(adcp(3).idx_s))
    warning('Number of bins do not agree')
end

% trim out high frequency data
f_max = 2;
for i = 1:length(adcp)
    idx_trim = adcp(i).f > f_max;
    adcp(i).f(idx_trim) = [];
    adcp(i).Puvw(idx_trim,:) = [];
end

% bin PSDs
% binning setup
nf = [4 8 16 16];
nf_bins = [0 1e-2 1e-1 1e0];

for i = 1:length(adcp)
    [Pb,fb,M,S_ci] = progressiveBin(adcp(i).f,adcp(i).Puvw(:,1),nf,nf_bins);
    adcp(i).fbin = fb;
    adcp(i).Pbin = repmat(Pb,[1 3]);
    adcp(i).Pbin_beam = repmat(Pb,[1 5]);
    adcp(i).Mbin = M;
    adcp(i).S_ci = S_ci;
    for j = 2:3
        adcp(i).Pbin(:,j) = progressiveBin(adcp(i).f,adcp(i).Puvw(:,j),nf,nf_bins);
    end
    for j = 1:5
        adcp(i).Pbin_beam(:,j) = progressiveBin(adcp(i).f,mean(adcp(i).P(:,adcp(i).idx_s,j),2),nf,nf_bins);
    end
    % adcp(i).Pbin(:,2) = adcp(i).Pbin(:,2)/sind(25); % test
end

% CPDFs of BL energy
for i = 1:length(adcp)
    fbi = adcp(i).fbin;
    nfb = length(fbi);
    
    % straight CPDF
    % adcp(i).CPDF = cumtrapz(fbi,sum(adcp(i).Pbin,2)); % all components
    adcp(i).CPDF = cumtrapz(fbi,sum(adcp(i).Pbin(:,[1]),2)); % just u (to minimize noise effects)
    adcp(i).f_80 = adcp(i).fbin(find(adcp(i).CPDF./adcp(i).CPDF(end)>=0.8,1));

end

%% figure
% set up plot
dep = 3;
lw = 1;
fs = 11;
alpha = 0.5;

pad = [0 0 .18 .12];
shift = [0.01 0.05];

fig_size = [8.5 8];
fig = figure(10); clf
setFigureSize(fig,fig_size);
ax = axes(fig,'position',axgridpos(1,1,1,pad,shift));
hold(ax,'on')
box(ax,'on')

% velocity PSDs
clear h_psd

% vertical frequencies
plot(mean(ctd(dep).N_ms)*[1 1],10.^[-6 4],'-','linewidth',1*lw,'color',0.5*[1 1 1])
% plot(adcp(dep).f_80*[1 1],10.^[-6 4],':','linewidth',0.9*lw,'color',0.3*[1 1 1])
Ntxt = text(mean(ctd(dep).N_ms),3e-4,'N/(2\pi)','fontsize',fs);

% shaded regions
for j = [3 1]
    plotConfidenceRegion(ax,adcp(dep).fbin,adcp(dep).S_ci(:,1).*adcp(dep).Pbin(:,j),adcp(dep).S_ci(:,2).*adcp(dep).Pbin(:,j),vel_clrs(j,:),alpha);
end
% P_dep = sum(adcp(dep).Pbin(:,[1 3]),2);
% plotConfidenceRegion(ax,adcp(dep).fbin,adcp(dep).S_ci(:,1).*P_dep,adcp(dep).S_ci(:,2).*P_dep,vel_clrs(1,:),alpha);

% CPDFs
yyaxis(ax,'right')
h_psd(2) = plot(ax,adcp(dep).fbin,adcp(dep).CPDF/adcp(dep).CPDF(end),'-','color',0.*[1 1 1],'linewidth',1.2*lw);
set(ax,'ytick',0:0.25:1,'ycolor','k','fontsize',fs-1)
ylim(ax,[0 1])
ylabel(ax,'cumulative PDF','fontsize',fs)
yyaxis(ax,'left')
set(ax,'ycolor','k')

% lines
for j = [3 1]
    h_psd(j) = plot(ax,0,0,'-','color',vel_clrs(j,:),'linewidth',1.5*lw);
    plot(ax,adcp(dep).fbin,adcp(dep).Pbin(:,j),'-','color',vel_clrs(j,:),'linewidth',lw)
    set(ax,'xscale','log','yscale','log')
end
% h_psd(1) = plot(ax,0,0,'-','color',vel_clrs(1,:),'linewidth',1.5*lw);
% plot(ax,adcp(dep).fbin,P_dep,'-','color',vel_clrs(1,:),'linewidth',lw)
% set(ax,'xscale','log','yscale','log')
    
lgd = legend(ax,h_psd([1 3 2]),{'u','w','cpdf'},'fontsize',fs-1);
% lgd = legend(ax,h_psd,{'P(u)+P(w)','cpdf'},'fontsize',fs-1);
lgd.Position(2) = 0.6;

yyaxis(ax,'right')
set(ax,'ytick',0:0.2:1)
yyaxis(ax,'left')

ylabel(ax,'\Psi_{u} [m^2/s^2 cps^{-1}]','fontsize',fs)
set(ax,'xtick',10.^(-6:1:2))
set(ax,'ytick',10.^(-6:2:2))
xlabel(ax,'f [cps]','fontsize',fs)
grid(ax,'off')
xlim(ax,[10^-4 2.5])
ylim(ax,10.^[-4.3 1.6])
set(ax,'FontSize',fs)

pnl_lbl = addPanelLabels(ax,fs);
pnl_lbl.String = 'i) velocity spectra';
pnl_lbl.FontWeight = 'bold';
