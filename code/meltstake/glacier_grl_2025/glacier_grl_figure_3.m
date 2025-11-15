% Script to plot boundary layer details for each deployment. This should
% answer zero-th order questions about what the boundary layer forcing
% looks like.
%
% KJW
% 12 Nov 2025
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

load platform_comparison_data\ms_psd.mat
load glacier_clrs.mat

dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

%% load hobo, T, polly ctd data
clear hobo T ctd
for i = 1:3
    % hobo
    hobo_load = load(fullfile(raw_dir,dep_tbl.Folder{i},'hobo','hobo.mat'));
    hobo(i) = hobo_load.hobo(2);
    % solos
    T_load = load(fullfile(raw_dir,dep_tbl.Folder{i},'rbr','T.mat'));
    T(i) = struct('time',T_load.T(2).time,'T',T_load.T(2).values);
    % polly ctd
    ctd_load = load(fullfile(proc_dir,dep_tbl.Folder{i},'ctd.mat'));
    ctd(i) = ctd_load.ctd(strcmp({ctd_load.ctd.vessel},'polly'));
end

%% average outer TS conditions
dz = 2;
ms_depth = nan(3,1);
for i = 1:3
    ms_depth(i) = mean(seg_tbl.depth(seg_tbl.Number==dep_tbl.Number(i)));
    idxd = abs(ctd(i).depth-ms_depth(i)) <= dz;
    ctd(i).T_ms = max(ctd(i).T(idxd,:),[],'all','omitnan');
    ctd(i).S_ms = max(ctd(i).S(idxd,:),[],'all','omitnan');
end

%% prep adcp data
% where to sample for PSDs and velocities
y_s = [0.12 0.21] + 1e-6*[-1 1]; % distance from wall 
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

% bin PSDs
% binning setup
nf = [4 8 16 32];
nf_bins = [0 1e-2 1e-1 1e0];

for i = 1:length(adcp)
    [Pb,fb,M,S_ci] = progressiveBin(adcp(i).f,adcp(i).Puvw(:,1),nf,nf_bins);
    adcp(i).fbin = fb;
    adcp(i).Pbin = repmat(Pb,[1 3]);
    adcp(i).Mbin = M;
    adcp(i).S_ci = S_ci;
    for j = 2:3
        adcp(i).Pbin(:,j) = progressiveBin(adcp(i).f,adcp(i).Puvw(:,j),nf,nf_bins);
    end
end

% CPDFs of BL energy
for i = 1:length(adcp)
    fbi = adcp(i).fbin;
    nfb = length(fbi);
    
    % straight CPDF
    adcp(i).CPDF = cumtrapz(fbi,sum(adcp(i).Pbin,2));
    
    % CPDF assuming remaining variance follows f^-5/3
    CPDF_turb = nan(nfb,1);
    
    %fig = figure(i+1); axi = axes(fig); hold on; set(axi,'xscale','log','yscale','log')

    Pb_all = sum(adcp(i).Pbin,2);
    for j = 1:nfb
        Pj_0 = Pb_all(j);
        fj_0 = fbi(j);
        Aj = Pj_0*fj_0^(5/3);
        Pj = [Pb_all(1:j-1); Aj*fbi(j:end).^(-5/3)];
        CPDF_turb(j) = trapz(fbi,Pj);
        if ~mod(j,5)
            %plot(axi,fbi,Pj)
        end
    end
    adcp(i).CPDF_turb = CPDF_turb;
end

%% figure
% set up plot
lw = 0.8;
fs = 10;
alpha = 0.5;

fig_size = [15 12.5]*1.;
fig = figure(1); clf
setFigureSize(fig,fig_size);
clear ax

pad = [.02 .06 .16 .08];
shift = [-0.05 0.01];

f_slope = [1e-1 1e0];
A0 = 8e-4;
A = A0*f_slope(1)^(5/3);

% create axes
for i = 1:9
    ax(i) = axes(fig,'position',axgridpos(3,3,i,pad,shift),'fontsize',fs-1);
    box(ax(i),'on')
    hold(ax(i),'on')
end
% shift middle row up
for i = 4:6
    ax(i).Position(2) = ax(i).Position(2) + .04;
end
% shrink bottom row in the x dir a little to make room for labels
for i = 7:9
    ax(i).Position(3) = ax(i).Position(3) - .015;
end

% velocity PSDs
clear h_psd
for i = 1:length(adcp)
    % CPDFs
    yyaxis(ax(i+3),'right')
    h_psd(4) = plot(ax(i+3),adcp(i).fbin,adcp(i).CPDF/adcp(i).CPDF(end),'--','color',0.5*[1 1 1]);
    h_psd(5) = plot(ax(i+3),adcp(i).fbin,hannFilter(adcp(i).CPDF_turb,10)/adcp(i).CPDF(end),'-','color',0.5*[1 1 1]);
    set(ax(i+3),'ytick',[],'ycolor','k')
    ylim(ax(i+3),[0 1])
    yyaxis(ax(i+3),'left')
    set(ax(i+3),'ycolor','k')

    % shaded regions
    for j = 1:3
        plotConfidenceRegion(ax(i+3),adcp(i).fbin,adcp(i).S_ci(:,1).*adcp(i).Pbin(:,j),adcp(i).S_ci(:,2).*adcp(i).Pbin(:,j),vel_clrs(j,:),alpha);
    end
    % lines
    for j = 1:3
        h_psd(j) = plot(ax(i+3),0,0,'-','color',vel_clrs(j,:),'linewidth',1.5*lw);
        plot(ax(i+3),adcp(i).fbin,adcp(i).Pbin(:,j),'-','color',vel_clrs(j,:),'linewidth',lw)
        set(ax(i+3),'xscale','log','yscale','log')
    end
    % slope line
    plot(ax(i+3),f_slope,A*f_slope.^(-5/3),'k-')
    
end
lgd1 = legend(ax(6),h_psd,{'u','v','w','CPDF','+f^{-5/3} est.'},'fontsize',fs-1);
lgd1.Position(1) = sum(ax(6).Position([1 3]))+.01;

% wind roses
vbins = [0 .05 .1 .25 .4];
dir_lbls = {{'+w', '', '-w', '+u'},{'+w', '', '-w', ''},{'+w', '-u', '-w', ''}};
options = {'axes',[],'labels',{'+w', '-u', '-w', '+u'},'freqlabelangle','none',...
    'vwinds',[0 .05 .1 .2 .3],'titlestring','','lablegend','','cmap',cmocean('-mat'),...
    'legendtype',0,'ndirections',24,'gridcolor',0.45*[1 1 1]};
for i = 1:3
    % plot data
    ui = -adcp(i).uvw(:,1);
    wi = adcp(i).uvw(:,3);
    spd = sqrt(ui.^2 + wi.^2);
    dir = atan2(wi,ui)*180/pi;
    % plotting
    options{2} = ax(i);
    options{4} = {'','','',''};%dir_lbls{i};
    WindRose(dir,spd,options);
    text(ax(i),.02,0.5,'+u','fontsize',fs-2,'units','normalized')
    text(ax(i),1-.1,0.5,'-u','fontsize',fs-2,'units','normalized')
    text(ax(i),0.5,0.95,'+w','fontsize',fs-2,'units','normalized','horizontalalignment','center')
    text(ax(i),0.5,0.06,'-w','fontsize',fs-2,'units','normalized','horizontalalignment','center')
end
% legend
idx_patch = zeros(length(vbins),1);
pos = 0;
for i = 1:length(ax(3).Children)
    if isa(ax(3).Children(i),'matlab.graphics.primitive.Patch')
        pos = pos + 1;
        idx_patch(pos) = i;
    end
    if prod(idx_patch)
        break
    end
end
h_wr = ax(3).Children(idx_patch);
wr_lbls = cell(length(vbins),1);
for i = 1:length(vbins)-1
    wr_lbls{i} = sprintf('[%.2f,%.2f)',round(vbins(i),2),round(vbins(i+1),2));
end
wr_lbls{end} = sprintf('>=%.2f',vbins(end));
lgd2 = legend(ax(3),h_wr,wr_lbls,'fontsize',fs-2);
lgd2.Position(1) = sum(ax(3).Position([1 3]))+.01;
lgd2.Title.String = '(u^2+w^2)^{1/2}';
lgd2.Title.FontSize = fs-2;


% T-S plots
S0 = [26.7 27.8 28];
T0 = [6 6.9 7.1];
fmix = [12/100 2/100];

ms = 4;
for i = 1:3
    % mixing lines
    T_ml = [T0(i)*(1-fmix(1)) T0(i) T0(i)*(1-fmix(2))-90*fmix(2)];
    S_ml = [S0(i)*(1-fmix(1)) S0(i) S0(i)*(1-fmix(2))];
    plot(ax(6+i),S_ml,T_ml,'k-','linewidth',lw)
    % TS points
    idxt = hobo(i).time >= dep_tbl.Start(i) & hobo(i).time <= dep_tbl.End(i);
    ti = hobo(i).time(idxt); ti = hours(ti-ti(1));
    Si = hannFilter(hobo(i).S(idxt),1);
    Ti = hannFilter(hobo(i).T_cal(idxt),1);
    plot(ax(6+i),Si,Ti,'k-')
    scatter(ax(6+i),Si,Ti,ms,ti,'filled')
    cmocean('-rain',ax(6+i))
end
cbar = colorbar(ax(9),'position',cbarpos(ax(9),.01,.02));
cbar.Label.String = 'time [hr]';
cbar.Label.FontSize = fs-1;

% dep labels
for i = 1:3
    text(ax(i),0.5,1.02,sprintf('%d (%d m, %d min)',i,round(ms_depth(i)),round(dep_tbl.Duration(i)*24*60)),'fontsize',fs,'fontweight','bold','horizontalalignment','center','verticalalignment','bottom','units','normalized')
end

% axes labels
% clear
linkaxes(ax(4:6))
% linkaxes(ax(7:9))
for i = [5 6]
    set(ax(i),'yticklabel',{})
end
% psds
ylabel(ax(4),'\Psi_{u} [m^2/s^2 cps^{-1}]','fontsize',fs)
for i = 4:6
    set(ax(i),'xtick',10.^(-6:2:2))
    set(ax(i),'ytick',10.^(-6:2:2))
    xlabel(ax(i),'f [cps]','fontsize',fs)
    grid(ax(i),'off')
end
xlim(ax(4),[10^-4 6])

% TS
for i = 7:9
    xlabel(ax(i),'S [psu]','fontsize',fs)
    clim(ax(i),[0 3])
end
ylabel(ax(7),'T [\circC]','fontsize',fs)





