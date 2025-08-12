% Script to make a summary figure from the glacier meltstake deployments
% during the July cruise in 2024
%
% KJW
% 8 Aug 2024

clear

data_path = 'F:AK_202406/data/meltstake';

% deployments
deps = {'ms02_20240712_2020','ms02_20240715_2001','ms01_20240717_0050'};

% time limits to plot
time_lims = {[datetime(2024,7,12,21,10,0) datetime(2024,7,12,21,40,0)],...
             [datetime(2024,7,15,22,21,0) datetime(2024,7,15,22,58,0)],...
             [datetime(2024,7,17,1,12,0) datetime(2024,7,17,1,55,0)]};

% adv beams to plot ice distance
adv_bn = [1, 3, 2];

% load data
adcp_all = {};
adv_all = {};
T_all = {};

ndeps = length(deps);
for i = 1:ndeps
    load(fullfile(data_path,deps{i},'adcp','adcp.mat'))
    load(fullfile(data_path,deps{i},'adv','adv.mat'))
    load(fullfile(data_path,deps{i},'rbr','T.mat'))

    adcp_all{i} = adcp;
    adv_all{i} = adv;
    T_all{i} = T;
end

adcp = [adcp_all{1} adcp_all{2} adcp_all{3}];
adv = [adv_all{1} adv_all{2} adv_all{3}];
T = [T_all{1}; T_all{2}; T_all{3}];
clear adcp_all adv_all T_all

%% plot
lw = 1;
fs = 11;
set(0,'defaulttextinterpreter','latex')

pad = [.05 .03 .1 .07];

panel_lbls = {'u','w','T (Solos)','ADV echosounder'};

r_max = 0.85;

% melt rate guesses
d0 = [150 130 125]/1000;
mr = [1.5 1.8 2.5];

% T colors
T_clr = [0.1 0.8 0.7]/1.;
T_clrs = T_clr'*[1 0.6 0.2];
T_lgd_lbls = {'near','mid','far'};

% create and size figure
figsize = [16 16];
fig = figure(1); clf
set(fig,'units','centimeters'); 
fpos = get(fig,'position');
fpos(3:4) = figsize;
set(fig,'position',fpos);
set(fig,'paperunits','centimeters')
set(fig,'papersize',figsize)

clear ax
% loop through deps
for i = 1:ndeps
    % time limits
    t1 = time_lims{i}(1)-seconds(30);
%     t2 = time_lims{i}(2);
    t2 = t1 + minutes(30);
    idx_adcp_att = adcp(i).attitude.time>=t1 & adcp(i).attitude.time<=t2;
    idx_adcp_brst = adcp(i).burst.time>=t1 & adcp(i).burst.time<=t2;
    idx_adv = adv(i).pck_time>=t1 & adv(i).pck_time<=t2;
    idx_T = T(i,1).time>=t1 & T(i,1).time<=t2;

    % adcp range limit
    idx_adcp_range = adcp(i).burst.range<=r_max;

    % deployment depth
    d = mean(adcp(i).attitude.pressure(idx_adcp_att));

    % adcp u and w
    for j = 1:2
        ax(i,j) = axes(figure(1),'position',axgridpos(4,ndeps,4*(i-1)+j,pad,'flip'));
        pcolor(adcp(i).burst.time(idx_adcp_brst),adcp(i).burst.range(idx_adcp_range),adcp(i).burst.vel_ice(idx_adcp_brst,idx_adcp_range,j)')
        shading flat
        cmocean('bal')
        clim(0.3*[-1 1])
        
    end

    % T
    ax(i,3) = axes(figure(1),'position',axgridpos(4,ndeps,4*(i-1)+3,pad,'flip'));
    hold on
    box on
    for j = 1:3
        plot(T(i,j).time(idx_T),T(i,j).values(idx_T),'linewidth',lw,'color',T_clrs(:,j))
        ylim([4 7.5])
    end

    % adv distance
    ax(i,4) = axes(figure(1),'position',axgridpos(4,ndeps,4*(i-1)+4,pad,'flip'));
    hold on
    pcolor(adv(i).pck_time(idx_adv),adv(i).pck_dist(1,:)/1000,adv(i).pck_amp(idx_adv,:,adv_bn(i))')
    shading flat
    clim([110 195])
    plot([t1 t2],[d0(i) d0(i)+days(t2-t1)*mr(i)],'k--','linewidth',lw)
    ylim([100 250]/1000)
    text(ax(i,4),0.5,0,sprintf('%.1f m/day',mr(i)),'units','normalized','horizontalalignment','center',...
        'verticalalignment','bottom','backgroundcolor','w')

    % linkaxes
    linkaxes(ax(i,:),'x')
    xlim(ax(i,1),[t1+seconds(30) t2])

    % check pitch bias
    figure(i+1); clf
    plot(adcp(i).attitude.time(idx_adcp_att),adcp(i).attitude.roll(idx_adcp_att),'k')

    % titles
    title(ax(i,1),{sprintf('Deployment %d',i),sprintf('%d m',round(d))},'fontsize',fs)

    % panels labels
    for j = 1:4
        text(ax(i,j),.03,0.99,panel_lbls{j},'units','normalized','verticalalignment','top',...
            'fontsize',fs,'FontWeight','bold')
    end

end

% axes labels
ylabel(ax(1,1),{'range','[m]'},'fontsize',fs)
ylabel(ax(1,2),{'range','[m]'},'fontsize',fs)
ylabel(ax(1,3),{'temperature','[$^\circ$C]'},'fontsize',fs)
ylabel(ax(1,4),{'distance','[m]'},'fontsize',fs)

% velocity colorbar
cbar = colorbar(ax(3,2),'position',cbarpos(ax(3,2),.005,.02));
cbar.Label.String = 'm/s';
cbar.Label.FontSize = fs-2;
cbar.Position(4) = 2*cbar.Position(4)+pyi;

% T legend
legend(ax(3,3),T_lgd_lbls,'location','southeast','fontsize',fs-3)