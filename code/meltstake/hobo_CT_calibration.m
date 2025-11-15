% Script to calibrate Hobo conductivity and temperature data. The Hobo time
% axis should be calibrated first using "hobo_time_calibration.m"
%
% 10 Oct 2025

clear

% load
raw_dir = 'F:/meltstake/data/raw';
load glacier_grl_2025\glacier_clrs.mat

% choose deployment
ms_tbl = loadMSInfo(27);
save_data = 1;

% specify C calibration times (when ROV CTD and Hobos should be ~same)
switch ms_tbl.Number
    case 26
        C_cal_time = [datetime(2024,7,12,20,17,5.5) datetime(2024,7,12,20,19,45.5)];
    case 27
        C_cal_time = [datetime(2024,7,15,19,58,17) datetime(2024,7,15,20,0,30)];
    case 28
        C_cal_time = [datetime(2024,7,17,0,47,42) datetime(2024,7,17,0,51,38)];
end

% hobos, solos, and rov ctd
load(fullfile(raw_dir,ms_tbl.Folder{1},'ctd','ctd.mat'))
load(fullfile(raw_dir,ms_tbl.Folder{1},'hobo','hobo.mat'))
load(fullfile(raw_dir,ms_tbl.Folder{1},'rbr','T.mat'))


%% plot stuff
lw = 1;
fs = 11;

% ctd
fig_ctd = figure(1); clf
clear ax_ctd
fld_nums = [4 2 1 5];
idxt_ctd = ctd.values(:,4) > 1;
TLIM = extrema(ctd.time(idxt_ctd));
% TLIM(2) = datetime(2024,7,15,20,1,10);
for i = 1:length(fld_nums)
    ax_ctd(i) = axes(fig_ctd,'position',axgridpos(length(fld_nums),1,i,[0 .05 .1 .08],[.04 0]));
    plot(ax_ctd(i),ctd.time(idxt_ctd),ctd.values(idxt_ctd,fld_nums(i)),'linewidth',lw)
    ylabel(ax_ctd(i),{ctd.channels{fld_nums(i)},sprintf('[%s]',ctd.units{fld_nums(i)})},'fontsize',fs)
end
linkaxes(ax_ctd,'x')
title(ax_ctd(1),'CTD','fontsize',fs)
xlim(ax_ctd(1),TLIM)

% hobos
fig_hobo = figure(2); clf
clear ax_hobo
for i = 1:2
    ax_hobo(i) = axes(fig_hobo,'position',axgridpos(2,1,i,[0 .05 .1 .08],[0 0]));
    hold on
    box on
end
for i = 1:length(hobo)
    idxt_hobo = hobo(i).time>=(TLIM(1)-minutes(0)) & hobo(i).time<=(TLIM(2)+minutes(0));
    plot(ax_hobo(1),hobo(i).time(idxt_hobo),hobo(i).C(idxt_hobo)/1000,'color',colors(i),'linewidth',lw)
    plot(ax_hobo(2),hobo(i).time(idxt_hobo),hobo(i).T(idxt_hobo),'color',colors(i),'linewidth',lw)
end
linkaxes(ax_hobo,'x')
ylabel(ax_hobo(1),'conducitivity [mS/cm]','fontsize',fs)
ylabel(ax_hobo(2),'temperature [^\circC]','fontsize',fs)
title(ax_hobo(1),'Hobos and CTD (black)','fontsize',fs)
% add ctd measurements
plot(ax_hobo(1),ctd.time(idxt_ctd),ctd.values(idxt_ctd,1),'k-','linewidth',lw/2)
plot(ax_hobo(2),ctd.time(idxt_ctd),ctd.values(idxt_ctd,2),'k-','linewidth',lw/2)

% solos
fig_sol = figure(3); clf
ax_sol = axes(fig_sol);
hold on; box on
for i = 1:length(T)
    idxt_sol = T(i).time>=TLIM(1) & T(i).time<=TLIM(2);
    plot(ax_sol,T(i).time(idxt_sol),T(i).values(idxt_sol),'color',colors(i),'linewidth',lw)
end
title(ax_sol,'solos and CTD (black)','fontsize',fs)
ylabel(ax_sol,'temperature [^\circC]','fontsize',fs)
% add ctd measurements
plot(ax_sol,ctd.time(idxt_ctd),ctd.values(idxt_ctd,2),'k-','linewidth',lw/2)

% compare solo and hobo T
fig_cmp1 = figure(4); clf
clear ax_cmp1
idx_sol = T(1).time>=ms_tbl.Start & T(1).time<=ms_tbl.End;
for i = 1:3
    ax_cmp1(i) = axes(fig_cmp1,'position',axgridpos(3,1,i,[0 .05 .1 .08],[0 0]));
    hold on; box on
    for j = i:(i+1)
        if j>3
            continue
        end
        % plot solos around hobo
        plot(ax_cmp1(i),T(j).time(idx_sol),T(j).values(idx_sol),'color',colors(j),'linewidth',lw)
    end
    % plot hobo
    idx_hob = hobo(i).time>=ms_tbl.Start & hobo(i).time<=ms_tbl.End;
    plot(ax_cmp1(i),hobo(i).time(idx_hob),hobo(i).T(idx_hob),'k-','linewidth',lw)
    ylabel(ax_cmp1(i),sprintf('temp (hobo %d)',i),'fontsize',fs)
end
linkaxes(ax_cmp1,'x')
title(ax_cmp1(1),'solo-hobo T comparison','fontsize',fs)

% profile of mean T
y = [25 30 39 46 54 61]-16;
fig_cmp2 = figure(5); clf
ax_cmp2 = axes(fig_cmp2);
hold on; box on
sol_mean = [nan nan nan];
sol_std = [nan nan nan];
hob_mean = [nan nan nan];
hob_std = [nan nan nan];
for i = 1:3
    % solo
    sol_mean(i) = mean(T(i).values(idx_sol));
    sol_std(i) = std(T(i).values(idx_sol));
    % hobo
    idx_hob = hobo(i).time>=ms_tbl.Start & hobo(i).time<=ms_tbl.End;
    hob_mean(i) = mean(hobo(i).T(idx_hob));
    hob_std(i) = std(hobo(i).T(idx_hob));
    errorbar(ax_cmp2,y(2*(i-1)+1),sol_mean(i),sol_std(i),'ko','markerfacecolor',colors(i),'linewidth',lw,'markersize',8)
    errorbar(ax_cmp2,y(2*(i-1)+2),hob_mean(i),hob_std(i),'ksquare','markerfacecolor',colors(i),'linewidth',lw,'markersize',8)
end
xlabel(ax_cmp2,'y [cm]','fontsize',fs)
ylabel(ax_cmp2,'temperature [^\circC]','fontsize',fs)
title(ax_cmp2,'mean temperature profile','fontsize',fs)


%% scatterplot of temperature and calibration
fig_cmp3 = figure(6); clf
clear ax_cmp3
ms = 4;
for i = 1:3
    ax_cmp3(i) = axes(fig_cmp3,'position',axgridpos(1,3,i,[.05 0 .08 .12],[0 0]));
    hold on; box on
    % interpolate solo T from adjacent sensors
    if i<3
        T_solo = (T(i).values(idx_sol)+T(i+1).values(idx_sol))/2;
    else
        T_solo = T(i).values(idx_sol);
    end
    idx_hob = hobo(i).time>=ms_tbl.Start & hobo(i).time<=ms_tbl.End;
%     T_slow = hannFilter(T_solo,2*10);
    T_slow = exponentLag(T_solo,0.5,2,0.1);
    T_interp = interp1(T(i).time(idx_sol),T_solo,hobo(i).time(idx_hob));
    T_interp_slow = interp1(T(i).time(idx_sol),T_slow,hobo(i).time(idx_hob));

    % calibrate
    hobo(i).calibration = struct();
    [B,BINT] = regress(T_interp_slow,[ones(size(T_interp_slow)) hobo(i).T(idx_hob)]);
    hobo(i).T_cal = B(2)*hobo(i).T+B(1);
    hobo(i).calibration.T_method = 'linear';
    hobo(i).calibration.T_coeffs = B;
    hobo(i).calibration.T_coeff_ints = BINT;
    plot(ax_cmp1(i),hobo(i).time(idx_hob),hobo(i).T_cal(idx_hob),'k:','linewidth',lw)
    plot(ax_cmp2,y(2*i),mean(hobo(i).T_cal(idx_hob)),'k^','markerfacecolor',colors(i),'linewidth',lw,'markersize',8)

    % plot
    plot(ax_cmp3(i),[0 10],[0 10],'k-')
    plot(ax_cmp3(i),T_interp,hobo(i).T(idx_hob),'k.','markersize',ms)
    plot(ax_cmp3(i),T_interp_slow,hobo(i).T(idx_hob),'r.','markersize',ms)
    plot(ax_cmp3(i),T_interp_slow,hobo(i).T_cal(idx_hob),'.','color',colors(5),'markersize',ms)

    xlim(ax_cmp3(i),extrema(T_interp)+0.5*[-1 1]); ylim(ax_cmp3(i),extrema(T_interp)+0.5*[-1 1])
    xlabel(ax_cmp3(i),'solo T [^\circC]','fontsize',fs)
    ylabel(ax_cmp3(i),sprintf('hobo %d T [^\\circC]',i),'fontsize',fs)

end

%% scatterplot of conductivity and calibration
fig_cmp4 = figure(7); clf
clear ax_cmp4
ms = 6;
for i = 1:3
    ax_cmp4(i) = axes(fig_cmp4,'position',axgridpos(1,3,i,[.05 0 .08 .12],[0 0]));
    hold on; box on

    % trim data
    idx_hob = hobo(i).time>=C_cal_time(1) & hobo(i).time<=C_cal_time(2);
    idx_ctd = ctd.time>=C_cal_time(1) & ctd.time<=C_cal_time(2);

    % slow down CTD
    C_slow = interp1(ctd.time(idx_ctd),hannFilter(ctd.values(idx_ctd,1),16),hobo(i).time(idx_hob));
    
    % calibrate
    [B,BINT] = regress(C_slow,hobo(i).C(idx_hob)/1e3);
    hobo(i).C_cal = B*hobo(i).C/1e3;
    hobo(i).calibration.C_method = 'scale';
    hobo(i).calibration.C_coeffs = B;
    hobo(i).calibration.C_coeffs_ints = BINT;
    hobo(i).calibration.C_cal_time = C_cal_time;

    % plot
    plot(ax_hobo(1),hobo(i).time(idx_hob),hobo(i).C_cal(idx_hob),'--','color',colors(i))
    % scatter
    plot(ax_cmp4(i),[20 30],[20 30],'k-')
    plot(ax_cmp4(i),interp1(ctd.time,ctd.values(:,1),hobo(i).time(idx_hob)),hobo(i).C(idx_hob)/1e3,'k.','markersize',ms)
    plot(ax_cmp4(i),C_slow,hobo(i).C(idx_hob)/1e3,'r.','markersize',ms)
    plot(ax_cmp4(i),C_slow,hobo(i).C_cal(idx_hob),'.','color',colors(5),'markersize',ms)

    xlim(ax_cmp4(i),extrema(C_slow)+0.5*[-1 1]); ylim(ax_cmp4(i),extrema(C_slow)+0.5*[-1 1])
    xlabel(ax_cmp4(i),'ctd C [mS/cm]','fontsize',fs)
    ylabel(ax_cmp4(i),sprintf('hobo %d C [mS/cm]',i),'fontsize',fs)

end

%% compute and compare available salinity
C_lag = 2; % s, slow down conductivity to match temperature response time
for i = 1:length(hobo)
    hobo(i).calibration.C_lag = C_lag;
    hobo(i).S = gsw_SP_from_C(exponentLag(hobo(i).C_cal,1,2,0.1),hobo(i).T_cal,0);
end

fig_cmp5 = figure(8); clf
clear ax_cmp5
ms = 6;
for i = 1:3
    ax_cmp5(i) = axes(fig_cmp5,'position',axgridpos(1,3,i,[.05 0 .08 .12],[0 0]));
    hold on; box on

    % trim data
    idx_hob = hobo(i).time>=C_cal_time(1) & hobo(i).time<=C_cal_time(2);
    idx_ctd = ctd.time>=C_cal_time(1) & ctd.time<=C_cal_time(2);

    % slow down CTD
    S_slow = interp1(ctd.time(idx_ctd),hannFilter(ctd.values(idx_ctd,5),16),hobo(i).time(idx_hob));

    % plot
    plot(ax_cmp5(i),[20 30],[20 30],'k-')
    plot(ax_cmp5(i),interp1(ctd.time,ctd.values(:,5),hobo(i).time(idx_hob)),hobo(i).S(idx_hob),'k.','markersize',ms)
    plot(ax_cmp5(i),S_slow,hobo(i).S(idx_hob),'r.','markersize',ms)

    xlim(ax_cmp5(i),extrema(S_slow)+0.5*[-1 1]); ylim(ax_cmp5(i),extrema(S_slow)+0.5*[-1 1])
    xlabel(ax_cmp5(i),'ctd S [PSU]','fontsize',fs)
    ylabel(ax_cmp5(i),sprintf('hobo %d S [PSU]',i),'fontsize',fs)
end

%% get rid of bad data due to bubbles/sed
fig_QC = figure(9); clf
ax_QC = axes(fig_QC);
hold on; box on
for i = 1:length(hobo)
    plot(hobo(i).time,hobo(i).C)
end

%% T-S plots
fig_TS = figure(10); clf
clear ax_TS
ms = 8;
max_S = 0;
max_T = 0;
for i = 1:length(hobo)
    ax_TS(i) = axes(fig_TS,'position',axgridpos(1,3,i,[.05 0 .08 .12],[-.02 0]));
    hold on; box on
    % ctd
    plot(ctd.values(idx_ctd,5),ctd.values(idx_ctd,2),'k.','markersize',ms)
    % hobos
    for j = setxor(1:length(hobo),i)
        idx_hob = hobo(j).time>=ms_tbl.Start & hobo(j).time<=ms_tbl.End;
        plot(hobo(j).S(idx_hob),hobo(j).T_cal(idx_hob),'.','color',0.5*[1 1 1],'markersize',ms)
    end
    idx_hob = hobo(i).time>=C_cal_time(1) & hobo(i).time<=ms_tbl.End;
    scatter(hobo(i).S(idx_hob),hobo(i).T_cal(idx_hob),ms,minutes(hobo(i).time(idx_hob)-min(hobo(i).time(idx_hob))),'o','filled')
    max_S = max(max_S,max(hobo(i).S(idx_hob)));
    max_T = max(max_T,max(hobo(i).T_cal(idx_hob)));
    xlabel(ax_TS(i),'S [psu]','fontsize',fs)
    ylabel(ax_TS(i),'T [^\circC]','fontsize',fs)
    title(ax_TS(i),sprintf('hobo %d',i),'fontsize',fs)

end

% mixing lines
f_sgd = 0.16;
f_mlt = .015;
max_S = max_S + 0.1;
max_T = max_T + 0.05;
for i = 1:length(hobo)
    % discharge line
    plot(ax_TS(i),max_S*[1-f_sgd 1],max_T*[1-f_sgd 1],'k-')
    % melt line
    plot(ax_TS(i),max_S*[1-f_mlt 1],max_T*[(1-f_mlt)-90/max_T*f_mlt 1],'k-')
end
cbar = colorbar(ax_TS(end),'position',cbarpos(ax_TS(end),.01,.02),'fontsize',fs);
cbar.Label.String = 'time [min]';
cbar.Label.FontSize = fs;
linkaxes(ax_TS)

%% save updated hobo data
if save_data
%     save(fullfile(raw_dir,ms_tbl.Folder{1},'hobo','hobo.mat'),'hobo')
end