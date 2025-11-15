% Script to compile velocity observations from RHIBs near the glacier face,
% ultimately to try and bridge the 2024 meltstake velocity observations and
% the 2018 mooring observations.
%
% KJW
% 4 Nov 2025
clear

load('G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\adcp_combo_QC.mat')
load('G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\term_combo.mat')
load('G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\ctd_combo.mat')

% just keep glacier data
idx_glac = false(length(adcp),1);
for i = 1:length(adcp)
    idx_glac(i) = adcp(i).nuc_time(1) > datenum(2024,7,1);
end
adcp = adcp(idx_glac);
term = term(idx_glac);
ctd = ctd(idx_glac);
ndeps = sum(idx_glac);

%% plot tracks and terminus positions
x0 = 661157;
y0 = 6301824;

fig1 = figure(1); clf
ax1 = axes(fig1);
hold on; box on

% terminus
for i = 1:length(term)
    plot(ax1,term(i).X-x0,term(i).Y-y0,'k-','linewidth',1)
end

% RHIB tracks
max_dist = 200; % m
for i = 1:length(adcp)
    idxi = nearTerminus(term(i),adcp(i).vessel_X,adcp(i).vessel_Y,max_dist);
    adcp(i).idx_near = idxi;
    plot(ax1,adcp(i).vessel_X-x0,adcp(i).vessel_Y-y0,'-','color',0.5*[1 1 1])
    plot(ax1,adcp(i).vessel_X(idxi)-x0,adcp(i).vessel_Y(idxi)-y0,'-','color',colors(i),'linewidth',1)
end

% RHIB ctd casts
for i = 1:length(adcp)
    idxi = nearTerminus(term(i),ctd(i).X,ctd(i).Y,max_dist);
    ctd(i).idx_near = idxi;
    plot(ax1,ctd(i).X(idxi)-x0,ctd(i).Y(idxi)-y0,'o','color',colors(i),'linewidth',1)
    plot(ax1,ctd(i).X(~idxi)-x0,ctd(i).Y(~idxi)-y0,'o','color',0.5*[1 1 1])
end
axis equal

%% plot adcp deployments
% for i = 1:ndeps
%     % create plotting structure
%     adcp_plt = struct('burst',struct());
%     adcp_plt.burst.time = datetime(adcp(i).nuc_time,'convertfrom','datenum')';
%     adcp_plt.burst.range = adcp(i).cell_depth';
%     adcp_plt.burst.vel = permute(adcp(i).vel_QC,[3 1 2]);
%     adcp_plt.burst.samplerate = 1/seconds(median(diff(adcp_plt.burst.time)));
%     adcp_plt.burst.cellsize = mean(diff(adcp(i).cell_depth));
% 
%     % figure
%     figi = figure(i+1); clf
%     axi = adcpQuickPlot(figi,adcp_plt,'vel',0.4,NaT,[0 80],1);
%     title(axi(1),sprintf('%s - %s',adcp_plt.burst.time(1),adcp_plt.burst.time(end)))
% end

%% calculate profiles
KE_profs = nan(size(adcp(1).vel_QC,1),ndeps);
mom_profs = nan(size(adcp(1).vel_QC,1),ndeps);
vKE_profs = nan(size(adcp(1).vel_QC,1),ndeps);
vmom_profs = nan(size(adcp(1).vel_QC,1),ndeps);
nmeas_profs = nan(size(adcp(1).vel_QC,1),ndeps);

ctd_depth = ctd(1).depth;
T_profs = nan(size(ctd_depth,1),ndeps,3);
S_profs = nan(size(ctd_depth,1),ndeps,3);

% loop through deps
for i = 1:ndeps
    % time-dependent energy and momentum
    KEi = squeeze(0.5*sum(adcp(i).vel_QC(:,:,adcp(i).idx_near).^2,2));
    momi = squeeze(vecnorm(adcp(i).vel_QC(:,:,adcp(i).idx_near),2,2));
    % sample density
    nmeas_profs(:,i) = sum(~isnan(KEi),2);
    idxi = nmeas_profs(:,i)/max(nmeas_profs(:,i))>0.25;
    % mean profiles
    KE_profs(:,i) = mean(KEi,2,'omitnan');
    mom_profs(:,i) = mean(momi,2,'omitnan');
    vKE_profs(:,i) = sqrt(2*KE_profs(:,i));
    % nan out undersampled data
    KE_profs(~idxi,i) = nan;
    mom_profs(~idxi,i) = nan;
    vKE_profs(~idxi,i) = nan;
    
    % ctd data
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
    T_profs(1:iz_min,i,1) = min(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    T_profs(1:iz_min,i,2) = mean(ctd(i).T(1:iz_min,:),2,'omitnan');
    T_profs(1:iz_min,i,3) = max(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    S_profs(1:iz_min,i,1) = min(ctd(i).S(1:iz_min,:),[],2,'omitnan');
    S_profs(1:iz_min,i,2) = mean(ctd(i).S(1:iz_min,:),2,'omitnan');
    S_profs(1:iz_min,i,3) = max(ctd(i).S(1:iz_min,:),[],2,'omitnan');
end

% mean across all deps
vKE_mean = mean(vKE_profs,2,'omitnan');
mom_mean = mean(mom_profs,2,'omitnan');
mom_std = std(mom_profs,0,2,'omitnan');
T_mean = squeeze(mean(T_profs,2,'omitnan'));
S_mean = squeeze(mean(S_profs,2,'omitnan'));

% plot
fig2 = figure(200); clf
clear ax2
for i = 3:-1:1
    ax2(i) = axes(fig2,'position',axgridpos(1,3,i,[.05 .05 .1 .1],[0 0]));
    hold on; box on; axis ij
end

for i = 1:ndeps
    plot(ax2(1),vKE_profs(:,i),adcp(i).cell_depth,'color',colors(i))
    plot(ax2(2),mom_profs(:,i),adcp(i).cell_depth,'color',colors(i))
    plot(ax2(3),nmeas_profs(:,i)/max(nmeas_profs(:,i)),adcp(i).cell_depth,'color',colors(i))
end
plot(ax2(1),vKE_mean,adcp(i).cell_depth,'k-','linewidth',1)
plot(ax2(2),mom_mean,adcp(i).cell_depth,'k-','linewidth',1)

linkaxes(ax2,'y')

fig3 = figure(201); clf
hold on; box on
plot(vKE_mean,adcp(i).cell_depth,'r-','linewidth',1)
plot(mom_mean,adcp(i).cell_depth,'k-','linewidth',1)
axis ij
xlabel('speed [m/s]')
ylabel('depth [m]')
xlim([0 0.4])

fig4 = figure(202); clf
axT = subplot(1,2,1); axS = subplot(1,2,2);
hold(axT,'on'); hold(axS,'on')
axis(axT,'ij'); axis(axS,'ij');
for i = 1:3
    plot(axT,T_mean(:,i),ctd_depth)
    plot(axS,S_mean(:,i),ctd_depth)
end

% melt rate
T_3eqn = interp1(ctd_depth,T_mean(:,2),adcp(1).cell_depth);
S_3eqn = interp1(ctd_depth,S_mean(:,2),adcp(1).cell_depth);
[m_3eqn,Tb,Sb] = solve3Eqn(mom_mean,T_3eqn,S_3eqn,adcp(1).cell_depth);

%% save locally
adcp_depth = adcp(1).cell_depth;
save platform_comparison_data\polly_data.mat mom_mean mom_std T_mean S_mean adcp_depth ctd_depth max_dist m_3eqn

%% functions
% % % % %
function idx = nearTerminus(term,vessel_X,vessel_Y,max_dist)
    nt = length(vessel_X);
    dist2 = nan(nt,1);
    for i = 1:nt
        Xi = vessel_X(i);
        Yi = vessel_Y(i);
        dist2(i) = min((term.X-Xi).^2 + (term.Y-Yi).^2);
    end
    idx = dist2 <= max_dist^2;
end