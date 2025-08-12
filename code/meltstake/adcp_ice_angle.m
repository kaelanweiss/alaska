% Script to try and calculate ice slope relative to the adcp using vertical
% beam ice intersections.
%
% KJW
% 28 Mar 2025

clear

raw_dir = 'F:/meltstake/data/raw';

ms_tbl = loadMSInfo(26:28,'manualwindows');
mcorr_tbl = loadMSInfo(26:28,'pitch_correction');

[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);

%% figure
% figsize = [22 10];
% pad = [.07 .05 .07 .12];
% shift = [0 0];

fig2 = figure(2); clf

pos = 0;
r_ice = nan(23,2); % beam 2, beam 4

% loop through deployments
for i = 1:ndeps
    seg_nums = ms_tbl.Window(ms_tbl.Number==dep_nums(i));

    % load adcp data
    load(fullfile(raw_dir,dep_names{i},'adcp/adcp.mat'))
    adcp.burst.vel = adcp.burst.vel(:,:,[2 4]);
    
    % loop through segments
    for j = 1:length(seg_nums)
        row_num = ms_tbl.Number == dep_nums(i) & ms_tbl.Window == j;
        ax2 = adcpQuickPlot(fig2,adcp,'vel',0.05*[-1 1],[ms_tbl.Start(row_num) ms_tbl.End(row_num)],[0 1],4);

        a = 0;

    end
end

%% analytical analysis of ice intersections and implied ice slope

% tan(beta) is given by the following expression (beta>0 is undercut) 
% tan(beta) = (1/tan(25))*(r4-r2)/(r2+r4+2r0)
r0 = 0.16;
Tbeta_func = @(r2,r4) (r4-r2)./(r2+r4+2*r0)/tand(25);

% uncertainty propagation delta(tan(beta)):
dr = 0.03;
dT_func = @(r2,r4) (2*dr/tand(25))*sqrt((r2 + r0).^2 + (r4 + r0).^2)./(r2 + r4 + 2*r0).^2;

% plotting grid
r = linspace(0.35,0.8,100);
[R2,R4] = meshgrid(r,r);

% calculate quanitities
% tan(beta)
Tbeta = Tbeta_func(R2,R4);
dT = dT_func(R2,R4);
Tbeta_ci = cat(3,Tbeta-dT,Tbeta+dT);
% beta
beta = atand(Tbeta);
beta_ci = atand(Tbeta_ci);
dbeta = diff(beta_ci,1,3);

% plot
fig1 = figure(1); clf
ax1 = axes(fig1);
hold on; box on
pcolor(r,r,dbeta./abs(beta))
shading flat
[Cb,Hb] = contour(r,r,beta,'k');
[Cdb,Hdb] = contour(r,r,dbeta,0:0.5:10,'r');
clabel(Cb,Hb)
clabel(Cdb,Hdb,'color','r')
xlabel('r_2')
ylabel('r_4')
legend([Hb,Hdb],{'\beta','\delta\beta'},'fontsize',11,'autoupdate',false)
colormap gray
clim([0 1])
cbar = colorbar;
cbar.Label.String = '|\delta\beta/\beta|';
cbar.Label.FontSize = 11;

if ismember('r2_ice',mcorr_tbl.Properties.VariableNames)
    for i = 1:ndeps
        dep_idx = mcorr_tbl.Number == dep_nums(i);
        plot(mcorr_tbl.r2_ice(dep_idx),mcorr_tbl.r4_ice(dep_idx),'o','markerfacecolor',colors(i))
    end
end
%% 

Tbeta_seg = Tbeta_func(mcorr_tbl.r2_ice,mcorr_tbl.r4_ice);
dT_seg = dT_func(mcorr_tbl.r2_ice,mcorr_tbl.r4_ice);

Tbeta_seg_ci = cat(2,Tbeta_seg-dT_seg,Tbeta_seg+dT_seg);
beta_seg = atand(Tbeta_seg);
beta_seg_ci = atand(Tbeta_seg_ci);
dbeta_seg = diff(beta_seg_ci,1,2);

for i = 1:length(beta_seg)
    fprintf('%.1f (%.1f)\n',round(beta_seg(i),1),round(dbeta_seg(i),1))
end


