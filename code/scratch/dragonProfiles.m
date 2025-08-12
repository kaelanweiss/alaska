% script to organize, plot Dragon profiles from Xeitl Sit' 2023

clear

%deps = {'20230921_180214','20230922_163549','20230925_193312'};
deps = {'20230922_163549','20230925_193312'};
dir_raw = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\Leconte2309\data\raw\dragon';
dir_proc = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\Leconte2309\data\processed\dragon';

% load adcp and ctd
ndeps = length(deps);
for i = 1:ndeps
    load(fullfile(dir_raw,deps{i},'adcp','adcp.mat'))
    load(fullfile(dir_raw,deps{i},'ctd','ctd.mat'))
    if ~isfield(adcp,'bt')
        adcp.bt = struct();
    end
    eval(sprintf('adcp%d=adcp;ctd%d=ctd(1);',i,i))
end

ctd = eval(['[' sprintf('ctd%d ',1:ndeps) '];']);
adcp = eval(['[' sprintf('adcp%d ',1:ndeps) '];']);
eval(['clear ' sprintf('ctd%d ',1:ndeps)])
eval(['clear ' sprintf('adcp%d ',1:ndeps)])

%% plot timeseries of pressure
for i = 1:ndeps
    figure(i); clf; hold on
    plot(ctd(i).time,ctd(i).values(:,strcmpi(ctd(i).channels,'pressure'))-9.5,'k-')
    plot(adcp(i).attitude.time,adcp(i).attitude.pressure,':',color=colors(1))
    datetick('x','mmmdd HH:MM','keeplimits')
    grid
    axis ij
end

% times of dragon casts
            % 22 Sept
dcasts = {{{'22-Sep-2023 17:54:00','22-Sep-2023 17:57:28'},...
           {'22-Sep-2023 18:15:47','22-Sep-2023 18:19:19'},...
           {'22-Sep-2023 22:00:10','22-Sep-2023 22:01:19'},...
           {'22-Sep-2023 22:32:42','22-Sep-2023 22:35:50'},...
           {'22-Sep-2023 22:49:59','22-Sep-2023 22:53:36'}},...
           ... % 25 Sept
          {{'25-Sep-2023 20:24:51','25-Sep-2023 20:25:21'},...
           {'25-Sep-2023 20:34:57','25-Sep-2023 20:36:08'},...
           {'25-Sep-2023 20:38:12','25-Sep-2023 20:40:37'},...
           {'25-Sep-2023 20:42:48','25-Sep-2023 20:44:45'},...
           {'25-Sep-2023 20:46:39','25-Sep-2023 20:53:08'},...
           {'25-Sep-2023 20:57:17','25-Sep-2023 21:00:28'},...
           {'25-Sep-2023 21:07:16','25-Sep-2023 21:17:45'},...
           {'25-Sep-2023 21:22:57','25-Sep-2023 21:25:35'},...
          }};

            % 22 Sept
ucasts = {{{'22-Sep-2023 18:06:30','22-Sep-2023 18:10:10'},...
           {'22-Sep-2023 18:19:27','22-Sep-2023 18:23:11'},...
           {'22-Sep-2023 22:01:27','22-Sep-2023 22:04:07'},...
           {'22-Sep-2023 22:36:14','22-Sep-2023 22:42:46'},...
           {'22-Sep-2023 22:53:45','22-Sep-2023 22:57:17'}},...
           ... % 25 Sept
          {{'25-Sep-2023 20:25:45','25-Sep-2023 20:27:17'},...
           {'25-Sep-2023 20:36:10','25-Sep-2023 20:37:22'},...
           {'25-Sep-2023 20:40:39','25-Sep-2023 20:42:19'},...
           {'25-Sep-2023 20:44:47','25-Sep-2023 20:45:40'},...
           {'25-Sep-2023 20:53:06','25-Sep-2023 20:55:39'},...
           {'25-Sep-2023 21:00:28','25-Sep-2023 21:04:06'},...
           {'25-Sep-2023 21:18:38','25-Sep-2023 21:22:45'},...
           {'25-Sep-2023 21:24:48','25-Sep-2023 21:26:54'},...
          }};

nprofs = zeros(length(dcasts),1);
for i = 1:length(nprofs)
    nprofs(i) = length(dcasts{i})+length(ucasts{i});
end

%% fall rate, heading rate
for i = 1:ndeps
    % ctd
    dt = mean(diff(ctd(i).time)*86400);
    ctd(i).dzdt = [diff(ctd(i).values(:,strcmpi(ctd(i).channels,'pressure')))/dt; 0];
    % adcp
    dt = adcp(i).burst.samplerate^-1;
    adcp(i).attitude.dzdt = [diff(adcp(i).attitude.pressure)/dt; 0];

    dhdt =  [diff(adcp(i).attitude.heading)/dt; 0];
    adcp(i).attitude.dhdt = [diff(cleanAngleWrap(adcp(i).attitude.heading))/dt; 0];

end

for i = 1:ndeps
    figure(i+2); clf; hold on
    plot(ctd(i).time,hannFilter(ctd(i).dzdt,17),'k-')
    plot(adcp(i).attitude.time,hannFilter(adcp(i).attitude.dzdt,7),'-',color=colors(1))
    datetick('x','mmmdd HH:MM','keeplimits')
    grid

    figure(i+4); clf; hold on
    plot(adcp(i).attitude.time,adcp(i).attitude.dhdt,'-',color=colors(1))
    datetick('x','mmmdd HH:MM','keeplimits')
    grid
end

%% organize data into profs
profs(sum(nprofs)) = struct();
pnum = 1;
t_init = nan(length(profs),1);

for i = 1:length(nprofs) % loop through files
    casts = cat(2,dcasts{i},ucasts{i});
    for j = 1:length(casts)
        t1 = datenum(casts{j}{1});
        t2 = datenum(casts{j}{2});

        t_init(pnum) = t1;

        %%% adcp %%%
        adcpj = struct('cfg',adcp(i).cfg,'units',adcp(i).units);
        % attitude
        idx = adcp(i).attitude.time >= t1 & adcp(i).attitude.time <= t2;
        adcpj.attitude = struct();
        for fld = fields(adcp(i).attitude)'
            fld = fld{1};
            sz = size(adcp(i).attitude.(fld));
            adcpj.attitude.(fld) = adcp(i).attitude.(fld)(idx,:);
        end
        
        % burst
        idx = adcp(i).burst.time >= t1 & adcp(i).burst.time <= t2;
        adcpj.burst = struct();
        for fld = fields(adcp(i).burst)'
            fld = fld{1};
            sz = size(adcp(i).burst.(fld));
            nt = sz(1);
            if nt ~= length(adcp(i).burst.time)
                adcpj.burst.(fld) = adcp(i).burst.(fld);
            else
                adcpj.burst.(fld) = adcp(i).burst.(fld)(idx,:,:);
            end
        end
        adcpj.burst.z = interp1(adcpj.attitude.time,adcpj.attitude.pressure,adcpj.burst.time);

        % bt
        if ~isempty(fields(adcp(i).bt))
            idx = adcp(i).bt.time >= t1 & adcp(i).bt.time <= t2;
            adcpj.bt = struct();
            for fld = fields(adcp(i).bt)'
                fld = fld{1};
                sz = size(adcp(i).bt.(fld));
                nt = sz(1);
                if nt ~= length(adcp(i).bt.time)
                    adcpj.bt.(fld) = adcp(i).bt.(fld);
                else
                    adcpj.bt.(fld) = adcp(i).bt.(fld)(idx,:);
                end
            end
            adcpj.bt.z = interp1(adcpj.attitude.time,adcpj.attitude.pressure,adcpj.burst.time);
        end

        %%% ctd %%%
        ctdj = struct();
        idx = ctd(i).time >= t1 & ctd(i).time <= t2;
        nt = length(ctd(i).time);
        for fld = fields(ctd(i))'
            fld = fld{1};
            if size(ctd(i).(fld),1) ~= nt
                ctdj.(fld) = ctd(i).(fld);
            else
                ctdj.(fld) = ctd(i).(fld)(idx,:);
            end
        end

        % bin by depth
        dz = 0.2;
        [binned, zbins] = time_bin(ctdj.values(:,strcmpi(ctdj.channels,'depth')),ctdj.values',dz);
        ctdj.zbins = zbins';
        ctdj.binned = binned';

        % prof
        profs(pnum).adcp = adcpj;
        profs(pnum).ctd = ctdj;

        % label up/down cast
        if diff(adcpj.attitude.pressure([1 end])) > 0
            profs(pnum).cast = 'down';
        else
            profs(pnum).cast = 'up';
        end

        pnum = pnum + 1;
    end
end

% sort by start time
[~,sidx] = sort(t_init);
profs = profs(sidx);
for i = 1:length(profs)
    profs(i).num = i;
end


%% plot ctd
for i = 1:length(profs)
    figure(300+i); clf;
    subplot(121)
    plot(profs(i).ctd.binned(:,strcmpi(profs(i).ctd.channels,'temperature')),profs(i).ctd.zbins)
    xlabel('T [C]')
    axis ij
    title(sprintf('profile %d',i))

    subplot(122)
    plot(profs(i).ctd.binned(:,strcmpi(profs(i).ctd.channels,'salinity')),profs(i).ctd.zbins)
    xlabel('S [psu]')
    axis ij
    title(datestr(profs(i).ctd.time(1)))
end

%% compare bt and dzdt (this actually works really well)
fnum = 401;
for i = 1:length(profs)
    if strcmpi(profs(i).cast,'up') || ~isfield(profs(i).adcp,'bt')
        continue
    end
    figure(fnum); clf; hold on
    fnum = fnum +1;

    idxf = profs(i).adcp.bt.fom>1000;
    vbt = profs(i).adcp.bt.vel;
    vbt(idxf) = nan;
    vbtz = (vbt(:,2) - vbt(:,4))/(2*sind(25));
    vbtx = (vbt(:,2) + vbt(:,4))/(2*cosd(25));
    zbt = profs(i).adcp.bt.z;
    plot(zbt,vbtz)

    plot(zbt,vbt(:,[1 3])/cosd(25))

    plot(zbt,vbtx,'k--')

    plot(profs(i).adcp.attitude.pressure,profs(i).adcp.attitude.dzdt,'k')
    
    profs(i).adcp.bt.z = zbt;

   %datetick('x','mmmdd HH:MM','keeplimits')
end

%% plot bt dist
fnum = 501;
for i = 1:length(profs)
    if strcmpi(profs(i).cast,'up') || ~isfield(profs(i).adcp,'bt')
        continue
    end
    figure(fnum); clf; hold on

    idxf = profs(i).adcp.bt.fom>1000;
    dbt = profs(i).adcp.bt.distance;
    dbt(idxf) = nan;
    zbt = interp1(profs(i).adcp.attitude.time,profs(i).adcp.attitude.pressure,profs(i).adcp.bt.time);
    plot(zbt,dbt)
    fprintf('%d: %d\n',fnum,i)
   %datetick('x','mmmdd HH:MM','keeplimits')
   fnum = fnum + 1;
end

good_profs = [15 19 23];
dist_max = [20 40 70];
profs2 = profs(good_profs);

%% adcp processing
cor_min = 50;
fom_min = 2000;

inst2ice = [  1   0   0   0   0;... % inst x      ice x
              0  -1   0   0   0;... % inst y      ice z 
              0   0   0   0   1;... % inst z1 --> ice y (b5)
              0   0 0.5 0.5   0;... % inst z2     ice y (b1-4)
              0   0   1  -1   0];   % inst b5     error 

ice_lbls = {'u (right)','w (up)','v b5 (inwrd)','v b1-4','err (1,3-2,4)'};
err_lbls = {'v1-v5','v2-v5','v1-v2'};

for i = 1:length(profs2)
    adcp = profs2(i).adcp;
    vel = adcp.burst.vel;

    % qa cor
    qa_cor = adcp.burst.cor<cor_min;
    vel(qa_cor) = nan;

    % qa cor
    vel_bt = adcp.bt.vel;
    dist_bt = adcp.bt.distance;
    qa_fom = adcp.bt.fom>fom_min;

    vel_bt(qa_fom) = nan;
    dist_bt(qa_fom) = nan;


    % interpolate bt values onto burst time grid
    t = adcp.burst.time;
    tbt = adcp.bt.time;

    vel_bt = interp1(tbt,vel_bt,t);
    dist_bt = interp1(tbt,dist_bt,t);

    % calculate beam 5 bt vel
    vel_bt5 = sum(vel_bt(:,[2 4]),2)/(2*cosd(25));
    vel_bt = [vel_bt vel_bt5];

    % combine burst and bt vel
    vel = permute(vel,[1 3 2]);
    vel = vel - repmat(vel_bt,[1 1 size(vel,3)]);
    
    % beam to instrument
    vel_inst = nan*vel;
    b2i = adcp.burst.beam2xyz;
    vel_inst(:,1,:) = b2i(1,1)*vel(:,1,:) + b2i(1,3)*vel(:,3,:);
    vel_inst(:,2,:) = b2i(2,2)*vel(:,2,:) + b2i(2,4)*vel(:,4,:);
    vel_inst(:,3,:) = b2i(3,1)*vel(:,1,:) + b2i(3,3)*vel(:,3,:);
    vel_inst(:,4,:) = b2i(4,2)*vel(:,2,:) + b2i(4,4)*vel(:,4,:);
    vel_inst(:,5,:) = vel(:,5,:);

    % instrument to ice
    vel_ice = nan*vel;
    vel_ice(:,1,:) = inst2ice(1,1)*vel_inst(:,1,:);
    vel_ice(:,2,:) = inst2ice(2,2)*vel_inst(:,2,:);
    vel_ice(:,3,:) = inst2ice(3,5)*vel_inst(:,5,:);
    vel_ice(:,4,:) = inst2ice(4,3)*vel_inst(:,3,:) + inst2ice(4,4)*vel_inst(:,4,:);
    vel_ice(:,5,:) = inst2ice(5,3)*vel_inst(:,3,:) + inst2ice(5,4)*vel_inst(:,4,:);

    vel = permute(vel,[1 3 2]);
    vel_inst = permute(vel_inst,[1 3 2]);
    vel_ice = permute(vel_ice,[1 3 2]);
    adcp.burst.vel_inst = vel_inst;
    adcp.burst.vel = vel;
    adcp.burst.vel_ice = vel_ice;

    % bin in z
    dz = 1;
    [vel_ice_binned,zbins] = time_bin(adcp.burst.z,permute(vel_ice,[2 3 1]),dz,'omitnan');
    vel_ice_binned = permute(vel_ice_binned,[3 1 2]);
    adcp.burst.vel_ice_binned = vel_ice_binned;
    adcp.burst.zbins = zbins';

    profs2(i).adcp = adcp;

    adcpQuickPlot(figure(600+i),adcp,'vel_inst',0.3*[-1 1],[0 Inf],[0 dist_max(i)],1);
    adcpQuickPlot(figure(650+i),adcp,'vel',0.25*[-1 1],[0 Inf],[0 dist_max(i)],1);
    adcpQuickPlot(figure(700+i),adcp,'vel_ice',0.25*[-1 1],[0 Inf],[0 dist_max(i)],1);

end
    
%% nice plots
% 5 panels, T, S, u, v, w
nf = 1;
mf = 5;
pdx = 0.02;
pdy = 0.05;
pdx_o = 0.1;
pdy_o = 0.1;

fs = 12;

bt_color = 'k';
bt_lw = 1.1;

T = cell(1,length(profs2));
S = cell(1,length(profs2));

for i = 1:length(profs2)
    clear ax
    figure(1000+i); clf
    
    % ctd
    % T
    ax(1) = axes(figure(1000+i),'position',axgridpos(nf,mf,1,pdx,pdy,pdx_o,pdy_o));
    plot(profs2(i).ctd.binned(:,2),profs2(i).ctd.zbins,'.-',color=colors(1))
    xlabel('T [C]',fontsize=fs)
    ylabel('depth (m)',fontsize=fs)
    axis ij
    grid
    T{i} = profs2(i).ctd.binned(:,2);
    title(sprintf('profile %d',profs2(i).num),fontsize=fs)

    % S
    ax(2) = axes(figure(1000+i),'position',axgridpos(nf,mf,2,pdx,pdy,pdx_o,pdy_o));
    plot(profs2(i).ctd.binned(:,5),profs2(i).ctd.zbins,'.-',color='r')
    xlabel('S [psu]',fontsize=fs)
    %ylabel('depth (m)',fontsize=fs)
    axis ij
    grid
    S{i} = profs2(i).ctd.binned(:,5);
    title([datestr(profs2(i).adcp.burst.time(1)) 'Z'],fontsize=fs)

    % adcp
    for j = 1:3
        ax(j+2) = axes(figure(1000+i),'position',axgridpos(nf,mf,j+2,pdx,pdy,pdx_o,pdy_o));
        pcolor(profs2(i).adcp.burst.range,profs2(i).adcp.burst.zbins,profs2(i).adcp.burst.vel_ice_binned(:,:,j))
        shading flat
        cmocean('bal')
        axis ij
        clim(0.25*[-1 1])
        xlabel('range [m]',fontsize=fs)
        title(ice_lbls{j},fontsize=fs)

        if j == 1
            hold on
            plot(profs2(i).adcp.bt.distance(:,1),profs2(i).adcp.bt.z,'-',linewidth=bt_lw,color=bt_color)
            plot(profs2(i).adcp.bt.distance(:,3),profs2(i).adcp.bt.z,'-',linewidth=bt_lw,color=bt_color)
        elseif j == 2
            hold on
            plot(profs2(i).adcp.bt.distance(:,2),profs2(i).adcp.bt.z,'-',linewidth=bt_lw,color=bt_color)
            plot(profs2(i).adcp.bt.distance(:,4),profs2(i).adcp.bt.z,'-',linewidth=bt_lw,color=bt_color)
        end
        
        xlim([1 dist_max(i)])
    end
    
    cbar = colorbar(ax(5));
    cbar.Position = [sum(ax(5).Position([1 3]))+0.005 ax(5).Position(2) 0.008 sum(ax(1).Position([2 4]))-ax(5).Position(2)];
    cbar.Label.String = 'vel [m/s]';
    cbar.Label.FontSize = fs;

    linkaxes(ax,'y')

end

%% TS plot
for i = 1:length(profs2)
    figure(2000+i); clf; hold on
    for j = 1:length(profs)
%         plot(S{j},T{j},'.',color=0.5*[1 1 1])
        plot(profs(j).ctd.binned(:,5),profs(j).ctd.binned(:,2),'.',color=0.5*[1 1 1])
    end
    plot(S{i},T{i},'k',linewidth=1)
    scatter(S{i},T{i},18,profs2(i).ctd.zbins,'filled')

    xlabel('S [psu]',fontsize=fs)
    ylabel('T [C]',fontsize=fs)
    title(sprintf('profile: %d | %sZ',profs2(i).num,datestr(profs2(i).adcp.burst.time(1))),fontsize=fs)
    grid
    box on
    cbar = colorbar;
    cbar.Label.String = 'depth [m]';
    cbar.Label.FontSize = fs;

    xlim([24 28.8])
    ylim([3.5 7.3])
    clim([0 90])

end


