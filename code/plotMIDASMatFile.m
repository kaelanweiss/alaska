clear; 

% folder = 'C:/Users/ROSS/Documents/Alaska2022/Spider_QueenofHearts/Nortek/';
% fname = 'Data.235.00002.ad2cp.00000_';

% folder = 'E:\Leconte_2022\Sunsetberg\Dragon_20220822_0340\103041_Data.234.00000\';
% fname = 'Data.234.00000.ad2cp.00000_';

% folder = 'E:\Leconte_2022\Dragon_Oysterberg\103041_Data.234.00001\';
% fname = 'Data.234.00001.ad2cp.00000_';

% folder = 'C:\Users\ROSS\Documents\Alaska2022\Spider_test_20220824\';
% fname = 'Data.236.00000.ad2cp.00000_';

% folder = 'C:\Users\ROSS\Documents\Alaska2022\Spider_teapartyberg_20220824\Nortek\';
% fname = 'Data.236.00001.ad2cp.00000_';

folder = 'F:\Alaska2022\data\raw\iceberg_surveys\20220825_singingflower\dragon\Nortek\';
fname = 'polly_test_202.237.00000.ad2cp.00000_';

filenames = {};
for i = 1:4
    filenames{end+1} = [folder fname num2str(i) '.mat'];
end
data = cell(length(filenames),1);
config = cell(length(filenames),1);

% load files
nt = NaN(length(filenames),1);
fprintf('Loading\n');
for i = 1:length(filenames)
    fprintf('file %d/%d...\n',i,length(filenames));
    ld_vars = load(filenames{i});
    data{i} = ld_vars.Data;
    config{i} = ld_vars.Config;
    nt(i) = length(data{i}.BurstHR_MatlabTimeStamp);
end
nc = config{1}.bursthr_nCells;
nb = 4;
nc_echo = size(data{1}.Echo11000_0kHz_Amp1000_0kHz,2);
nt_tot = sum(nt);
nt_vec = cumsum([0; nt]);
%%
burst = struct();
burst.vel = NaN(nc,nt_tot,nb);
burst.cor = NaN(nc,nt_tot,nb);
burst.amp = NaN(nc,nt_tot,nb);
burst.hdg = NaN(nt_tot,1);
burst.pitch = NaN(nt_tot,1);
burst.roll = NaN(nt_tot,1);
burst.pressure = NaN(nt_tot,1);

% bt = struct();
% bt.vel = NaN(nt_tot,4);
% bt.dist = NaN(nt_tot,4);

echo = struct();
echo.amp = NaN(nt_tot,nc_echo);
echo.d = data{1}.Echo11000_0kHz_Range;

d = data{1}.BurstHR_Range;
time = NaN(nt_tot,1);

fprintf('Extracting\n');
for j = 1:length(filenames)
    fprintf('file %d/%d...\n',j,length(filenames));
    tlim = [nt_vec(j)+1 nt_vec(j+1)];
    time(tlim(1):tlim(2)) = data{j}.BurstHR_MatlabTimeStamp;
    % velocities
    for i = 1:4
        if i<5
            burst.vel(:,tlim(1):tlim(2),i) = data{j}.(['BurstHR_VelBeam' num2str(i)])';
            burst.cor(:,tlim(1):tlim(2),i) = data{j}.(['BurstHR_CorBeam' num2str(i)])';
            burst.amp(:,tlim(1):tlim(2),i) = data{j}.(['BurstHR_AmpBeam' num2str(i)])';
        else
            burst.vel(:,tlim(1):tlim(2),i) = data{j}.(['IBurst_VelBeam' num2str(i)])';
            burst.cor(:,tlim(1):tlim(2),i) = data{j}.(['IBurst_CorBeam' num2str(i)])';
            burst.amp(:,tlim(1):tlim(2),i) = data{j}.(['IBurst_AmpBeam' num2str(i)])';
        end
    end
    
    % orientation
    burst.hdg(tlim(1):tlim(2)) = data{j}.BurstHR_Heading;
    burst.pitch(tlim(1):tlim(2)) = data{j}.BurstHR_Pitch;
    burst.roll(tlim(1):tlim(2)) = data{j}.BurstHR_Roll;
    burst.pressure(tlim(1):tlim(2)) = data{j}.BurstHR_Pressure;

%     % bt
%     for i = 1:4
%         nt_bti = length(data{j}.(['BurstBT_VelBeam' num2str(i)]));
%         tlim2 = min([tlim(1)+nt_bti-1 tlim(2)]);
%         bt.vel(tlim(1):tlim2,i) = data{j}.(['BurstBT_VelBeam' num2str(i)]);
%         bt.dist(tlim(1):tlim2,i) = data{j}.(['BurstBT_DistanceBeam' num2str(i)]);
%     end
    
    % echo
    nt_eci = size(data{j}.Echo11000_0kHz_Amp1000_0kHz,1);
    tlim2 = min([tlim(1)+nt_eci-1 tlim(2)]);
    echo.amp(tlim(1):tlim2,:) = data{j}.Echo11000_0kHz_Amp1000_0kHz(1:(tlim2-tlim(1)+1),:);
end

burst.roll(burst.roll<0) = burst.roll(burst.roll<0) + 360;

%save([folder 'adcp.mat'],'burst','echo','-v7.3')

t_idx = burst.pressure>0;
%% plot vel
dlim = [0 0.8];
tlim = [datenum(2022,8,24,19,26,0) datenum(2022,8,24,19,37,0)];
clim = {.04*[-1 1],[40 100],[20 90]};
flds = {'vel','cor','amp'};
for i = 1:3
    ax = [];
    figure(i); clf;
    for j = 1:nb
        ax(j) = subplot(nb,1,j);
        pcolor(time(t_idx),d,burst.(flds{i})(:,t_idx,j));
        axis ij;
        shading flat;
        %xlim(tlim);
        ylim(dlim);
        if j==1
            title(flds{i});
        end
        colorbar;
        caxis(clim{i});
        datetick('x','mmmdd HH:MM','keeplimits');
        ylabel('depth [m]');
        if i==1
            cmocean('balance');
        elseif i==3
            cmocean('amp');
        end
    end
    linkaxes(ax,'x');
end

%% plot echo
figure(4); clf;
didx = echo.d>=0 & echo.d<=5;
pcolor(time(1:10:end),echo.d(didx),echo.amp(1:10:end,didx)');
%xlim(tlim);
datetick('x','HH:MM','keeplimits');
shading flat;
axis ij;
%ylim(dlim);
title('echosounder');
ylabel('depth [m]');
colorbar;

%% plot orientation
flds = {'hdg','pitch','roll','pressure'};
figure(5); clf;
ax = [];
for i = 1:4
    ax(i) = subplot(4,1,i);
    plot(time(t_idx),burst.(flds{i})(t_idx));
    %xlim(tlim);
    datetick('x','HH:MM','keeplimits');
    title(flds{i});
    grid;
end
linkaxes(ax,'x');

