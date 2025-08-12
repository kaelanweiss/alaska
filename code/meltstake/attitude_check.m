% Script to check if meltstake is consitently melting upward while
% deployed. This is support work for the melt rate paper.
% 
% KJW
% 22 Jan 2024

clear

% collect adcp attitude data
% mstbl = readtable('G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx');
% depnames = table2cell(mstbl(:,3));
% empty_deps = isnan(mstbl.Number);
% depnames = depnames(~empty_deps);
ms_tbl = loadMSInfo;
depnames = ms_tbl.Folder;
ndeps = length(depnames);
raw_dir = 'F:/meltstake/data/raw';

att(ndeps) = struct;
for i = 1:ndeps
    fprintf('%d/%d\n',i,ndeps)
    load(fullfile(raw_dir,depnames{i},'adcp','adcp.mat'))
    flds = fieldnames(adcp.attitude);
    for j = 1:length(flds)
        att(i).(flds{j}) = adcp.attitude.(flds{j});
    end
    rolli = unwrap(att(i).roll*pi/180)*180/pi;
    rolli = (rolli-270)+360*(mean(rolli)<0);
    att(i).roll = rolli;
    att(i).time = datetime(att(i).time,'convertfrom','datenum');
end

%% smoothing and roll rate
diff_vec = [1 0 -1];
hw = 2*60; % sec
for i = 1:ndeps
    dt = seconds(diff(att(i).time([1 2])));
    att(i).roll_smth = hannFilter(att(i).roll,round(hw/dt));
    att(i).roll_rate = [nan; conv(att(i).roll_smth,diff_vec,'valid')/(2*dt); nan];
    att(i).roll_rate(1) = diff(att(i).roll_rate(1:2))/dt;
    att(i).roll_rate(end) = diff(att(i).roll_rate([end-1 end]))/dt;
end

%% some statistics
roll_diff = nan(ndeps,1);
roll_diff_dt = nan(ndeps,1);
for i = 1:ndeps
    roll_diff(i) = diff(att(i).roll([1 end]));
    roll_diff_dt(i) = roll_diff(i)/seconds(diff(att(i).time([1 end])));
end

%%
figure(1); clf
hold on
for i = 1:ndeps
    plot(att(i).time,att(i).roll)
end
grid on

figure(2); clf
hold on
for i = 1:ndeps
    ti = minutes(att(i).time-att(i).time(1));
    plot(ti,att(i).roll_smth)
    text(ti(1),att(i).roll_smth(1),num2str(i))
end
grid on

figure(3); clf
hold on
for i = 1:ndeps
    ti = minutes(att(i).time-att(i).time(1));
    plot(ti,att(i).roll_rate*3600)
end
grid on