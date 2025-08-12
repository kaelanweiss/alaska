% Script to find some nice CTD casts from the RHIB to provide context for
% our iceberg measurements and Eric's model
%
% KJW
% 3 Apr 2023
clear

% load
ctd_file = 'G:\Shared drives\Ice-ocean-interactions\data\LeConte2208\data\RHIB\proc\CTD\0824\compiled_casts_0824.nc';
finfo = ncinfo(ctd_file);
ctd = struct;
for i = 1:length(finfo.Variables)
    ctd.(finfo.Variables(i).Name) = ncread(ctd_file,finfo.Variables(i).Name);
end
ctd.datetime = datetime(ctd.time,'ConvertFrom','epochtime');

% N2
rho = ctd.rho;
rho0 = mean(rho,'all','omitnan');
for i = 1:size(rho,2)
    rho(:,i) = hannFilter(rho(:,i),7);
end
g = 9.8;
dz = [diff(ctd.z); nan];
drho = nan*rho;
drho(2:end-1,:) = rho(3:end,:)-rho(1:end-2,:);
drdz = drho./repmat(2*dz,[1 size(rho,2)]);
N2 = -(g/rho0)*drdz;
for i = 1:size(rho,2)
    N2(:,i) = hannFilter(N2(:,i),1);
end

gps_dir = dir('F:\alaska2022\data\RHIB\mat\polly\*.mat');
for i = 1:length(gps_dir)
    gps(i) = load(fullfile(gps_dir(i).folder,gps_dir(i).name));
end

load F:\alaska2022\data\iceberg_surveys\proc\20220824_singingflower\dragon\profiles.mat

meanT = mean(ctd.CT,2,'omitnan');
meanS = mean(ctd.SP,2,'omitnan');
meanTS = [meanT meanS];

%% figures
figure(1); clf
clear ax

zmax = -14;
idxz = ctd.z>=zmax;

flds = {'CT','SP'};
flds_drag = {'temperature','salinity'};
ls = {'-','--'};
pnums = [1 11];
for i = 1:2
    ax(i) = axes(figure(1),'position',axgridpos(1,2,i,0.08,0.08,.1,.1));
    hold on
    plot(ctd.(flds{i})(idxz,:),ctd.z(idxz))
    grid on
    box on
    xlabel(flds_drag{i})
    ylabel('z (m)')
    for j = 1:length(pnums)
        plot(profs(pnums(j)).ctd.(flds_drag{i}),-profs(pnums(j)).ctd.depth,'k','linewidth',0.8,'linestyle',ls{j})
    end
    plot(meanTS(idxz,i),ctd.z(idxz),'k-','linewidth',2);
    ylim([zmax 0])
end

