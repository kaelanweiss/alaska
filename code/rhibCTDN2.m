% script to calculate buoyancy frequency from rhib ctd casts
% 
% KJW
% 17 May 2023

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

%%
t_cast = min(ctd.datetime);
idxz = ctd.z >= -7 & ctd.z <= -3;
idxt = t_cast >= datetime('25-Aug-2022 00:42:00');

%% 
figure(1); clf
subplot(1,2,1)
plot(rho(idxz,idxt)-1000,ctd.z(idxz))
xlabel('density (kg/m^3)')
ylabel('depth (m)')

subplot(1,2,2)
plot(N2(idxz,idxt),ctd.z(idxz))
xlabel('N^2')


