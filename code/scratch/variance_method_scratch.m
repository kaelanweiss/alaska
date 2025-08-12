clear

load F:/alaska2022/data/iceberg_surveys/mat/20220824_singingflower/spider/adcp.mat
load F:/alaska2022/data/iceberg_surveys/mat/20220824_singingflower/spider/adcp_exclude.mat
% load F:/AK_202305/adcp/103041_AK2/adcp.mat
% load F:/AK_202305/sections.mat

[Nt,nc,nb] = size(adcp.burst.vel);

% Coordinates:
%   x along beam 1 and -3 (right when facing ice)
%   y along beam 5 (positive into ice)
%   z along beam 2 and -4
%
% <u'v'> = C(<u1'^2>-<u3'^2>)
% <v'w'> = C(<u2'^2>-<u4'^2>)
% C = 2sin(25)cos(25) = 0.766

phi = 25;
C = 2*sind(phi)*cosd(phi);

% clean up velocity
cor_min = 70;
qa_cor = adcp.burst.cor<cor_min;
vel = adcp.burst.vel;
qa = qa_cor | repmat(exclude_idx,[1 nc nb]);
vel(qa) = nan;

% half overlapping mean windows
T = 10*60;
tsec = 86400*(adcp.burst.time-adcp.burst.time(1));
dt = mode(diff(tsec));
nt = 2*round(T/dt/2);
nw = fix((2*Nt-nt)/nt);

%% calculate mean and variance for each window
vel_mu = nan([nw nc nb]);
vel_var = nan([nw nc nb]);
vel_n = nan([nw nc nb]);
for i = 1:nw
    slice = (1:nt)+(i-1)*nt/2;
    vel_mu(i,:,:) = mean(vel(slice,:,:),'omitnan');
    vel_var(i,:,:) = var(vel(slice,:,:),'omitnan');
    vel_n(i,:,:) = sum(~isnan(vel(slice,:,:)));
end

uv = C*(vel_var(:,:,1)-vel_var(:,:,3));
vw = C*(vel_var(:,:,2)-vel_var(:,:,4));

figure(3); clf
for i = 1:nb
    % mean
    axes(figure(3),'position',axgridpos(nb,2,i,.05,.05,.05,.05,'flip'));
    plot(adcp.burst.range,vel_mu(:,:,i),'k');
    % std
    axes(figure(3),'position',axgridpos(nb,2,i+nb,.05,.05,.05,.05,'flip'));
    plot(adcp.burst.range,vel_var(:,:,i),'k');
end

figure(4); clf
subplot(1,2,1)
plot(adcp.burst.range,uv,'k')
xlim([0 0.8])
ylim(5e-4*[-1 1])
grid
title('u''v''')

subplot(1,2,2)
plot(adcp.burst.range,vw,'k')
xlim([0 0.8])
ylim(8e-4*[-1 1])
grid
title('v''w''')
