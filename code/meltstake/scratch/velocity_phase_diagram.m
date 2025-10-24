% Script to play around with plotting BL velocity as a phase diagram
%
% that was fun, looks cool, probably not useful
%
% KJW
% 23 Oct 2025

clear
addpath ..

% load
dep_tbl = loadMSInfo(26:28);
dep_num = 27;

load(fullfile('F:meltstake/data/raw',dep_tbl.Folder{dep_tbl.Number==dep_num},'adcp','adcp.mat'))

seg_tbl = loadMSInfo(dep_num,'segments');

%% time and space indexing
idxt = adcp.burst.time>=seg_tbl.Start(1) & adcp.burst.time<=seg_tbl.End(end);
idxr = adcp.burst.range<=min(seg_tbl.rmax);
time = adcp.burst.time(idxt);

% outer velocity time series
U_raw = squeeze(mean(adcp.burst.vel_ice(idxt,idxr,[1 2]),2,'omitnan'));
U_raw = U_raw(:,1) + 1i*U_raw(:,2);

% filter a little bit
k_hann = 1*adcp.burst.samplerate;
U = hannFilter(U_raw,k_hann);


figure(1); clf
hold on; box on
plot(time,real(U_raw),'linewidth',1)
plot(time,imag(U_raw),'linewidth',1)
plot(time,real(U),'linewidth',1)
plot(time,imag(U),'linewidth',1)

% time derivative
dU_dt = [nan; U(3:end)-U(1:end-2); nan]/(2/adcp.burst.samplerate);

% forcing timescale
T_f = U./dU_dt;

% rotary forcing timescale
phi = unwrap(angle(U));
dphi_dt = [nan; phi(3:end)-phi(1:end-2); nan]/(2/adcp.burst.samplerate);
dphi_dt_no0 = dphi_dt;
dphi_dt_no0(abs(dphi_dt)<0.01) = nan;

%% phase diagrams
figure(2); clf
clear ax
for i = 1:8
    ax(i) = axes(figure(2),'position',axgridpos(4,2,i,[.05 .05 .08 .08],[0 0]));
    hold on; box on
end

plot(ax(1),real(U),real(dU_dt),'k.')
plot(ax(2),imag(U),real(dU_dt),'k.')

plot(ax(3),real(U),imag(dU_dt),'k.')
plot(ax(4),imag(U),imag(dU_dt),'k.')

plot(ax(5),real(U),abs(dU_dt),'k.')
plot(ax(6),imag(U),abs(dU_dt),'k.')

plot(ax(7),real(U),angle(dU_dt),'k.')
plot(ax(8),imag(U),angle(dU_dt),'k.')

xlabel(ax(7),'u')
xlabel(ax(8),'w')

ylabel(ax(1),'real(dU/dt)')
ylabel(ax(3),'imag(dU/dt)')
ylabel(ax(5),'abs(dU/dt)')
ylabel(ax(7),'angle(dU/dt)')

