% Script to determine if the spectral peak in the u velocity during
% deployment 26 is actual flow or platform motion.
%
% KJW
% 5 Dec 2025
clear

% load
seg_tbl = loadMSInfo(28,'segments');
load(sprintf('F:/meltstake/data/raw/%s/adcp/adcp.mat',seg_tbl.Folder{1}))

%% indexing
seg_num = 10;
idx_att = adcp.attitude.time>=seg_tbl.Start(1) & adcp.attitude.time<=seg_tbl.End(end);
idx_vel = adcp.burst.time>=seg_tbl.Start(1) & adcp.burst.time<=seg_tbl.End(end);
idx_r = adcp.burst.range < min(seg_tbl.rmax);

% data
t_att = adcp.attitude.time(idx_att);
p = detrend(adcp.attitude.pressure(idx_att),4);
ptc = adcp.attitude.pitch_corrected(idx_att);
rol = adcp.attitude.roll_corrected(idx_att);

ptc(ptc<-180) = ptc(ptc<-180)+360;

t_vel = adcp.burst.time(idx_vel);
uvw = squeeze(mean(adcp.burst.vel_ice(idx_vel,idx_r,[1 3 2]),2,'omitnan'));
T = seconds(diff(extrema(t_vel)));
nt = sum(idx_vel);

% remove mean from velocity
for i = 1:3
    uvw(:,i) = uvw(:,i) - mean(uvw(:,i),'omitnan');
end

% fill in nans
idx_nan = isnan(uvw(:,1));
for i = 1:3
    uvw(idx_nan,i) = interp1(t_vel(~idx_nan),uvw(~idx_nan,i),t_vel(idx_nan),'linear');
end

% change in depth due to roll if right screw is stuck and left is free
L = 0.42; % (m) width between screw tips
H = 0.36; % (m) height between screw tips and ADCP pressure sensor
dz_rol = (L/2)*sind(rol) - H*(cosd(rol)-1);
eta = p+dz_rol; % corrected pressure
% calculate velocity effects ...

% PSDs
[Puvw,f] = psd_matlab(uvw,adcp.burst.samplerate,'notaper');
Pp = psd_matlab(eta,adcp.burst.samplerate,'notaper');
Pptc = psd_matlab(detrend(ptc),adcp.burst.samplerate,'notaper');
Prol = psd_matlab(rol,adcp.burst.samplerate,'notaper');

P = [Puvw Pp Pptc Prol]./repmat([var(uvw) var(eta) var(detrend(ptc)) var(rol)],[length(f) 1]);

%% cross spectrum of p and u
% cross correlations
k = 3000;
Rup = laggedCrossCorr(uvw(:,1),eta,k);
Rwp = laggedCrossCorr(uvw(:,3),eta,k);
Ruw = laggedCrossCorr(uvw(:,1),uvw(:,3),k);
Ruu = laggedCrossCorr(uvw(:,1),uvw(:,1),k);
Rww = laggedCrossCorr(uvw(:,3),uvw(:,3),k);
Rpp = laggedCrossCorr(eta,eta,k);

tau = (-k:k)/adcp.burst.samplerate;
fup = (-k:k)/(2*max(tau));

% FFTs of xcorrs
Fup = fftshift(fft(Rup))/length(tau);
Fwp = fftshift(fft(Rwp))/length(tau);
Fuw = fftshift(fft(Ruw))/length(tau);
Fuu = fftshift(fft(Ruu))/length(tau);
Fww = fftshift(fft(Rww))/length(tau);
Fpp = fftshift(fft(Rpp))/length(tau);

% FFTs of individual data
fS = ((-nt/2+1):(nt/2))'/T;
Fu = fftshift(fft(uvw(:,1)));
Fp = fftshift(fft(p));
Fw = fftshift(fft(uvw(:,3)));
Suu = conj(Fu).*Fu/T;
Sww = conj(Fw).*Fw/T;
Spp = conj(Fp).*Fp/T;
Sup = conj(Fu).*Fp/T;
Swp = conj(Fw).*Fp/T;
Suw = conj(Fu).*Fw/T;

Aup = abs(Fup); Phiup = atan2(-imag(Fup),real(Fup));
Awp = abs(Fwp); Phiwp = atan2(-imag(Fwp),real(Fwp));
Auw = abs(Fuw); Phiuw = atan2(-imag(Fuw),real(Fuw));

% coherence
% Cup = Aup.^2./(Fuu.*Fpp);
nb = 5;
[Supb,fb] = progressiveBin(fS,Sup,[nb],[0]);
Suub = progressiveBin(fS,Suu,[nb],[0]);
Sppb = progressiveBin(fS,Spp,[nb],[0]);
Cup = conj(Supb).*Supb./(Suub.*Sppb);

alpha = 0.01;
M = 2*nb;
qF = finv(1-alpha,2,2*M-2);
Gcrit = qF/(M-1+qF);
Gcrit = 1-alpha^(2/(M-2));


%% some plots
clear ax1
figure(1); clf
ax1(1) = axes(figure(1));
hold on
plot(t_att,p,'-','color',0.5*[1 1 1])
plot(t_att,dz_rol,'r-')
plot(t_att,eta,'k-')
ylabel('pressure [dbar]')

figure(2); clf
ax1(2) = axes(figure(2));
hold on
plot(t_att,ptc,'k-')
plot(t_att,rol,'r-')
legend({'pitch','roll'})
ylabel('pitch/roll [deg]')

figure(3); clf
ax1(3) = axes(figure(3));
plot(t_vel,uvw)
legend({'u','v','w'})
ylabel('vel [m/s]')

linkaxes(ax1,'x')

figure(4); clf
plot(f,P)
set(gca,'xscale','log','yscale','log')
legend({'u','v','w','pres','ptc','rol'})
ylabel('variance-normalized PSD')

figure(5); clf
hold on
plot(tau,Rup)
plot(tau,Rwp)
ylabel('cross-correlation')
legend({'u-pres','w-pres'})

figure(6); clf
hold on
plot(fup,real(Fup))
plot(fup,imag(Fup))
plot(fup,Aup,'k-')
xlim(0.1*[-1 1])
ylabel('u-pres cross spectra')

figure(7); clf
hold on
plot(fup,real(Fwp))
plot(fup,imag(Fwp))
plot(fup,Awp,'k-')
xlim(0.1*[-1 1])
ylabel('w-pres cross spectra')

figure(8); clf
hold on
plot(fup,real(Fuw))
plot(fup,imag(Fuw))
plot(fup,Auw,'k-')
xlim(0.1*[-1 1])
ylabel('u-w cross spectra')

figure(12); clf
hold on
plot(fup,Phiuw*180/pi,'color',colors(3))
plot(fup,Phiwp*180/pi,'color',colors(2))
plot(fup,Phiup*180/pi,'color',colors(1))
legend({'u-w','w-pres','u-pres'})
ylabel('angle spectrum')

figure(15); clf
hold on
plot(fb,Cup)
set(gca,'xscale','log')
ylabel('Coh^2')


