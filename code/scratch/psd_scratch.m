% Script to test some spectral analysis tools and concepts
%
% KJW
% 23 Nov 2022

%% adcp data
try
    adcp;
catch
    load F:/alaska2022/data/iceberg_surveys/mat/20220824_teaparty/spider/adcp.mat
end
nt = length(adcp.burst.time);

beams = {'right','up','left','down'};

fs = round((86400*diff(adcp.burst.time(1:2)))^-1);
r_idx = 15;
t0_idx = 65000;
tf_idx = 93840;
vel = squeeze(adcp.burst.vel(t0_idx:tf_idx,r_idx,:));
vel_nan = vel; vel_nan(squeeze(adcp.burst.cor(t0_idx:tf_idx,r_idx,:))<50) = nan;
nfft = 2^nextpow2(size(vel,1));
if length(size(vel))>=3
    Svel = squeeze(mean(fft(detrend(vel.*repmat(hann(size(vel,1)),[1 size(vel,2)])),nfft,1),2));
else
    Svel = fft(detrend(vel.*repmat(hann(size(vel,1)),[1 size(vel,2)])),nfft,1);
end
%Svel = fft(vel,nfft,1);
fvel = fft_freq(1/8,nfft,nfft);
ifmax = find(fvel<=fs/2,1,'last');
Svel = Svel(2:ifmax,:);
fvel = fvel(2:ifmax);

figure(1); clf
for i = 1:4
    subplot(4,1,i)
    axlist(i) = gca;
    pcolor(adcp.burst.time(1:1:end),adcp.burst.range(1:60),adcp.burst.vel(1:1:end,1:60,i)')
    shading flat
    cmocean('balance')
    clim(.02*[-1 1])
    title(beams{i})
    datetick('x','keeplimits')
end
linkaxes(axlist,'x')
%xlim(extrema(adcp.burst.time(t0_idx:tf_idx)))

figure(2); clf; hold on
for i = 1:4
    %plot(adcp.burst.time(t0_idx:tf_idx),meanFilter(vel(:,i)-mean(vel(:,i)),1),'-')
    plot(adcp.burst.time(t0_idx:tf_idx),vel_nan(:,i))
end
title(sprintf('bin: %d, range: %.2fm',r_idx,adcp.burst.range(r_idx)))
xlim(extrema(adcp.burst.time(t0_idx:tf_idx)))
datetick('x','keeplimits')
grid
legend(beams)

figure(3); clf; hold on
for i = 1:size(Svel,2)  
    plot(fvel,meanFilter(abs(Svel(:,i)).^2,21).*(fvel'.^(0))*(10^(i*1.5)))
end
grid on
box on
xlabel('f')
ylabel('|S|^2 (unnormalized)')
title(sprintf('bin: %d, range: %.2fm',r_idx,adcp.burst.range(r_idx)))
legend(beams)
set(gca,'yscale','log','xscale','log')


%%% subfunctions %%%
function [S,f] = sum_ft(s,dt,n)
    N = length(s);
    t = cumsum(dt*ones(1,N))-dt;
    if n > N
        n = N;
    end
    f = (0:n-1)/(N*dt);
    S = nan*f;
    for k = 1:n
        hk = exp(-2*pi*1i*(k-1)*t/(N*dt));
        S(k) = sum(s.*hk);
    end
    S = (2/N)*S; % not sure about this normalization
end

function freq = fft_freq(dt,n,N)
    freq = (0:n-1)/(dt*N);
end