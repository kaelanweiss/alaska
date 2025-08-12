% script to compute heat and salt transfer coefficients based on Kader and 
% Yaglom used in the 3-equation parameterization of Jenkins (1991)

clear

% constants
Cd = 2.5e-3;
Pr = 13.8;
Sc = 2432;

% range of Re
Re = 10.^(1:7);

% coeffs (divided by sqrt(Cd)*U)
GammaT = (2.12*log(sqrt(Cd)*Re) + 12.5*Pr^(2/3) - 8.68).^-1;
GammaS = (2.12*log(sqrt(Cd)*Re) + 12.5*Sc^(2/3) - 8.68).^-1;

% plot
% GammaT
figure(1); clf
subplot(2,2,1)
plot(Re,GammaT,'ko-')
xlabel('Re')
ylabel('\Gamma_T = \gamma_T/C_D^{1/2}U')
ylim([0 max(ylim)])
set(gca,'xscale','log')
% GammaS
subplot(2,2,2)
plot(Re,GammaS,'ko-')
xlabel('Re')
ylabel('\Gamma_S = \gamma_S/C_D^{1/2}U')
ylim([0 max(ylim)])
set(gca,'xscale','log')
% Thermal stanton
subplot(2,2,3)
plot(Re,sqrt(Cd)*GammaT,'ko-')
xlabel('Re')
ylabel('St_T = C_D^{1/2}\Gamma_T = \gamma_T/U')
ylim([0 max(ylim)])
set(gca,'xscale','log')
% Haline stanton
subplot(2,2,4)
plot(Re,sqrt(Cd)*GammaS,'ko-')
xlabel('Re')
ylabel('St_S = C_D^{1/2}\Gamma_S = \gamma_S/U')
ylim([0 max(ylim)])
set(gca,'xscale','log')