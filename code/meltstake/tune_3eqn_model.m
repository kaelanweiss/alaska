% Script to calculate the turbulent transfer coefficients in the 3-equation
% model based on the statistical relationship between observed and
% predicted melt rates. From what I can tell, the tuning can be applied
% after solving the 3-equation system because the system can be expressed
% independent of the Stanton number (the ratio of heat to salt transfer is
% what is specified). Therefore, the tuning statistics can be performed on
% the set of predicted melt rates without needing to resolve the system
% after each iteration.
%
% KJW
% 11 Sep 2026
clear

% load data and keep track of iceberg/glacier deps
ms_tbl = loadMSInfo('segments');
inc = ms_tbl.include;
idx_berg = ms_tbl.Number <= 19 & inc;
idx_glac = false(size(ms_tbl,1),3);
for i = 1:3
    idx_glac(:,i) = ms_tbl.Number == (i+25) & inc;
end

% load model results
mod_glac = load('../../data/melt_models_glacier.mat');
mod_berg = load('../../data/melt_models.mat');
mod_glac.m_mean = mod_glac.m_mean(any(idx_glac(ms_tbl.Start>datetime(2024,7,1),:),2),:);
mod_berg.m_mean = mod_berg.m_mean(idx_berg(ms_tbl.Start<datetime(2024,7,1)),:);

% cm/hr to m/dy
cmh2md = 0.24;

% extract melt rate and confidence intervals from spreadsheet
[m,~] = msTable2Vector(ms_tbl.m);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);

% convert melt rate to m/day
m = m*cmh2md;
m_ci_l = m_ci_l*cmh2md;
m_ci_h = m_ci_h*cmh2md;

% glacier melt rates
m_glac = m(inc & any(idx_glac,2));
m_glac_ci = m_ci_l(inc & any(idx_glac,2));
m_mod = mod_glac.m_mean(:,1);

% method 1: melt ratio with bootstrap (this scales the uncertainty by m_mod
% which is weird)
r = m(inc)./[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)];
r_ci = m_ci_l(inc)./[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)];
r_glac = m(inc & any(idx_glac,2))./mod_glac.m_mean(:,1);
r_glac_ci = m_ci_l(inc & any(idx_glac,2))./mod_glac.m_mean(:,1);
% [coeff, coeff_ci] = bootstrapMean(r_glac,.05,10000,'x_var',r_glac_ci);
n_bs = 1e4;
alpha = .05;
[coeff_r, coeff_r_ci] = bootstrapMean(r_glac,alpha,n_bs,'x_var',m_glac_ci.^2);

% linear least squares of a y = mx model
X = m_mod;
Y = m_glac;
W = diag(1./m_glac_ci.^2);
coeff_lsq = (X'*W*X)\(X'*W*Y);
coeff_lsq_sigma2 = inv(X'*W*X);

% print results
fprintf('Mean of ratios: %.2f [%.2f, %.2f]\n',round(coeff_r,2),round(coeff_r_ci(1),2),round(coeff_r_ci(2),2))
fprintf('Least squares: %.2f [%.2f, %.2f]\n',round(coeff_lsq,2),round(coeff_lsq-sqrt(coeff_lsq_sigma2),2),round(coeff_lsq+sqrt(coeff_lsq_sigma2),2))

% 
% % tuning
% T_stanton = 0.0011;
% S_stanton = 0.000031;
% r_glac = m(inc & any(idx_glac,2))./mod_glac.m_mean(:,1);
% r_glac_ci = m_ci_l(inc & any(idx_glac,2))./mod_glac.m_mean(:,1);
% [coeff, coeff_ci] = bootstrapMean(r_glac,.05,10000,'x_var',r_glac_ci);
% fprintf('Tuning: %.2f [%.2f, %.2f]\n',round(coeff,2),round(coeff_ci(1),2),round(coeff_ci(2),2))
% fprintf('Gamma_T*Cd^1/2: %.5f [%.5f, %.5f]\n',round(coeff*T_stanton,5),round(coeff_ci(1)*T_stanton,5),round(coeff_ci(2)*T_stanton,5))
% fprintf('Gamma_S*Cd^1/2: %.7f [%.7f, %.7f]\n',round(coeff*S_stanton,7),round(coeff_ci(1)*S_stanton,7),round(coeff_ci(2)*S_stanton,7))
