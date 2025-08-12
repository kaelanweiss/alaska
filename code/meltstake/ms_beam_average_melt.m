% Script to convert individual beam melt rates with uncertainties to an
% average melt rate. Just for OSM results so this is a rough draft.
%
% KJW
% 6 Feb 2024

clear

ms_tbl = loadMSInfo(26:28,'manualwindows');
pc_tbl = loadMSInfo(26:28,'pitch_correction');

n = size(ms_tbl,1);
mi = nan(n,3);
mi_std = mi;

for i = 1:3
    [mi(:,i),mi_std(:,i)] = msTable2Vector(ms_tbl.(['m' num2str(i)]));
end

% average
m = mean(mi,2,'omitnan');

% add uncertainty in quadrature
m_ci = sqrt(sum(mi_std.^2,2,'omitnan')./sum(~isnan(mi_std),2).^2); % error prop through definition of mean
m_std = std(mi,0,2,'omitnan'); % stdev of measurements

m_ci_quad = sqrt(m_ci.^2 + 0*m_std.^2); % add error prop and stdev in quadrature

% apply correction due to meltstake pitch (if available)
[m_bias,dm_bias] = msTable2Vector(pc_tbl.melt_bias);
for i = 1:n
   if ~isnan(m_bias(i))
       m(i) = m(i) - m_bias(i);
       m_ci(i) = sqrt(m_ci(i)^2 + dm_bias(i)^2);
   end
end

% print
for i = 1:n
    fprintf('%.2f (%.2f)\n',round(m(i),2),round(m_ci_quad(i),2))
end
