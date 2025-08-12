% Script to define realistic confidence intervals on melt rate
% observations. The Confidence intervals take into account:
%   1) linear regression
%   2) disagreement between beams
%   3) missing beams
%   4) ice obliquity
%   5) uncertainty from the tank experiment
%
% KJW
% 8 Apr 2024

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');
pc_tbl = loadMSInfo(1:19,'pitch_correction');

[m,m_ci] = msTable2Vector(ms_tbl.m);
[m1,m1_ci] = msTable2Vector(ms_tbl.m1);
[m2,m2_ci] = msTable2Vector(ms_tbl.m2);
[m3,m3_ci] = msTable2Vector(ms_tbl.m3);

% calculate differences between beams and set a characteristic difference
% to apply to missing beams (this may be double counting since the std of
% melt rate between beams is already incorporated in m_ci)
m_diff = abs([m1 m2 m3] - repmat(m,[1 3]));
m_diff_mean = mean(m_diff,'all','omitnan');
m_diff(isnan(m_diff)) = m_diff_mean;
m_diff_max = max(m_diff,[],2,'omitnan');

% set additional uncertainty parameters
ice_angle = 30;
tank_discrepancy_scale = .11; % see "tank_test_melt_rate.m"

% regression uncertainty, disagreement between beams (including missing beams), ice obliquity,
% uncertainty from tank experiment
m_ci = [sqrt(m_ci.^2 + ...
            m_diff_max.^2 + ...
            (m*(1-cosd(ice_angle))).^2 + ...
            (tank_discrepancy_scale*m).^2) ...
        sqrt(m_ci.^2 + ...
            m_diff_max.^2)];

for i = 1:size(m_ci,1)
    fprintf('%.2f (%.2f)\n',round(m_ci(i,1),2),round(m_ci(i,2),2))
end


