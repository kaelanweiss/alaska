function [val,sigma] = msTable2Vector(col)
% Function to convert a table column formatted as "mean (std)" to usable
% numeric values.
%
% [val,sigma] = msTable2Vector(col)
%
% KJW
% 6 Feb 2024

% preallocate
n = length(col);
val = nan(n,1);
sigma = nan(n,1);

% loop through rows
for i = 1:n
    stri = col{i};
    if isempty(stri)
        val(i) = nan;
        sigma(i) = nan;
    else
        C = strsplit(stri,' ');
        val(i) = str2double(C{1});
        sigma(i) = str2double(C{2}(2:end-1));
    end
end