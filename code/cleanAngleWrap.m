function out = cleanAngleWrap(in)
% Function to minimize angle wrapping (in a heading time series, e.g.).
%
% out = cleanAngleWrap(in)
%
% KJW
% 27 Sep 2022

cuts = min(in):max(in);
performance = nan*cuts;

for i = 1:length(cuts)
    tmp = in;
    if cuts(i) < 0
        tmp(tmp<cuts(i)) = tmp(tmp<cuts(i))+360;
    else
        tmp(tmp>=cuts(i)) = tmp(tmp>=cuts(i))-360;
    end
    performance(i) = max(abs(diff(tmp)));
end
[~,cut_idx] = min(performance);
cut_best = cuts(cut_idx);

out = in;
if cut_best<0
    out(out<cut_best) = out(out<cut_best)+360;
else
    out(out>=cut_best) = out(out>=cut_best)-360;
end