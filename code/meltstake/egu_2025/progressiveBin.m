function [Pbin,fbin,M] = progressiveBin(f,P,nf,nf_stops)
% function to perform progressive binning of PSDs, inputs should be column
% vectors (f and P)
%
% KJW
% 12 Apr 2025

% useful values
n_stops = length(nf);
df = diff(f(1:2));
% preallocate
Pbin = [];
fbin = [];
M = []; % degrees of freedom

% flip columns to rows
f = f';
P = P';

% loop through bin widths
nf_stops(end+1) = Inf;
for i = 1:n_stops
    idxf = f>=nf_stops(i) & f<nf_stops(i+1);
    [Pi,fi] = time_bin(f(idxf)',P(:,idxf),nf(i)*df);
    Mi = 2*nf(i)*ones(size(fi));
    Pbin = [Pbin; Pi'];
    fbin = [fbin; fi'+0.5*diff(fi(1:2))];
    M = [M; Mi'];
end
