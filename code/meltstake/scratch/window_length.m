T = (1:(1/60):240)'; % record length (min)

T1 = 20;
T2 = 30;

n = nan*T;
L = nan*T;

for i = 1:length(T)

    if T(i) <= T1
        n(i) = 1;
        L(i) = T(i);
        continue
    end
    
    nj = floor(2*T(i)/T2-1):ceil(2*T(i)/T1-1);
    Lj = 2*T(i)./(nj+1);
    idx = find((Lj>=T1) & (Lj<=T2),1,'first');

    if ~isempty(idx)
        n(i) = nj(idx);
        L(i) = Lj(idx);
    end

end

figure(1); clf
plot(T,L,'.-')
yyaxis right
plot(T,n,'.-')
sum(isnan(n))