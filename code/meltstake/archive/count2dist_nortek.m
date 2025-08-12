function SampleDist = count2dist_nortek(SampleNumber,SoundSpeed)
% COPYING FUNCTION FROM NORTEK FOLKS
% CountToDist(SampleNumber, SoundSpeed)    // mm
% {
%     VertDist = 25.3   // mm
%     HorzDist = 71.7   // mm
%     CountToDist = SoundSpeed/480.0
%     Dist2 = HorzDist*HorzDist + VertDist^2
%     MinDist = sqrt(Dist2)
%     Const = HorzDist^2 - VertDist^2
%  
%   TotalDist = CountToDist*SampleNumber
%     if (TotalDist > MinDist)
%         SampleDist = 0.5*(TotalDist^2 - Dist2)/(TotalDist - VertDist)
%  
%   return SampleDist;  // mm
% }
%  
% Regards,
% Nery

VertDist = 25.3;
HorzDist = 71.7;
CountToDist = SoundSpeed/480.0; % implies sample period of 1.1 us

Dist2 = HorzDist^2 + VertDist^2;
MinDist = sqrt(Dist2);
Const = HorzDist.^2 - VertDist.^2;

TotalDist = CountToDist.*SampleNumber;

SampleDist = 0.5*(TotalDist.^2 - Dist2)./(TotalDist - VertDist);% + sqrt(Const);

idx = ~(TotalDist > MinDist);
SampleDist(idx) = nan;