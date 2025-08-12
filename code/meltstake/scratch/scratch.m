roll = -7;
pitch = 4;

CR = cos(roll*pi/180);
SR = sin(roll*pi/180);
CP = cos(pitch*pi/180);
SP = sin(pitch*pi/180);

V = eye(5);

       %   (u   w  v5 v1-4  err)
Rroll = [   1   0   0   0   0;... % u
            0  CR  SR   0   0;... % w'
            0 -SR  CR   0   0;... % v5'
            0 -SR   0  CR   0;... % v1-4'
            0   0   0   0   1];   % err

        %   (u   w  v5 v1-4  err)
Rpitch = [  CP  SP   0   0   0;... % u'
           -SP  CP   0   0   0;... % w'
             0   0   1   0   0;... % v5
             0   0   0   1   0;... % v1-4
             0   0   0   0   1];   % err

R = Rpitch*Rroll;