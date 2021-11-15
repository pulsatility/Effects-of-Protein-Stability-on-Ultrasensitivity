function dydt = stabilization_ode(t, y, option, param)


%% -------------------------- PARAMETERS MAPPING ----------------------------------%%

Xtot  =  param.Xtot;
Ytot  =  param.Ytot;
k0  =  param.k0;
k1f = param.k1f;
k1b = param.k1b;
k1c = param.k1c;
k2f = param.k2f;
k2b = param.k2b;
k2c = param.k2c;
k3 =  param.k3;
k4 =  param.k4;

%% ------------------------- STATE NAME MAPPING----------------------------%%

R   =   y(1);
Rp  =   y(2);
RX  =   y(3);
RpY =   y(4); 

X = Xtot - RX;
Y = Ytot - RpY;

%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%R
dydt(1) = k0 - k3 * R - k1f * X * R + k1b * RX + k2c * RpY ;

%Rp
dydt(2)	=  k1c * RX - k4 * Rp - k2f * Rp * Y + k2b * RpY  ;

%RX
dydt(3) = k1f * R * X - k1b * RX - k1c * RX - k3 * RX; 

%RpY
dydt (4) = k2f * Rp * Y - k2b * RpY - k2c * RpY - k4 * RpY;  
%

end