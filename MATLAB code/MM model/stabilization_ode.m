function dydt = stabilization_ode(t, y, option, param)


%% -------------------------- PARAMETERS MAPPING ----------------------------------%%

X   = param.X;
Y = param.Y;
k0 = param.k0;
k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
k4 = param.k4;
Km1 = param.Km1;
Km2 = param.Km2;

%% ------------------------- STATE NAME MAPPING----------------------------%%

R     = y(1);
Rp     = y(2);

%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%R
dydt(1) = k0 - k1*X*R/(Km1+R) + k2*Y*Rp/(Km2+Rp) - k3*R;

%Rp
dydt(2)	=  k1*X*R/(Km1+R) - k2*Y*Rp/(Km2+Rp) - k4*Rp;

%

end