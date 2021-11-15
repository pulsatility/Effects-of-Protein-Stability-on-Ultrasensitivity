%----------Reset the system for each run----------------------------%
clear all
clc
close all

tic

%% -------------------------- PARAMETERS ----------------------------------%%
param.X = 0;        %Kinase (concentration)
param.Y = 1;        %Phosphatase (concentration)
param.k0	= 1;    %Rate constant of synthesis of R (concentration/time)
param.k1	= 10;   %Catalytic rate constant for phosphorylation (1/time)
param.k2	= 10;   %Catalytic rate constant for dephosphorylation (concentration/time)
param.k3	= 0.01; %Degradation rate constant of R (1/time)
param.k4	= 0.01; %Degradation rate constant of Rp(1/time)
param.Km1 = 10;     %Michaelis constant for phosphorylation (concentration)
param.Km2 = 10;     %Michaelis constant for dephosphorylation

%Re-assign to temporary parameters
k0 = param.k0;
k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
k4 = param.k4;
Km1 = param.Km1;
Km2 = param.Km2;

%% ------------------------- INITIAL CONDITION -----------------------------%%
init_R = 100;
init_Rp = 0;
y0 = [init_R, init_Rp];

%% ------------------------ DOSE RESPONSE SIMULATIONS-----------------------%%
doseincrement= 1.01;

tmin=0;
tmax=100000;
tspan = [tmin:tmax];

LRC_R_MinArray = [];
HillCoef_RArray = [];
LRC_Rp_MaxArray = [];
HillCoef_RpArray = [];
LRC_Rtot_MaxMinArray = [];
HillCoef_RtotArray = [];


%%Loop through different fold changes for a parameter (k4 as an example here)
parameter_name = genvarname('k4');
for n = [0.1,0.33,1,3.3,10]
    n
    param.k4 = n*k4;
    
    param.X = 0;
    
    Xdose= [];
    RssArray= [];
    RpssArray= [];
    RtotssArray= [];
    flux_k1_Array = [];
    flux_k2_Array = [];
    flux_k3_Array = [];
    flux_k4_Array = [];
    
    
    while param.X <= 1000
        % ------------------------- RUN SIMULATION--------------------------------%
        [t,y] = ode15s('stabilization_ode',tspan, y0, [], param); %Note: use @stabilization_ode does not work for passing argument
        
        %Simulation result;
        R    = y(:,1);
        Rp   = y(:,2);
        
        Rtot = R + Rp;
        
        % ------------------------- PLOT TIME-COuRSE RESULTS----------------------------------%%     
%         figure(1)
%         plot(t,R,'LineWidth',2);
%         xlabel('time');
%         ylabel('R');
%         hold on
%         box on
%         pbaspect([1.5 1 1])
%         set(gca,'Fontsize',28);
%         
%         figure(2)
%         plot(t,Rp,'LineWidth',2);
%         xlabel('time');
%         ylabel('Rp');
%         hold on
%         box on
%         pbaspect([1.5 1 1])
%         set(gca,'Fontsize',28);
%         
%         figure(3)
%         plot(t,Rtot,'LineWidth',2);
%         xlabel('time');
%         ylabel('Rtot');
%         hold on
%         box on
%         pbaspect([1.5 1 1])
%         set(gca,'Fontsize',28);
        
        %Record steady-state values
        RssArray= [RssArray,R(end)];
        RpssArray= [RpssArray,Rp(end)];
        RtotssArray= [RtotssArray,Rtot(end)];
        
        Xdose= [Xdose,param.X];
        
        flux_k1 = param.X*param.k1*R(end)/(param.Km1+R(end));
        flux_k1_Array = [flux_k1_Array, flux_k1];
        
        flux_k2 = param.Y*param.k2*Rp(end)/(param.Km2+Rp(end));
        flux_k2_Array = [flux_k2_Array, flux_k2];
        
        flux_k3 = param.k3*R(end);
        flux_k3_Array = [flux_k3_Array, flux_k3];
        
        flux_k4 = param.k4*Rp(end);
        flux_k4_Array = [flux_k4_Array, flux_k4];
        
        if param.X == 0
            param.X = 0.01;
        else
            param.X = param.X*doseincrement;
        end
        
    end %END OF WHILE LOOP
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%Plot Dose Response%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    loglog(Xdose,RssArray,'LineWidth',2);
    xlabel('X');
    ylabel('Rss');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    figure(5)
    loglog(Xdose,RpssArray,'LineWidth',2);
    xlabel('X');
    ylabel('Rpss');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    figure(6)
    loglog(Xdose,RtotssArray,'LineWidth',2);
    xlabel('X');
    ylabel('Rtotss');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    figure(40)
    loglog(Xdose,flux_k1_Array,'LineWidth',2);
    xlabel('X');
    ylabel('flux k1');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);

    figure(50)
    loglog(Xdose,flux_k2_Array,'LineWidth',2);
    xlabel('X');
    ylabel('flux k2');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    figure(60)
    loglog(Xdose,flux_k3_Array,'LineWidth',2);
    xlabel('X');
    ylabel('flux k3');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    figure(70)
    loglog(Xdose,flux_k4_Array,'LineWidth',2);
    xlabel('X');
    ylabel('flux k4');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for R%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of R
    LRC=[];
    for i= 2:1:length(RssArray)   %skipping the first dose where X = 0
        if i==length(RssArray)    %stop when reaching the last dose becasue no next dose is available to calculate the delta
            break
        end
        delta_R=RssArray(i+1)-RssArray(i);
        PerInc_R=delta_R/RssArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    LRC_R_Min = min(LRC)
    LRC_R_MinArray = [LRC_R_MinArray,LRC_R_Min];
    
    figure(7)
    lineR = semilogx(Xdose(2:(end-1)),LRC,'LineWidth',2);
    xlabel('X');
    ylabel('R Response Coefficient');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    %Hill coefficient of R
    Rmax = max(RssArray);
    R10 = Rmax * 0.1;
    R90 = Rmax * 0.9;
   
    for i= length(RssArray):-1:1
        if RssArray(i) >= R10
            X10=Xdose(i);
            break
        end
    end
    
    for i= length(RssArray):-1:1
        if RssArray(i) >= R90
            X90 = Xdose(i);
            break
        end
    end
    
    HillCoef_R = log(81)/log(X90/X10)
    HillCoef_RArray = [HillCoef_RArray, HillCoef_R];
    
    figure(7)
    hline = refline(0,HillCoef_R);
    hline.Color = get(lineR, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for Rp%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of Rp
    LRC=[];
    for i= 2:1:length(RpssArray)   %skipping the first dose where X = 0
        if i==length(RpssArray)    %stop when reaching the last dose becasue no next dose avaiable to calculate the delta
            break
        end
        delta_R=RpssArray(i+1)-RpssArray(i);
        PerInc_R=delta_R/RpssArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    LRC_Rp_Max = max(LRC)
    LRC_Rp_MaxArray = [LRC_Rp_MaxArray,LRC_Rp_Max];
    
    figure(8)
    lineRp = semilogx(Xdose(2:(end-1)),LRC,'LineWidth',2);
    xlabel('X');
    ylabel('Rp Response Coefficient');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    
    %Hill coefficient of Rp
    Rpmax = max(RpssArray);
    Rp10 = Rpmax * 0.1;
    Rp90 = Rpmax * 0.9;
    
    for i= 1:1:length(RpssArray)
        if RpssArray(i) >= Rp10
            X10=Xdose(i);
            break
        end
    end
    
    for i= 1:1:length(RpssArray)
        if RpssArray(i) >= Rp90
            X90 = Xdose(i);
            break
        end
    end
    
    HillCoef_Rp = log(81)/log(X90/X10)
    HillCoef_RpArray = [HillCoef_RpArray, HillCoef_Rp];
    
    figure(8)
    hline = refline(0,HillCoef_Rp);
    hline.Color = get(lineRp, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for Rtot%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of Rtot
    LRC=[];
    for i= 2:1:length(RtotssArray)   %skipping the first dose where X = 0
        if i==length(RtotssArray)    %stop when reaching the last dose becasue no next dose avaiable to calculate the delta
            break
        end
        delta_R=RtotssArray(i+1)-RtotssArray(i);
        PerInc_R=delta_R/RtotssArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    figure(9)
    lineRtot = semilogx(Xdose(2:(end-1)),LRC,'LineWidth',2);
    xlabel('X');
    ylabel('Rtot Response Coefficient');
    hold on
    box on
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',28);
    

    %Hill coefficient of Rtot
    [Rtotmax,RtotmaxIndex] = max(RtotssArray);
    [Rtotmin,RtotminIndex] = min(RtotssArray);
    Rtot10 = (Rtotmax-Rtotmin) * 0.1 + Rtotmin;
    Rtot90 = (Rtotmax-Rtotmin) * 0.9 + Rtotmin;
    
    %For monotonic increasing curve
    if RtotmaxIndex > RtotminIndex
        for i= 1:1:length(RtotssArray)
            if RtotssArray(i) >= Rtot10
                X10=Xdose(i);
                break
            end
        end

        for i= 1:1:length(RtotssArray)
            if RtotssArray(i) >= Rtot90
                X90 = Xdose(i);
                break
            end
        end
        
        LRC_Rtot_Max = max(LRC)
        LRC_Rtot_MaxMinArray = [LRC_Rtot_MaxMinArray,LRC_Rtot_Max];
        
    end
    
    %For monotonic decreasing curve
    if RtotmaxIndex < RtotminIndex
        for i= length(RtotssArray):-1:1
            if RtotssArray(i) >= Rtot10
                X10=Xdose(i);
                break
            end
        end

        for i= length(RtotssArray):-1:1
            if RtotssArray(i) >= Rtot90
                X90 = Xdose(i);
                break
            end
        end
        
        LRC_Rtot_Min = min(LRC)
        LRC_Rtot_MaxMinArray = [LRC_Rtot_MaxMinArray,LRC_Rtot_Min];
        
    end
    
    HillCoef_Rtot = log(81)/log(X90/X10)
    HillCoef_RtotArray = [HillCoef_RtotArray, HillCoef_Rtot];
    
    figure(9)
    hline = refline(0,HillCoef_Rtot);
    hline.Color = get(lineRtot, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
 
end %END OF n FOR LOOP
    
toc