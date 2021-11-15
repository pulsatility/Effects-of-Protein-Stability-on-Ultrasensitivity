%----------Reset the system for each run----------------------------%
clear all
clc
close all

tic

%% -------------------------- PARAMETERS ----------------------------------%%
param.Xtot = 0;
param.Ytot = 1;
param.k0 = 1;
param.k1f = 10;
param.k1b = 90;
param.k1c = 10;
param.k2f = 10;
param.k2b = 90;
param.k2c = 10;
param.k3 = 0.01;
param.k4 = 0.01;

Km1 = (param.k1b+param.k1c)/param.k1f;
Km2 = (param.k2b+param.k2c)/param.k2f;

%Re-assign to temporary parameters
Ytot = param.Ytot;
k0 = param.k0;
k1f = param.k1f;
k1b = param.k1b;
k1c = param.k1c;
k2f = param.k2f;
k2b = param.k2b;
k2c = param.k2c;
k3 = param.k3;
k4 = param.k4;

%% ------------------------- INITIAL CONDITION-----------------------------%%
init_R = 100;
init_Rp = 0;
init_RX = 0; 
init_RpY = 0;
y0 = [init_R, init_Rp, init_RX, init_RpY];

%% ------------------------ DOSE RESPONSE SIMULATIONS-----------------------%%
doseincrement= 1.01

tmin=0;
tmax=100000;
tspan = [tmin:tmax];

LRC_R_MinArray = [];
HillCoef_RArray = [];

LRC_RRX_MaxArray = [];
HillCoef_RRXArray = [];

LRC_Rp_MaxArray = [];
HillCoef_RpArray = [];

LRC_RpRpY_MaxArray = [];
HillCoef_RpRpYArray = [];

LRC_Rtotfree_MaxMinArray = [];
HillCoef_RtotfreeArray = [];

LRC_Rtot_MaxMinArray = [];
HillCoef_RtotArray = [];


%%Loop through different fold changes for a parameter (k4 as an example here)
parameter_name = genvarname('k4');
for n = [0.1,0.33,1,3.3,10]
    n
    param.k4 = n*k4;
    
    param.Xtot = 0;
    
    Xdose= [];
    RssArray= [];
    RpssArray= [];
    RXssArray= [];
    RpYssArray= [];
    RtotssArray= [];
    
    %% ------------------------- RUN SIMULATION--------------------------------%%
    while param.Xtot <= 1e4
        % ------------------------- RUN SIMULATION--------------------------------%
        [t,y] = ode15s('stabilization_ode',tspan, y0, [], param); %Note: use @stabilization_ode does not work for passing argument
        
        %Simulation result;
        R    = y(:,1);
        Rp   = y(:,2);
        RX   = y(:,3);
        RpY  = y(:,4); 
        
        Rtot = R + Rp + RX + RpY;
  
        % ------------------------- PLOT TIME-COURSE RESULTS----------------------------------%%
        
%             figure(1)
%             plot(t,R);
%             xlabel('time');
%             ylabel('R');
%             hold on
% 
%             figure(2)
%             plot(t,RX);
%             xlabel('time');
%             ylabel('RX');
%             hold on
%         
%             figure(3)
%             plot(t,Rp);
%             xlabel('time');
%             ylabel('Rp');
%             hold on
% 
%             figure(4)
%             plot(t,RpY);
%             xlabel('time');
%             ylabel('RpY');
%             hold on
%         
%             figure(5)
%             plot(t,R+Rp);
%             xlabel('time');
%             ylabel('Rfreetot');
%             hold on
%             
%             figure(6)
%             plot(t,Rtot);
%             xlabel('time');
%             ylabel('Rtot');
%             hold on
        
        %Record steady-state values
        RssArray= [RssArray,R(end)];
        RpssArray= [RpssArray,Rp(end)];
        RXssArray= [RXssArray,RX(end)];
        RpYssArray= [RpYssArray,RpY(end)];
        RtotssArray= [RtotssArray,Rtot(end)];
        
        Xdose= [Xdose,param.Xtot];
        
        if param.Xtot == 0
            param.Xtot = 0.1;
        else
            param.Xtot = param.Xtot*doseincrement;
        end
    end  %END OF WHILE LOOP
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%Plot Dose Response%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xmin = 0.1;
    xmax = 1e4;
    
    figure(40)
    loglog(Xdose,RssArray,'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rss');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    figure(50)
    loglog(Xdose,RssArray+RXssArray, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rss+RXss');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    figure(60)
    loglog(Xdose,RpssArray, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rpss');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    figure(70)
    loglog(Xdose,RpssArray+RpYssArray, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rpss+RpYss');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    figure(80)
    loglog(Xdose,RssArray+RpssArray, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rtotfreess');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    figure(90)
    loglog(Xdose,RtotssArray, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rtotss');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for R%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of R
    LRC=[];
    for i= 2:1:length(RssArray)   %skipping the first dose where X = 0
        if i==length(RssArray)    %stop when reaching the last dose becasue no next dose available to calculate the delta
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
    
    figure(100)
    lineR = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('R Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
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
    
    figure(100)
    hline = refline(0,HillCoef_R);
    hline.Color = get(lineR, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for R+RX%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of R+RX
    LRC=[];
    Rss_RXssArray = RssArray+RXssArray;
    for i= 2:1:length(Rss_RXssArray)   %skipping the first dose where X = 0
        if i==length(Rss_RXssArray)    %stop when reaching the last dose becasue no next dose available to calculate the delta
            break
        end
        delta_R = Rss_RXssArray(i+1) - Rss_RXssArray(i);
        PerInc_R = delta_R/Rss_RXssArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    LRC_R_Min = min(LRC)
    LRC_R_MinArray = [LRC_R_MinArray,LRC_R_Min];
    
    figure(105)
    lineRRX = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('R+RX Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    %Hill coefficient of R+RX
    Rmax = max(Rss_RXssArray);
    Rmin = min(Rss_RXssArray);
    R10 = Rmax - (Rmax-Rmin)*0.9;
    R90 = Rmax - (Rmax-Rmin)*0.1;
    
    for i= length(Rss_RXssArray):-1:1
        if Rss_RXssArray(i) >= R10
            i
            X10=Xdose(i);
            break
        end
    end
    
    for i= length(Rss_RXssArray):-1:1
        if Rss_RXssArray(i) >= R90
            X90 = Xdose(i);
            break
        end
    end
    
    HillCoef_RRX = log(81)/log(X90/X10)
    HillCoef_RRXArray = [HillCoef_RRXArray, HillCoef_RRX];
    
    figure(105)
    hline = refline(0,HillCoef_RRX);
    hline.Color = get(lineRRX, 'Color');
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
    
    figure(110)
    lineRp = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rp Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
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
    
    figure(110)
    hline = refline(0,HillCoef_Rp);
    hline.Color = get(lineRp, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for Rp+RpY%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of Rp+RpY
    LRC=[];
    Rpss_RpYssArray = RpssArray+RpYssArray;
    for i= 2:1:length(Rpss_RpYssArray)   %skipping the first dose where X = 0
        if i==length(Rpss_RpYssArray)    %stop when reaching the last dose becasue no next dose available to calculate the delta
            break
        end
        delta_R = Rpss_RpYssArray(i+1) - Rpss_RpYssArray(i);
        PerInc_R = delta_R/Rpss_RpYssArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    LRC_R_Min = min(LRC)
    LRC_R_MinArray = [LRC_R_MinArray,LRC_R_Min];
    
    figure(115)
    lineRpRpY = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rp+RpY Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    
    %Hill coefficient of Rp+RpY
    Rpmax = max(Rpss_RpYssArray);
    Rp10 = Rmax * 0.1;
    Rp90 = Rmax * 0.9;
    
    for i= length(Rpss_RpYssArray):-1:1
        if (Rpss_RpYssArray)(i) >= Rp10
            X10=Xdose(i);
            break
        end
    end
    
    for i= length(Rpss_RpYssArray):-1:1
        if (Rpss_RpYssArray)(i) >= Rp90
            X90 = Xdose(i);
            break
        end
    end
    
    HillCoef_RpRpY = log(81)/log(X90/X10)
    HillCoef_RpRpYArray = [HillCoef_RpRpYArray, HillCoef_RpRpY];
    
    figure(115)
    hline = refline(0,HillCoef_RpRpY);
    hline.Color = get(lineRpRpY, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate and plot LRC and Hill coefficient for Rtotfree = R+Rp%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LRC of Rtotfree
    LRC=[];
    RtotfreessArray = RssArray + RpssArray;
    for i= 2:1:length(RtotfreessArray)   %skipping the first dose where X = 0
        if i==length(RtotfreessArray)    %stop when reaching the last dose becasue no next dose avaiable to calculate the delta
            break
        end
        delta_R=RtotfreessArray(i+1)-RtotfreessArray(i);
        PerInc_R=delta_R/RtotfreessArray(i);
        
        delta_X=Xdose(i+1)-Xdose(i);
        PerInc_X=delta_X/Xdose(i);
        
        LRC(i-1) = PerInc_R/PerInc_X;
    end
    
    figure(130)
    lineRtotfree = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rtotfree Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    
    %Hill coefficient of Rtotfree
    [Rtotfreemax,RtotfreemaxIndex] = max(RtotfreessArray);
    [Rtotfreemin,RtotfreeminIndex] = min(RtotfreessArray);
    Rtotfree10 = (Rtotfreemax-Rtotfreemin)* 0.1 + Rtotfreemin;
    Rtotfree90 = (Rtotfreemax-Rtotfreemin)* 0.9 + Rtotfreemin;
    
    %For monotonic increasing curve
    if RtotfreemaxIndex > RtotfreeminIndex
        
        for i= 1:1:length(RtotfreessArray)
            if RtotfreessArray(i) >= Rtotfree10
                X10=Xdose(i);
                break
            end
        end

        for i= 1:1:length(RtotfreessArray)
            if RtotfreessArray(i) >= Rtotfree90
                X90 = Xdose(i);
                break
            end
        end
        
        LRC_Rtotfree_Max = max(LRC)
        LRC_Rtotfree_MaxMinArray = [LRC_Rtotfree_MaxMinArray,LRC_Rtotfree_Max];
        
    end
    
    %For monotonic decreasing curve
    if RtotfreemaxIndex < RtotfreeminIndex
        
        for i= length(RtotfreessArray):-1:1
            if RtotfreessArray(i) >= Rtotfree10
                X10=Xdose(i);
                break
            end
        end

        for i= length(RtotfreessArray):-1:1
            if RtotfreessArray(i) >= Rtotfree90
                X90 = Xdose(i);
                break
            end
        end
        
        LRC_Rtotfree_Min = min(LRC)
        LRC_Rtotfree_MaxMinArray = [LRC_Rtotfree_MaxMinArray,LRC_Rtotfree_Min];
        
    end
    
    HillCoef_Rtotfree = log(81)/log(X90/X10)
    HillCoef_RtotfreeArray = [HillCoef_RtotfreeArray, HillCoef_Rtotfree];
    
    figure(130)
    hline = refline(0,HillCoef_Rtotfree);
    hline.Color = get(lineRtotfree, 'Color');
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
    
    figure(140)
    lineRtot = semilogx(Xdose(2:(end-1)),LRC, 'LineWidth',2);
    xlabel('Xtot');
    ylabel('Rtot Response Coefficient');
    hold on
    box on
    xlim([xmin xmax])
    pbaspect([1.5 1 1])
    set(gca,'Fontsize',18);
    
    
    %Hill coefficient of Rtot
    [Rtotmax,RtotmaxIndex] = max(RtotssArray);
    [Rtotmin,RtotminIndex] = min(RtotssArray);
    Rtot10 = (Rtotmax-Rtotmin)* 0.1 + Rtotmin;
    Rtot90 = (Rtotmax-Rtotmin)* 0.9 + Rtotmin;
    
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
    
    figure(140)
    hline = refline(0,HillCoef_Rtot);
    hline.Color = get(lineRtot, 'Color');
    hline.LineStyle = '--';
    hline.LineWidth = 2;
 
end %END OF n FOR LOOP

toc