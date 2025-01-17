function main
clear all;
close all;

import java.awt.Robot;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import java.io.File;
%% Input Data
evalin('base', 'load Elec_Price.txt');
evalin('base', 'load Elec_SellPrice.txt');
evalin('base', 'load P_PV1.txt');
evalin('base', 'load P_PV2.txt');
evalin('base', 'load P_CL1.txt');
evalin('base', 'load P_NL1.txt');
evalin('base', 'load P_CL2.txt');
evalin('base', 'load P_NL2.txt');

% P_PV2 = P_PV2.*1;

%% Parameters
alpha = 1135.6;
beta = 301.03;
gamma = 0.39;
M_power = 20;

PV1_Max=30;
PV2_Max=40;

t = 1/4;                    % Time interval [h]
Num_var = 24/t;    % Number of variables
P_grid_max = 80;            % Maximum grid power [kW]
P_conv_max = 100;           % Maximum AC-DC converter power [kW]
Eff_conv = 0.95;            % Converter efficiency

P_BESS1_max = 10;           % Maximum BESS1 power [kW]
P_BESS2_max = 10;           % Maximum BESS2 power [kW]
Eff_BESS1 = 0.9;             % BESS efficiency
Eff_BESS2 = 0.95;             % BESS efficiency
SOC1_max = 90;              % Maximum SOC [%]
SOC1_min = 10;              % Minimum SOC [%]
SOC2_max = 90;              % Maximum SOC [%]
SOC2_min = 10;              % Minimum SOC [%]
CAP1 = 85;                 % Capacity 1 [kWh]
CAP2 = 100;                % Capacity 2 [kWh]
SOC1_init = 50;             % Initial SOC [%]
SOC2_init = 50;             % Initial SOC [%]

P_EV1_max = 11;           % ?? ???? [kW]
P_EV2_max = 11;           % ?? ???? [kW]
Eff_EV1 = 0.9;             % BESS(PCS*???) ??
Eff_EV2 = 0.95;             % BESS(PCS*???) ??
EV1SOC_max = 90;              % ?? ?? SOC [%]
EV1SOC_min = 10;              % ?? ?? SOC [%]
EV2SOC_max = 90;              % ?? ?? SOC [%]
EV2SOC_min = 10;              % ?? ?? SOC [%]
EV_CAP1 = 46;                 % ???? [kWh]
EV_CAP2 = 68;                % ???? [kWh]
EV1_SOCT = 80;
EV2_SOCT = 70;
EV1SOC_init = 30;
EV2SOC_init = 40;
Ta1=8;
Td1=17;
Ta2=7;
Td2=20;

Pen_P_CL1 = 1000*ones(Num_var,1);
Pen_P_NL1 = 500*ones(Num_var,1);
Pen_P_CL2 = 1000*ones(Num_var,1);
Pen_P_NL2 = 500*ones(Num_var,1);
Pen_PV1 = 100*ones(Num_var,1);
Pen_PV2 = 100*ones(Num_var,1);

%% GUI
fig = uifigure('Name', 'Modular EMS GUI', 'Position', [100, 100, 1100, 800]);

F_S=18;
set(fig, 'DefaultAxesFontSize', F_S);

% Display image above the figures
% Create an axes object
ax = axes('Parent', fig, 'Position', [0.3, 0.3, 0.4, 0.4]);
% Load and display the image
img = imshow('Legend.png', 'Parent', ax);
% Adjust the axes size to make the image bigger
set(ax, 'Units', 'pixels', 'Position', [300, 370, 500, 120]);

equipmentStatusPanel = uipanel(fig, 'Title', 'Components', 'Position', [40, 470, 170, 320], 'FontSize', F_S);

% Adjust positions of input fields
uilabel(fig, 'Position', [50, 730, 150, 22], 'Text', 'Grid:', 'FontSize', F_S);
gridStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 730, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 700, 150, 22], 'Text', 'PV AC:', 'FontSize', F_S);
pvACStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 700, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 670, 150, 22], 'Text', 'PV DC:', 'FontSize', F_S);
pvDCStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 670, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 640, 150, 22], 'Text', 'Micro-turbine:', 'FontSize', F_S);
cdgStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 640, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 610, 150, 22], 'Text', 'BESS AC:', 'FontSize', F_S);
bessACStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 610, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 580, 150, 22], 'Text', 'BESS DC:', 'FontSize', F_S);
bessDCStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 580, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 550, 150, 22], 'Text', 'EV AC:', 'FontSize', F_S);
evACStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 550, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 520, 150, 22], 'Text', 'EV DC:', 'FontSize', F_S);
evDCStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 520, 100, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [50, 490, 150, 22], 'Text', 'ILC:', 'FontSize', F_S);
ilcStatusField = uicheckbox(fig, 'Text', 'Yes', 'Position', [160, 490, 100, 22], 'FontSize', F_S);

calculateButton = uibutton(fig, 'Position', [50, 450, 150, 30], 'Text', 'Run Optimization', 'ButtonPushedFcn', @(btn, event) calculate_and_plot(), 'FontSize', F_S);



optimizationParametersPanel = uipanel(fig, 'Title', 'Parameters', 'Position', [220, 470, 830, 320], 'FontSize', F_S);
% Input field for P_grid_max
uilabel(fig, 'Position', [250, 730, 150, 22], 'Text', 'Power (kW):', 'FontSize', F_S);
PGridMaxField = uieditfield(fig, 'numeric', 'Value', P_grid_max, 'Position', [350, 730, 50, 22], 'FontSize', F_S);

% Panel to enclose "Optimal Cost" and result label
resultPanel = uipanel(fig, 'Title', '', 'Position', [650, 450, 250, 40], 'FontSize', F_S);

uilabel(resultPanel, 'Position', [10, 10, 150, 22], 'Text', 'Optimal Cost(₩):', 'FontSize', F_S);
resultLabel = uilabel(resultPanel, 'Position', [150, 10, 250, 22], 'Text', '', 'FontSize', F_S);

% Input field for Micro-turbine
uilabel(fig, 'Position', [250, 640, 150, 22], 'Text', 'γ(₩/kW2):', 'FontSize', F_S);
gammaField = uieditfield(fig, 'numeric', 'Value', gamma, 'Position', [350, 640, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 640, 150, 22], 'Text', 'β(₩/kW):', 'FontSize', F_S);
betaField = uieditfield(fig, 'numeric', 'Value', beta, 'Position', [520, 640, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [580, 640, 150, 22], 'Text', 'α(₩):', 'FontSize', F_S);
alphaField = uieditfield(fig, 'numeric', 'Value', alpha, 'Position', [670, 640, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [730, 640, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
powerField = uieditfield(fig, 'numeric', 'Value', M_power, 'Position', [830, 640, 50, 22], 'FontSize', F_S);

% Input field for PVs
uilabel(fig, 'Position', [250, 700, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
PV1maxField = uieditfield(fig, 'numeric', 'Value', PV1_Max, 'Position', [350, 700, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [250, 670, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
PV2maxField = uieditfield(fig, 'numeric', 'Value', PV2_Max, 'Position', [350, 670, 50, 22], 'FontSize', F_S);

% Input fields for AC BESS parameters
uilabel(fig, 'Position', [250, 610, 150, 22], 'Text', 'SoC L(%):', 'FontSize', F_S);
SOC1MinField = uieditfield(fig, 'numeric', 'Value', SOC1_min, 'Position', [350, 610, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 610, 150, 22], 'Text', 'SoC H(%):', 'FontSize', F_S);
SOC1MaxField = uieditfield(fig, 'numeric', 'Value', SOC1_max, 'Position', [520, 610, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [580, 610, 150, 22], 'Text', 'Cap(kWh):', 'FontSize', F_S);
CAP1Field = uieditfield(fig, 'numeric', 'Value', CAP1, 'Position', [670, 610, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [730, 610, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
PBESS1MaxField = uieditfield(fig, 'numeric', 'Value', P_BESS1_max, 'Position', [830, 610, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [900, 610, 150, 22], 'Text', 'Eff(%):', 'FontSize', F_S);
EffBESS1Field = uieditfield(fig, 'numeric', 'Value', Eff_BESS1, 'Position', [980, 610, 50, 22], 'FontSize', F_S);


% Input fields for DC BESS parameters
uilabel(fig, 'Position', [250, 580, 150, 22], 'Text', 'SoC L(%):', 'FontSize', F_S);
SOC2MinField = uieditfield(fig, 'numeric', 'Value', SOC2_min, 'Position', [350, 580, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 580, 150, 22], 'Text', 'SoC H(%):', 'FontSize', F_S);
SOC2MaxField = uieditfield(fig, 'numeric', 'Value', SOC2_max, 'Position', [520, 580, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [580, 580, 150, 22], 'Text', 'Cap(kWh):', 'FontSize', F_S);
CAP2Field = uieditfield(fig, 'numeric', 'Value', CAP2, 'Position', [670, 580, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [730, 580, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
PBESS2MaxField = uieditfield(fig, 'numeric', 'Value', P_BESS2_max, 'Position', [830, 580, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [900, 580, 150, 22], 'Text', 'Eff(%):', 'FontSize', F_S);
EffBESS2Field = uieditfield(fig, 'numeric', 'Value', Eff_BESS2, 'Position', [980, 580, 50, 22], 'FontSize', F_S);

% Input fields for AC EV parameters
uilabel(fig, 'Position', [250, 550, 150, 22], 'Text', 'SOC T(%):', 'FontSize', F_S);
EV1SOCTarField = uieditfield(fig, 'numeric', 'Value', EV1_SOCT, 'Position', [350, 550, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 550, 150, 22], 'Text', 'Arrive(h):', 'FontSize', F_S);
EV1ARRField = uieditfield(fig, 'numeric', 'Value', Ta1, 'Position', [520, 550, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [580, 550, 150, 22], 'Text', 'Depart(h):', 'FontSize', F_S);
EV1DEPField = uieditfield(fig, 'numeric', 'Value', Td1, 'Position', [670, 550, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [730, 550, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
EV1MaxField = uieditfield(fig, 'numeric', 'Value', P_EV1_max, 'Position', [830, 550, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [900, 550, 150, 22], 'Text', 'Eff(%):', 'FontSize', F_S);
EffEV1Field = uieditfield(fig, 'numeric', 'Value', Eff_EV1, 'Position', [980, 550, 50, 22], 'FontSize', F_S);


% Input fields for DC EV parameters
uilabel(fig, 'Position', [250, 520, 150, 22], 'Text', 'SoC T(%):', 'FontSize', F_S);
EV2SOCTarField = uieditfield(fig, 'numeric', 'Value', EV2_SOCT, 'Position', [350, 520, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 520, 150, 22], 'Text', 'Arrive(h):', 'FontSize', F_S);
EV2ARRField = uieditfield(fig, 'numeric', 'Value', Ta2, 'Position', [520, 520, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [580, 520, 150, 22], 'Text', 'Depart(h):', 'FontSize', F_S);
EV2DEPField = uieditfield(fig, 'numeric', 'Value', Td2, 'Position', [670, 520, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [730, 520, 150, 22], 'Text', 'Power(kW):', 'FontSize', F_S);
EV2MaxField = uieditfield(fig, 'numeric', 'Value', P_EV2_max, 'Position', [830, 520, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [900, 520, 150, 22], 'Text', 'Eff(%):', 'FontSize', F_S);
EffEV2Field = uieditfield(fig, 'numeric', 'Value', Eff_EV2, 'Position', [980, 520, 50, 22], 'FontSize', F_S);


% Input fields for AC-DC converter parameters
uilabel(fig, 'Position', [250, 490, 150, 22], 'Text', 'Cap(kW):', 'FontSize', F_S);
PConvMaxField = uieditfield(fig, 'numeric', 'Value', P_conv_max, 'Position', [350, 490, 50, 22], 'FontSize', F_S);

uilabel(fig, 'Position', [430, 490, 150, 22], 'Text', 'Eff(%):', 'FontSize', F_S);
EffConvField = uieditfield(fig, 'numeric', 'Value', Eff_conv, 'Position', [520, 490, 50, 22], 'FontSize', F_S);


% Move the graphs to a lower position to avoid overlap
ax3 = uiaxes(fig, 'Position', [50, 20, 500, 200]);
title(ax3, 'Generation');
ylabel(ax3, 'kW');
xlabel(ax3, 'Time of day (15 min)');

ax4 = uiaxes(fig, 'Position', [550, 20, 500, 200]);
title(ax4, 'Load');
ylabel(ax4, 'kW');
xlabel(ax4, 'Time of day (15 min)');


ax1 = uiaxes(fig, 'Position', [50, 220, 500, 200]);
title(ax1, 'SOC');
ylabel(ax1, '%');

ax2 = uiaxes(fig, 'Position', [550, 220, 500, 200]);
title(ax2, 'Load');
ylabel(ax2, 'kW');


    function calculate_and_plot()
      
        Elec_Price = evalin('base', 'Elec_Price');
        Elec_SellPrice = evalin('base', 'Elec_SellPrice');
        P_PV1 = evalin('base', 'P_PV1');
        P_PV2 = evalin('base', 'P_PV2');
        P_CL1 = evalin('base', 'P_CL1');
        P_NL1 = evalin('base', 'P_NL1');
        P_CL2 = evalin('base', 'P_CL2');
        P_NL2 = evalin('base', 'P_NL2');
        
        
        
        
        % Get input data from fields
        Grid_status = gridStatusField.Value;
        PV_AC_status = pvACStatusField.Value;
        PV_DC_status = pvDCStatusField.Value;
        CDG_status = cdgStatusField.Value;
        BESS_AC_status = bessACStatusField.Value;
        BESS_DC_status = bessDCStatusField.Value;
        EV_AC_status = evACStatusField.Value;
        EV_DC_status = evDCStatusField.Value;
        ILC_status = ilcStatusField.Value;
        
        % Get alpha, beta, gamma values
        alpha = alphaField.Value;
        beta = betaField.Value;
        gamma = gammaField.Value;
        M_power = powerField.Value;
        
        % Get PV max values
        PV1_Max = PV1maxField.Value;
        PV2_Max = PV2maxField.Value;
        
        % Get P_conv_max and Eff_conv values
        P_conv_max = PConvMaxField.Value;
        Eff_conv = EffConvField.Value;
        
        % Get AC BESS parameters
        P_BESS1_max = PBESS1MaxField.Value;
        Eff_BESS1 = EffBESS1Field.Value;
        SOC1_max = SOC1MaxField.Value;
        SOC1_min = SOC1MinField.Value;
        CAP1 = CAP1Field.Value;
        
        % Get DC BESS parameters
        P_BESS2_max = PBESS2MaxField.Value;
        Eff_BESS2 = EffBESS2Field.Value;
        SOC2_max = SOC2MaxField.Value;
        SOC2_min = SOC2MinField.Value;
        CAP2 = CAP2Field.Value;
        
                % Get AC EV parameters
                P_EV1_max = EV1MaxField.Value;
                Eff_EV1 = EffEV1Field.Value;
                EV1_SOCT = EV1SOCTarField.Value;
                Ta1 = EV1ARRField.Value;
                Td1 = EV1DEPField.Value;
        
                % Get DC EV parameters
                P_EV2_max = EV2MaxField.Value;
                Eff_EV2 = EffEV2Field.Value;
                EV2_SOCT = EV2SOCTarField.Value;
                Ta2 = EV2ARRField.Value;
                Td2 = EV2DEPField.Value;
        
        % Get Maximum Grid Power
        P_grid_max = PGridMaxField.Value;
        
        % Call different modules (Assuming these are defined elsewhere in your code)
        P_PV1 = PV_AC_Module(Num_var, P_PV1, PV1_Max, PV_AC_status);
        P_PV2 = PV_DC_Module(Num_var, P_PV2, PV2_Max, PV_DC_status);
        [BAC_lb, BAC_ub, BAC_A, BAC_b, SOC1_init] = BESS_AC_Module(t, Num_var, Eff_BESS1, P_BESS1_max,  CAP1, SOC1_init, SOC1_min, SOC1_max, BESS_AC_status);
        [BDC_lb, BDC_ub, BDC_A, BDC_b, SOC2_init] = BESS_DC_Module(t, Num_var, Eff_BESS2, P_BESS2_max, CAP2, SOC2_init, SOC2_min, SOC2_max, BESS_DC_status);
        [EAC_lb, EAC_ub, EAC_A, EAC_b, EV1SOC_init, AC_EVP] = EV_AC_Module(t, Num_var, Eff_EV1, P_EV1_max, EV_CAP1, EV1SOC_init, EV1SOC_min, EV1SOC_max, Ta1, Td1, EV1_SOCT, Grid_status, EV_AC_status);
        [EDC_lb, EDC_ub, EDC_A, EDC_b, EV2SOC_init, DC_EVP] = EV_DC_Module(t, Num_var, Eff_EV2, P_EV2_max, EV_CAP2, EV2SOC_init, EV2SOC_min, EV2SOC_max, Ta2, Td2, EV2_SOCT, Grid_status, ILC_status, EV_DC_status);
        [G_lb, G_ub] = Grid_Module(Num_var, P_grid_max, P_CL1, P_NL1, P_CL2, P_NL2, P_PV1, P_PV2, Grid_status, ILC_status);
        [S_diesel, C_lb, C_ub, C_A, C_b] = CDG_Module(alpha, beta, gamma, M_power, Num_var, CDG_status);
        [I_lb, I_ub] = ILC_Module(Num_var, Eff_conv, P_conv_max, ILC_status);
        
        f = t * [Elec_Price; Elec_SellPrice; Pen_P_CL1; Pen_P_NL1; Pen_P_CL2; Pen_P_NL2; Pen_PV1; Pen_PV2; alpha * ones(Num_var, 1); S_diesel(1) * ones(Num_var, 1); S_diesel(2) * ones(Num_var, 1); S_diesel(3) * ones(Num_var, 1); S_diesel(4) * ones(Num_var, 1); S_diesel(5) * ones(Num_var, 1); S_diesel(6) * ones(Num_var, 1); S_diesel(7) * ones(Num_var, 1); S_diesel(8) * ones(Num_var, 1); S_diesel(9) * ones(Num_var, 1); S_diesel(10) * ones(Num_var, 1); zeros(Num_var * 10, 1)];
        intcon = 1 + Num_var * 3 : Num_var * 4; % Integer variables
        
        % Lower Bounds
        lb = [G_lb; C_lb; BAC_lb; BDC_lb; EAC_lb; EDC_lb; I_lb];
        
        % Upper Bounds
        ub = [G_ub; C_ub; BAC_ub; BDC_ub; EAC_ub; EDC_ub; I_ub];
        
        % Inequality Constraints
        A = [BAC_A; BDC_A; EAC_A; EDC_A; C_A];
        b = [BAC_b; BDC_b; EAC_b; EDC_b; C_b];
        
        % Equality Constraints
        Aeq1 = [eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), zeros(Num_var,Num_var*2), -eye(Num_var), zeros(Num_var,Num_var*1), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var), zeros(Num_var,Num_var*2), eye(Num_var), eye(Num_var), zeros(Num_var,Num_var*2), -eye(Num_var)/Eff_conv, -eye(Num_var)*Eff_conv];
        beq1 = P_CL1 + P_NL1 - P_PV1;       % AC MG ? ??? = ???
        Aeq2 = [zeros(Num_var,Num_var*4), eye(Num_var), eye(Num_var), zeros(Num_var,Num_var*1), -eye(Num_var), zeros(Num_var,Num_var*13), eye(Num_var), eye(Num_var), zeros(Num_var,Num_var*2), eye(Num_var), eye(Num_var), eye(Num_var), eye(Num_var)];
        beq2 = P_CL2 + P_NL2 - P_PV2;       % DC MG ? ??? = ???
        Aeq = [Aeq1; Aeq2];
        beq = [beq1; beq2];
        
        % Solve Mixed-integer Linear Programming
        [x,fval] = linprog(f, A, b, Aeq, beq, lb, ub);
        
        % Extract solution variables
        P_grid = x(1 : Num_var) + x(1 + Num_var : Num_var*2);
        P_Shed_CL1 = x(1 + Num_var*2 : Num_var*3)
        P_Shed_NL1 = x(1 + Num_var*3 : Num_var*4)
        P_Shed_CL2 = x(1 + Num_var*4 : Num_var*5);
        P_Shed_NL2 = x(1 + Num_var*5 : Num_var*6);
        P_Cur_PV1 = x(1 + Num_var*6 : Num_var*7);
        P_Cur_PV2 = x(1 + Num_var*7 : Num_var*8);

        AC_Side=sum(P_Shed_CL1)+sum(P_Shed_NL1)
        DC_Side=sum(P_Shed_CL2)+sum(P_Shed_NL2)

        AC_PV=sum(P_Cur_PV1)
        DC_PV=sum(P_Cur_PV2)
        
        P_diesel = x(1 + Num_var*8 : Num_var*9) + x(1 + Num_var*9 : Num_var*10) + x(1 + Num_var*10 : Num_var*11) + x(1 + Num_var*11 : Num_var*12) + x(1 + Num_var*12 : Num_var*13) + x(1 + Num_var*13 : Num_var*14) + x(1 + Num_var*14 : Num_var*15) + x(1 + Num_var*15 : Num_var*16) + x(1 + Num_var*16 : Num_var*17) + x(1 + Num_var*17 : Num_var*18) + x(1 + Num_var*18 : Num_var*19);
        P_BESS1dis = x(1 + Num_var*19 : Num_var*20);
        P_BESS1chg = x(1 + Num_var*20 : Num_var*21);
        P_BESS1 = P_BESS1dis + P_BESS1chg;
        P_BESS2dis = x(1 + Num_var*21 : Num_var*22);
        P_BESS2chg = x(1 + Num_var*22 : Num_var*23);
        P_BESS2 = P_BESS2dis + P_BESS2chg;
        
        P_EV1dis = x(1 + Num_var*23 : Num_var*24);
        P_EV1chg = x(1 + Num_var*24 : Num_var*25);
        P_EV1 = P_EV1dis + P_EV1chg;
        
        P_EV2dis = x(1 + Num_var*25 : Num_var*26);
        P_EV2chg = x(1 + Num_var*26 : Num_var*27);
        P_EV2 = P_EV2dis + P_EV2chg;
        
        P_conv = x(1 + Num_var*27 : Num_var*28) + x(1 + Num_var*28 : Num_var*29);
        SOC1 = SOC1_init*ones(Num_var,1) - (t*100/CAP1)*tril(ones(Num_var))*P_BESS1dis/Eff_BESS1 - (t*100/CAP1)*tril(ones(Num_var))*P_BESS1chg*Eff_BESS2;
        SOC2 = SOC2_init*ones(Num_var,1) - (t*100/CAP2)*tril(ones(Num_var))*P_BESS2dis/Eff_BESS2 - (t*100/CAP2)*tril(ones(Num_var))*P_BESS2chg*Eff_BESS2;
        EV1SOC = EV1SOC_init*ones(Num_var,1) - (t*100/EV_CAP1)*tril(ones(Num_var))*P_EV1dis/Eff_EV1 - (t*100/EV_CAP1)*tril(ones(Num_var))*P_EV1chg*Eff_EV1;
        EV2SOC = EV2SOC_init*ones(Num_var,1) - (t*100/EV_CAP2)*tril(ones(Num_var))*P_EV2dis/Eff_EV2 - (t*100/EV_CAP2)*tril(ones(Num_var))*P_EV2chg*Eff_EV2;
        
        EV1SOC=EV1SOC.*AC_EVP;
        EV2SOC=EV2SOC.*DC_EVP;
        
        resultLabel.Text = num2str(fval);
        display (num2str(fval))
        % Plot results
%         cla(ax1);
%         hold(ax1, 'on');
%         bar(ax1, [P_PV1, P_BESS1, P_EV1, -P_conv, P_grid, P_diesel],   'stacked');
%         hold(ax1, 'on');
%         plot(ax1, P_CL1+P_NL1,'k.');
%         
%         legend(ax1, {'PV', 'BESS', 'EV', 'ILC', 'Grid', 'Micro', 'Load'}, 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');
%         title(ax1, 'AC MG Power Balance');
%         ylabel(ax1, 'Power (kW)');
        % Define the colors for each component
colors = [
    [0.4660 0.6740 0.1880];   % PV - Green
    [0      0.4470 0.7410];   % BESS - Blue
    [0.8500 0.3250 0.0980];   % EV - Orange
    [0.9290 0.6940 0.1250];   % ILC - Yellow
    [1      1      1];        % Grid - Whilte
    [0.6    0.6    0.6];      % Micro - Gray
];

% Data for the bar graph
data = [P_PV1, P_BESS1, P_EV1, -P_conv, P_grid, P_diesel];

% Clear and hold the axes
cla(ax1);
hold(ax1, 'on');

% Create the bar graph
hBar = bar(ax1, data, 'stacked', 'FaceColor', 'flat');

% Apply the colors
for k = 1:length(hBar)
    hBar(k).CData = repmat(colors(k, :), size(hBar(k).YData, 1), 1);
end

% Hold the plot
hold(ax1, 'on');

% Plot additional data
plot(ax1, P_CL1+P_NL1,'k.','MarkerSize', 12);

% Add legend
% legend(ax1, {'PV', 'BESS', 'EV', 'ILC', 'Grid', 'Micro', 'Load'}, 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');

% Add title and labels
title(ax1, 'Power balance: AC side');
ylabel(ax1, 'Power (kW)');

%         cla(ax2);
%         hold(ax2, 'on');
%         bar(ax2, [P_PV2, P_BESS2, P_EV2, P_conv], 'stacked');
%         hold(ax2, 'on');
%         plot(ax2, P_CL2+P_NL2,'k.');
%         
%         %         legend(ax2, 'P_PV2', 'P_BESS2', 'P_conv', 'P_Load');
%         title(ax2, 'DC MG Power Balance');
%         ylabel(ax2, 'Power (kW)');
% Clear and hold the second axes
cla(ax2);
hold(ax2, 'on');

% Create the second bar graph
hBar2 = bar(ax2, [P_PV2, P_BESS2, P_EV2, P_conv], 'stacked', 'FaceColor', 'flat');

% Apply the colors to the second bar graph
% Note: Match indices to ensure the same quantities have the same colors
hBar2(1).CData = repmat(colors(1, :), size(hBar2(1).YData, 1), 1); % PV
hBar2(2).CData = repmat(colors(2, :), size(hBar2(2).YData, 1), 1); % BESS
hBar2(3).CData = repmat(colors(3, :), size(hBar2(3).YData, 1), 1); % EV
hBar2(4).CData = repmat(colors(4, :), size(hBar2(4).YData, 1), 1); % ILC

% Hold the plot
hold(ax2, 'on');

% Plot additional data for the second graph
plot(ax2, P_CL2 + P_NL2, 'k.','MarkerSize', 12);

% Add legend, title, and labels for the second graph
% legend(ax2, {'PV', 'BESS', 'EV', 'ILC', 'Load'}, 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');
title(ax2, 'Power balance: DC side');
ylabel(ax2, 'Power (kW)');        
        cla(ax3);
        plot(ax3, SOC1);
        hold(ax3, 'on');
        plot(ax3, SOC2);
        hold(ax3, 'on');
        plot(ax3, EV1SOC);
        hold(ax3, 'on');
        plot(ax3, EV2SOC);
        legend(ax3, {'BESS_1', 'BESS_1','EV_1', 'EV_2'}, 'Location', 'northwest', 'Orientation', 'horizontal', 'Box', 'off');
        title(ax3, 'SOC: BESS and EV');
        ylabel(ax3, 'SOC (%)');
        ylim(ax3,[0,130]);
        
        cla(ax4);
        plot(ax4, P_Shed_CL1);
        hold(ax4, 'on');
        plot(ax4, P_Shed_NL1);
        plot(ax4, P_Shed_CL2);
        plot(ax4, P_Shed_NL2);
        plot(ax4, P_Cur_PV1);
        plot(ax4, P_Cur_PV2);
        ylim(ax4,[0,20]);
        legend(ax4, {'C_1', 'N_1', 'C_2', 'N_2', 'P_1', 'P_2'}, 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');
        title(ax4, 'PV and Load Curtailments');
        ylabel(ax4, 'Power (kW)');
        
    end

end
