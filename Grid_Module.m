
function [G_lb, G_ub]= Grid_Module (Num_var,P_grid_max, P_CL1, P_NL1, P_CL2, P_NL2, P_PV1, P_PV2, Grid_status, ILC_status)

lb1 = zeros(Num_var,1);                    % P_gridBuy? ???? 0
lb2 = -P_grid_max*ones(Num_var,1);         % P_gridSell? ???? -P_grid_max
lb3 = zeros(Num_var,1);                    % AC side Level 1 load shedding
lb4 = zeros(Num_var,1);                    % AC side Level 2 load shedding
lb5 = zeros(Num_var,1);                    % DC side Level 1 load shedding
lb6 = zeros(Num_var,1);                    % DC side Level 2 load shedding
lb7 = zeros(Num_var,1);                    % AC PV curtailment
lb8 = zeros(Num_var,1);                    % DC PV curtailment

ub1 = P_grid_max*ones(Num_var,1);          % P_gridBuy? ???? P_grid_max
ub2 = zeros(Num_var,1);                    % P_gridSell? ???? 0
ub3 = zeros(Num_var,1);                    % AC side Level 1 load shedding
ub4 = zeros(Num_var,1);                    % AC side Level 2 load shedding
ub5 = zeros(Num_var,1);                    % DC side Level 1 load shedding
ub6 = zeros(Num_var,1);                    % DC side Level 2 load shedding
ub7 = zeros(Num_var,1);                    % AC PV curtailment
ub8 = zeros(Num_var,1);                    % DC PV curtailment

    
if (Grid_status ==0 || P_grid_max==0) % Islanded case no trading
    lb2 = zeros(Num_var,1);  
    ub1 = zeros(Num_var,1); 
    
    ub3 = P_CL1;                    % AC side Level 1 load shedding
    ub4 = P_NL1;                    % AC side Level 2 load shedding
    ub5 = P_CL2;                    % DC side Level 1 load shedding
    ub6 = P_NL2;                    % DC side Level 2 load shedding
    
    ub7 = P_PV1;
    ub8 = P_PV2;
end

if ILC_status ==0 % DC grid is islanded
    ub5 = P_CL2;                    % DC side Level 1 load shedding
    ub6 = P_NL2;
    ub8 = P_PV2;
end

G_lb= [lb1;lb2;lb3;lb4;lb5;lb6;lb7;lb8];
G_ub= [ub1;ub2;ub3;ub4;ub5;ub6;ub7;ub8];
end