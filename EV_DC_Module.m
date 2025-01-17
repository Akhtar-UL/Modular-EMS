function [EDC_lb, EDC_ub, EDC_A, EDC_b, EV2SOC_init,EVP] = EV_DC_Module(t, Num_var, Eff_EV2, P_EV2_max, EV_CAP2, EV2SOC_init, EV2SOC_min, EV2SOC_max, Ta2, Td2, EV2_SOCT, Grid_status, ILC_status, EV_DC_status)
Ta2=Ta2/t;
Td2=Td2/t;
EVP = ones(Num_var, 1);
EVP(1:Ta2-1) = 0; % Before arrival time, EV is not present
EVP(Td2+1:end) = 0; % After departure time, EV is not present

% Upper and lower bounds
lb1 = zeros(Num_var, 1);                     % P_EVdis >= 0
lb2 = -P_EV2_max * EVP;         % P_EVchg <= -P_EV1_max
ub1 = P_EV2_max * EVP;          % P_EVdis <= P_EV1_max
ub2 = zeros(Num_var, 1);                     % P_EVchg <= 0

if (EV_DC_status == 0 || EV_CAP2==0) 
    lb2 = zeros(Num_var, 1);
    ub1 = zeros(Num_var, 1);
    
    EV2SOC_init = 0;
    EV2SOC_min = 0;
    EV2_SOCT=0;
    
    EV_CAP2=1;
end
if (Grid_status == 0 || ILC_status==0)
    EV2_SOCT=EV2SOC_init;
end

EDC_lb = [lb1; lb2];
EDC_ub = [ub1; ub2];

A1 = -[zeros(Num_var,Num_var*25), (t*100/EV_CAP2)*tril(ones(Num_var))/Eff_EV2, (t*100/EV_CAP2)*tril(ones(Num_var))*Eff_EV2, zeros(Num_var,Num_var*2)];
b1 = (-EV2SOC_init + EV2SOC_max)*ones(Num_var,1);       % SOC1? max?? ??? ?
A2 = [zeros(Num_var,Num_var*25), (t*100/EV_CAP2)*tril(ones(Num_var))/Eff_EV2, (t*100/EV_CAP2)*tril(ones(Num_var))*Eff_EV2, zeros(Num_var,Num_var*2)];
b2 = (EV2SOC_init - EV2SOC_min).*ones(Num_var,1);       % SOC1? min?? ?? ?

% Define time window matrix
time_window = tril(ones(Num_var));
time_window(1:Td2-1, :) = 0;    % Before arrival time, no charging/discharging
time_window(Td2+1:end, :) = 0;  % After departure time, no charging/discharging

% Define EV presence vector
EVT = ones(Num_var, 1);
EVT(1:Td2-1) = 0;               % Before arrival time, EV is not present
EVT(Td2+1:end) = 0;             % After departure time, EV is not present

A3 = [zeros(Num_var,Num_var*25), (t*100/EV_CAP2)*time_window/Eff_EV2, (t*100/EV_CAP2)*time_window*Eff_EV2, zeros(Num_var,Num_var*2)];
b3 = (EV2SOC_init - EV2_SOCT).*EVT;       % SOC1? min?? ?? ?



%     % Adjust constraints to ensure target SOC at departure
%     A3 = zeros(Num_var, Num_var * 27);
%     A3(departure_idx*ones(Num_var, 1), Num_var * 23 + 1: Num_var * 23 + Num_var) = (t * 100 / EV_CAP1) / Eff_EV1 * tril(ones(Num_var)) + (t * 100 / EV_CAP1) * Eff_EV1 * tril(ones(Num_var));
%     b3 = (EV1_SOCT - EV1SOC_init) * ones(Num_var, 1);

EDC_A = [A1; A2; A3];
EDC_b = [b1; b2; b3];
end
