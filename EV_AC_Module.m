function [EAC_lb, EAC_ub, EAC_A, EAC_b, EV1SOC_init,EVP] = EV_AC_Module(t, Num_var, Eff_EV1, P_EV1_max, EV_CAP1, EV1SOC_init, EV1SOC_min, EV1SOC_max, Ta1, Td1, EV1_SOCT, Grid_status, EV_AC_status)
Ta1=Ta1/t;
Td1=Td1/t;
EVP = ones(Num_var, 1);
EVP(1:Ta1-1) = 0; % Before arrival time, EV is not present
EVP(Td1+1:end) = 0; % After departure time, EV is not present

% Upper and lower bounds
lb1 = zeros(Num_var, 1);                     % P_EVdis >= 0
lb2 = -P_EV1_max * EVP;         % P_EVchg <= -P_EV1_max
ub1 = P_EV1_max * EVP;          % P_EVdis <= P_EV1_max
ub2 = zeros(Num_var, 1);                     % P_EVchg <= 0

if (EV_AC_status == 0 || EV_CAP1==0)
    lb2 = zeros(Num_var, 1);
    ub1 = zeros(Num_var, 1);
    
    EV1SOC_init = 0;
    EV1SOC_min = 0;
    EV1_SOCT=0;
    
    EV_CAP1=1;
end

if Grid_status == 0
    EV1_SOCT=EV1SOC_init;
end

EAC_lb = [lb1; lb2];
EAC_ub = [ub1; ub2];

A1 = -[zeros(Num_var,Num_var*23), (t*100/EV_CAP1)*tril(ones(Num_var))/Eff_EV1, (t*100/EV_CAP1)*tril(ones(Num_var))*Eff_EV1, zeros(Num_var,Num_var*4)];
b1 = (-EV1SOC_init + EV1SOC_max)*ones(Num_var,1);       % SOC1? max?? ??? ?
A2 = [zeros(Num_var,Num_var*23), (t*100/EV_CAP1)*tril(ones(Num_var))/Eff_EV1, (t*100/EV_CAP1)*tril(ones(Num_var))*Eff_EV1, zeros(Num_var,Num_var*4)];
b2 = (EV1SOC_init - EV1SOC_min).*ones(Num_var,1);       % SOC1? min?? ?? ?

% Define time window matrix
time_window = tril(ones(Num_var));
time_window(1:Td1-1, :) = 0;    % Before arrival time, no charging/discharging
time_window(Td1+1:end, :) = 0;  % After departure time, no charging/discharging

% Define EV presence vector
EVT = ones(Num_var, 1);
EVT(1:Td1-1) = 0;               % Before arrival time, EV is not present
EVT(Td1+1:end) = 0;             % After departure time, EV is not present

A3 = [zeros(Num_var,Num_var*23), (t*100/EV_CAP1)*time_window/Eff_EV1, (t*100/EV_CAP1)*time_window*Eff_EV1, zeros(Num_var,Num_var*4)];
b3 = (EV1SOC_init - EV1_SOCT).*EVT;       % SOC1? min?? ?? ?



%     % Adjust constraints to ensure target SOC at departure
%     A3 = zeros(Num_var, Num_var * 27);
%     A3(departure_idx*ones(Num_var, 1), Num_var * 23 + 1: Num_var * 23 + Num_var) = (t * 100 / EV_CAP1) / Eff_EV1 * tril(ones(Num_var)) + (t * 100 / EV_CAP1) * Eff_EV1 * tril(ones(Num_var));
%     b3 = (EV1_SOCT - EV1SOC_init) * ones(Num_var, 1);

EAC_A = [A1; A2; A3];
EAC_b = [b1; b2; b3];
end
