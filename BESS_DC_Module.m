function [BDC_lb, BDC_ub, BDC_A, BDC_b, SOC2_init] = BESS_DC_Module (t, Num_var, Eff_BESS, P_BESS2_max, CAP2, SOC2_init, SOC2_min, SOC2_max, BESS_DC_status)
% Upper and lower bounds
lb3 = zeros(Num_var,1);                    % P_BESS2dis? ???? 0
lb4 = -P_BESS2_max*ones(Num_var,1);        % P_BESS2chg? ???? -P_BESS2_max
ub3 = P_BESS2_max*ones(Num_var,1);         % P_BESS2dis? ???? P_BESS2_max
ub4 = zeros(Num_var,1);                    % P_BESS2chg? ???? 0

if (BESS_DC_status ==0 || CAP2==0) 
    lb4=zeros(Num_var,1);
    ub3=zeros(Num_var,1);
    
    SOC2_init=0;
    SOC2_min=0;
    
    CAP2=1;
 
end

BDC_lb= [lb3;lb4];
BDC_ub= [ub3;ub4];

A3 = -[zeros(Num_var,Num_var*21), (t*100/CAP2)*tril(ones(Num_var))/Eff_BESS, (t*100/CAP2)*tril(ones(Num_var))*Eff_BESS, zeros(Num_var,Num_var*6)];
b3 = (-SOC2_init + SOC2_max)*ones(Num_var,1);       % SOC2? max?? ??? ?
A4 = [zeros(Num_var,Num_var*21), (t*100/CAP2)*tril(ones(Num_var))/Eff_BESS, (t*100/CAP2)*tril(ones(Num_var))*Eff_BESS, zeros(Num_var,Num_var*6)];
b4 = (SOC2_init - SOC2_min).*ones(Num_var,1);       % SOC2? min?? ?? ?
BDC_A = [A3; A4];
BDC_b = [b3; b4];
end