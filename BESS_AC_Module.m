function [BAC_lb, BAC_ub, BAC_A, BAC_b, SOC1_init] = BESS_AC_Module (t, Num_var, Eff_BESS, P_BESS1_max, CAP1, SOC1_init, SOC1_min, SOC1_max, BESS_AC_status)
% Upper and lower bounds
lb1 = zeros(Num_var,1);                    % P_BESS1dis? ???? 0
lb2 = -P_BESS1_max*ones(Num_var,1);        % P_BESS1chg? ???? -P_BESS1_max
ub1 = P_BESS1_max*ones(Num_var,1);         % P_BESS1dis? ???? P_BESS1_max
ub2 = zeros(Num_var,1);                    % P_BESS1chg? ???? 0

if (BESS_AC_status ==0 || CAP1 ==0)
    lb2=zeros(Num_var,1);
    ub1=zeros(Num_var,1);
    
    SOC1_init=0;
    SOC1_min=0;
    
    CAP1=1;
 
end

BAC_lb= [lb1;lb2];
BAC_ub= [ub1;ub2];

A1 = -[zeros(Num_var,Num_var*19), (t*100/CAP1)*tril(ones(Num_var))/Eff_BESS, (t*100/CAP1)*tril(ones(Num_var))*Eff_BESS, zeros(Num_var,Num_var*8)];
b1 = (-SOC1_init + SOC1_max)*ones(Num_var,1);       % SOC1? max?? ??? ?
A2 = [zeros(Num_var,Num_var*19), (t*100/CAP1)*tril(ones(Num_var))/Eff_BESS, (t*100/CAP1)*tril(ones(Num_var))*Eff_BESS, zeros(Num_var,Num_var*8)];
b2 = (SOC1_init - SOC1_min).*ones(Num_var,1);       % SOC1? min?? ?? ?
BAC_A = [A1; A2];
BAC_b = [b1; b2];
end