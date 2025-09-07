clear; close all; clc;

load("empc_params.mat");

%% MAIN LOGIC
% for each region
% for each inequality
%   normalize to have max gain e.g., = .7
%   determine number of required resistors on + and -
%   compute conductances

input_str = {"0.2*i_L","v_C","i_o*0.1","V_in","V_batt=5V"};
output_file = fopen("analog_design_empc_results_scaled_iL.txt", "w");

for r_idx = 1:3
    if r_idx <= 2
        fprintf(output_file,"############## Region %d ############## \n", r_idx);
        % LOAD AND GET GAINS
        Reg_conds = H_c{r_idx};
        % A*p <= b  --> A*p - b <= 0
        Reg_conds(:,end) = - Reg_conds(:,end);  
    else
        fprintf(output_file,"############## Sigma ############## \n");
        Reg_conds = [a_opt, b_opt];
    end
    % last gain uses V_batt as input, so cmpute g5 s.t. b = g_5*V_batt
    V_batt = 5; Reg_conds(:,end) = Reg_conds(:,end)./V_batt;
    for ineq = 1:size(Reg_conds,1)
        % gains = Reg_conds(ineq, :)./[20*10e-3, 1, 1, 1, 1];
        gains = Reg_conds(ineq, :)./[1, 1, 1, 1, 1];
        gains = 0.7*(gains./max(abs(gains)));
        for sign = [1, -1]
            if sign == 1
                sign_str = "+";
                idx_s = find(gains>0);
                gs = gains(idx_s)';
            else
                sign_str = "-";
                idx_s = find(gains<0);
                gs = -gains(idx_s)';
            end
            % fix resistor at ground
            Gs = zeros(length(gs)+1,1);
            RG = 10e3;  Gs(end) = 1/RG;
            Gs(1:end-1) = - (gs*ones(1,length(gs))-eye(length(gs)))\(Gs(end)*gs);
            Rs = 1./Gs;
            
            fprintf(output_file,"Condition %d (%s side):\n", ineq, sign_str);
            descr = ""; 
            if isscalar(Rs)
                fprintf(output_file,"ground (no inputs, compare with zero)\n");
            end
            for ii = 1:length(Rs)-1
                fprintf(output_file,input_str(idx_s(ii)) + " with R = " + num2str(Rs(ii)*1e-3) + "kOhm\n");
            end
        end
        fprintf(output_file,"\n");
    end
end            
fclose(output_file);