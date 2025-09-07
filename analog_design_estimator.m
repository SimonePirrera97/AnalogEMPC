clear all
clc

%% Estimator Analog Design

s = tf('s');

% parametri del Buck
V_out = 5;
f_sw = 500e3; T_sw = 1/f_sw;
L_range = 8.2e-6 * [0.8 , 1.2]; L_nominal = mean(L_range);
R_in = 7e-3;
C_range = (10e-6 - 0.041*10e-6) * (9 + 17) * [0.9 , 1.1]; 
C_nominal = mean(C_range);
Rcnom = 5e-3; 
Rc_min = 0.5*Rcnom;
Rc_max = 1.5*Rcnom;
R_sw = 10e-3;

% 
I_max = 15;
V_in_max = 75; V_in_nominal = 50;
RL_min = V_out/I_max;
RL_max = 2*L_range(1)*f_sw/(1 - V_out/V_in_max);
RL_range = [RL_min, RL_max]; 
RL_nominal = mean(RL_range); 

% Condizioni operative e parametri nominali
RL = RL_nominal; 
L = L_nominal; Co = C_nominal;
V_in = V_in_nominal;
D = V_out/V_in;
RLparRc = 1/(1/RL + 1/Rcnom);
RLnom = RL_nominal;

% FdT dello stimatore: io_estim = - G1*vo + iL 
G1 = (Co*(RLnom+Rcnom)*s+1)/(RLnom*(1+s*Co*Rcnom));
Estim = [-G1, 1];

%%

clc
% design dello stimatore.

k_s_io = .1;

[z_eps, p_eps, k_eps] = zpkdata(G1);
z_eps = z_eps{1}
p_eps = p_eps{1}
dc_gain = k_eps * z_eps/p_eps;

zpk(G1)
G1check = zpk(dc_gain*(1-s/z_eps)/(1-s/p_eps))


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Soluzione 1 %%%%%%%%%%%%%%%%%%%%%%%%%%

% Partitore con C0 sul percorso diretto. () per la vo. R1 è verso GND.
C0 = 10000e-12
R0 = - 1 / (z_eps * C0) 
R01 = - 1 /(p_eps * C0)
R1 = 1 / (1 / R01 - 1/R0)

% Rete non invertente sul (-) del secondo opamp per la vo. 
g_vo  =  dc_gain * k_s_io / (R1 / (R0 + R1))  
Rg0 = 10e3
Rg1 = g_vo * Rg0

% Partitore sul (+) del secondo opamp per la iL. R3 è verso GND.
R3 = 10e3
R2 = R3 * ((1 + g_vo - k_s_io) / k_s_io )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Soluzione 2 %%%%%%%%%%%%%%%%%%%%%%%%%%

C0 = 1000e-12
R0 = - 1 / (z_eps * C0) 
R1 = R0 * dc_gain * k_s_io
C1 = - 1/ (R1 * p_eps)
R0_prime = 1 / (C0/(C1*R1) - 1/R0);

g_il = 20 * 10e-3;
R3 = 10e3;
R2 = R3 * (g_il * 1/k_s_io * (1 + R1/(R0_prime^-1 + R0^-1)^-1) - 1 )