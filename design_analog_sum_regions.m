%% Sommatore generalizzato
clear
clc
% Prima regione r1

%%% (-) input
g1_m_r1 = 1/(20*10e-3) * 8.2198e-02; % i_L
%g1_m_r1 = 8.2198e-02; % i_L
g2_m_r1 = 2.8873e+00; % v_c \approx v_o
g3_m_r1 = 4.4956e-03; % v_in

%%% (+) input
g1_p_r1 = 8.2461e-01; % i_o
g2_p_r1 = 1.4957e+01/1.8; % Vbatt = 1.8V. Dal file .txt trovo il guadagno come 1.71936/1.8

%%% Scelta della topologia
An_r1 = g1_m_r1 + g2_m_r1 + g3_m_r1; 
Ap_r1 = g1_p_r1 + g2_p_r1;

An_r1 > Ap_r1 - 1  

%%%

Gf_r1 = 1/1e3;
Rf_r1 = 1/Gf_r1

%%% Negative gains

G1_m_r1 = Gf_r1 * g1_m_r1;
G2_m_r1 = Gf_r1 * g2_m_r1;
G3_m_r1 = Gf_r1 * g3_m_r1;

R1_m_r1 = 1/G1_m_r1
R2_m_r1 = 1/G2_m_r1
R3_m_r1 = 1/G3_m_r1

%%% Positive gains

G1_p_r1 = Gf_r1 * g1_p_r1;
G2_p_r1 = Gf_r1 * g2_p_r1;

R1_p_r1 = 1 / G1_p_r1
R2_p_r1 = 1 / G2_p_r1

%%% shunt resistance on the (+)

G_gnd_n_r1 = G1_p_r1 + G2_p_r1 - Gf_r1 - G1_m_r1 - G2_m_r1 - G3_m_r1;
R_gnd_n_r1 = 1/ G_gnd_n_r1

%%
clear
clc
% Seconda regione r2

%%% (-) input
g1_m_r2 = 1/(20*10e-3) *  6.8144e-02; % i_L
%g1_m_r2 = 6.8144e-02; % i_L
g2_m_r2 = 4.0212e+00; % v_c \approx v_o
g3_m_r2 = 2.0623e-03; % v_in

%%% (+) input
g1_p_r2 = 6.8252e-01; % i_o
g2_p_r2 = 2.0376e+01/1.8; % Vbatt = 1.8V. Dal file .txt trovo il guadagno come 1.71936/1.8

%%% Scelta della topologia
An_r2 = g1_m_r2 + g2_m_r2 + g3_m_r2; 
Ap_r2 = g1_p_r2 + g2_p_r2;

An_r2 > Ap_r2-1  

%%%

Gf_r2 = 1/1e3;
Rf_r2 = 1/Gf_r2

%%% Negative gains

G1_m_r2 = Gf_r2 * g1_m_r2;
G2_m_r2 = Gf_r2 * g2_m_r2;
G3_m_r2 = Gf_r2 * g3_m_r2;

R1_m_r2 = 1/G1_m_r2
R2_m_r2 = 1/G2_m_r2
R3_m_r2 = 1/G3_m_r2

%%% Positive gains

G1_p_r2 = Gf_r2 * g1_p_r2;
G2_p_r2 = Gf_r2 * g2_p_r2;

R1_p_r2 = 1 / G1_p_r2
R2_p_r2 = 1 / G2_p_r2

%%% shunt resistance on the (+)

G_gnd_n_r2 = G1_p_r2 + G2_p_r2 - Gf_r2 - G1_m_r2 - G2_m_r2 - G3_m_r2;
R_gnd_n_r2 = 1/ G_gnd_n_r2
