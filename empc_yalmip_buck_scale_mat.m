clc
clear variables
close all
yalmip('clear');

%% System (Buck CT)

nx = 2;
nu = 1;
ny = 1;

V_out = 5;
f_sw = 500e3; T_sw = 1/f_sw;
L_range = 8.2e-6 * [0.8 , 1.2]; L_nominal = mean(L_range);
R_Lf = 7e-3;
C_range = (10e-6 - 0.041*10e-6) * (9 + 17) * [0.9 , 1.1]; 
C_nominal = mean(C_range);

Rc_nominal = 50e-3; % 5e-3; ceramic or elettrolitic
Rc_min = 0.5*Rc_nominal;
Rc_max = 1.5*Rc_nominal;

I_max = 15;
V_in_max = 75; V_in_nominal = 50;
RL_min = V_out/I_max;
RL_max = 2*L_range(1)*f_sw/(1 - V_out/V_in_max);
RL_range = [RL_min, RL_max]; 
RL_nominal = mean(RL_range); 

V_in = V_in_nominal;
D = V_out/V_in;

RLparRc = 1/(1/RL_nominal + 1/Rc_nominal);

A = [-(RLparRc + R_Lf)/L_nominal, -RL_nominal/(L_nominal*(Rc_nominal+RL_nominal)); 
    RL_nominal/(C_nominal*(Rc_nominal+RL_nominal)), -1/(C_nominal*(Rc_nominal+RL_nominal))]; 
B1 = [RLparRc/L_nominal; -RL_nominal/(C_nominal*(Rc_nominal+RL_nominal))];
B2 = [1/L_nominal; 0];
C = [RLparRc, RL_nominal/(Rc_nominal+RL_nominal)];
D1 = -RLparRc;

data.T_sw = T_sw;

data.A = A;
data.B1 = B1;
data.B2 = B2;
data.C = C;
data.D1 = D1;

%% System discretization

Ad = expm(A*T_sw);
Bd1 = (expm(A*T_sw)-eye(nx))*inv(A)*B1;
Cd = C;

% Linearization
Bd = expm(A*(1-D)*T_sw)*T_sw*B2*V_in;
Bd2 = expm(A*T_sw)*(eye(nx)-expm(-A*D*T_sw))*inv(A)*B2;
bd = expm(A*T_sw)*(eye(nx)-expm(-A*D*T_sw))*inv(A)*B2*V_in;

%% MPC data
% Horizons
Np = 5;
Nc = 2;

% Weights
Qy = 1e2;
R = 1e-2;
Rd = 1;

% Input bounds
d_lb = 0;
d_ub = 1;

% Params bounds
% par = [x_in, io, v_in], x_in = [i_L, v_C]
i_L_b = [0, 80];
v_C_b = [0, 20];
io_b = [-5, 20];
v_in_b = V_in + [-35, 35];

% Params scaling
V_meas = 2;

ks_i_L = 5; %max(abs(i_L_b))/V_meas;
ks_v_C = 1; %max(abs(v_C_b))/V_meas;
ks_io = 10; %max(abs(io_b))/V_meas;
ks_v_in = 1;%max(abs(v_in_b))/V_meas;

% Input scaling
ks_d = 1; %10;

data.ks_d = ks_d;

%% MPC formulation

% Yalmip variables

% States
x = sdpvar(nx*ones(1,Np+1),ones(1,Np+1));
% x_in = sdpvar(nx,1);

% Inputs
d = sdpvar(nu*ones(1,Nc),ones(1,Nc)); % duty cycle

% Outputs
y = sdpvar(ny*ones(1,Np),ones(1,Np));

% Reference
% r = sdpvar(ny,1);
r = 5;

% Additional variables
% io = sdpvar(1,1); % Output current disturbance
% v_in = sdpvar(1,1); % Input voltage disturbance

% Params
par = sdpvar(4,1); % par = [x_in, io, v_in] (downscaled)

x_in = [par(1)*ks_i_L; par(2)*ks_v_C];
io = par(3)*ks_io;
v_in = par(4)*ks_v_in;

n_par = length(par);

% ===== Cost function =====
cost = 0;

% Input parametrization
ind_u = zeros(Np,1);
i = 0;
for k = 1:1:Np
	if k <= Nc
		i=i+1;
	end
	ind_u(k) = i;
end

% Output (tracking), input
for k = 1:1:Np
	cost = cost + (y{k}-r)'*Qy*(y{k}-r) + d{ind_u(k)}'*R*d{ind_u(k)};
end

% Input rate
for k = 2:1:Np
	cost = cost + (d{ind_u(k)}-d{ind_u(k-1)})'*Rd*(d{ind_u(k)}-d{ind_u(k-1)});
end

% ===== Constraints =====
constr = [];

% Initial condition
constr = [constr; x{1} == x_in];

% Prediction model
for k = 1:1:Np
	constr = [constr;
		x{k+1} == Ad*x{k} + Bd*(d{ind_u(k)}-D) + Bd1*io + Bd2*(v_in-V_in) + bd;
		y{k} == C*x{k} + D1*io];
end

% Input bounds
for i = 1:1:Nc
	constr = [constr; d_lb <= d{i} <= d_ub];
end

% Params bounds
constr = [constr; i_L_b(1) <= par(1)*ks_i_L <= i_L_b(2)]; % i_L (x1)
constr = [constr; v_C_b(1) <= par(2)*ks_v_C <= v_C_b(2)]; % v_C (x2)
constr = [constr; io_b(1) <= par(3)*ks_io <= io_b(2)]; % io
constr = [constr; v_in_b(1) <= par(4)*ks_v_in <= v_in_b(2)]; % v_in

% ===== Optimizer object (MPC controller) =====
params_in = par;

sol_out = d{1};

options = sdpsettings('verbose',1,'solver','quadprog');

mpc = optimizer(constr,cost,options,params_in,sol_out);

% Solve MP-QP for explicit MPC controller

[sol_mp,diagn,Z,value_fun,optimizer] = solvemp(constr,cost,[],par,d{1});
fprintf('\n====================\n')

% Merge polytopic regions

pu = mpt_mpsol2pu(sol_mp);
pu = pu.merge('primal');
fprintf('\n====================\n')
%pu = pu.merge('primal','optimal',true);

fprintf('\n====================\n')

% Linear separation of two sets of polyhedra (saturated regions)

ks_sigma = 1; % Separation function scaling

% Find saturated regions

pu_sat = pu.findSaturated('primal','min',d_lb,'max',d_ub);
pu_sat_l = pu_sat.Imin;
pu_sat_u = pu_sat.Imax;

% Approach via vertices
a = sdpvar(1,n_par);
b = sdpvar(1,1);
sep_fun_sdp = @(x) a*x+b;
e = -1e-03;

sep_constr = [];

for i = 1:1:length(pu_sat_l)
	V = pu.Set(pu_sat_l(i)).V;
	sep_constr = [sep_constr; sep_fun_sdp(V') <= e];
end
for i = 1:1:length(pu_sat_u)
	V = pu.Set(pu_sat_u(i)).V;
	sep_constr = [sep_constr; -sep_fun_sdp(V') <= e];
end

sep_cost = norm(a,2) + norm(b,2);

optimize(sep_constr,sep_cost);

a_opt = value(a)*ks_sigma;
b_opt = value(b)*ks_sigma;

%sep_fun = @(x) a_opt*x + b_opt; % Linear separation function

data.a_opt = a_opt;
data.b_opt = b_opt;

fprintf('\n====================\n')

% Remove saturated regions

pu_del = sort([pu_sat_l, pu_sat_u]);
fprintf('Remove saturated regions: %d region(s) removed\n',length(pu_del))
pu.remove(pu_del)
fprintf('\n====================\n')

% Generate EMPC controller

n_reg = length(pu.Set); % Number of non-saturated regions

% Extract H-representation of each region
H_c = pu.Set.forEach(@(x) x.H, 'UniformOutput', false);

% Rescale H rep. based on affine term
H_rep_aff = 5; % Desired value of the affine term in H rep.

for i=1:1:n_reg
	for j=1:1:size(H_c{i},1)
		if H_c{i}(j,end) ~= 0
			H_c{i}(j,:) = H_c{i}(j,:)*abs(H_rep_aff/H_c{i}(j,end));
		end
	end
end

% Extract affine function of each region (and downscale)
F_c = pu.Set.forEach(@(x) x.Functions('primal').F / ks_d, 'UniformOutput', false);
g_c = pu.Set.forEach(@(x) x.Functions('primal').g / ks_d, 'UniformOutput', false);

g = cell2mat(g_c)';

H = []; F = []; ni = [];
for i = 1:1:n_reg
	H = [H; H_c{i}];
	F = [F; F_c{i}];
	if i == 1
		ni = [0; size(H_c{i},1)];
	else
		ni = [ni; ni(end)+size(H_c{i},1)];
	end
end

data.n_reg = n_reg;
data.H = H;
data.F = F;
data.g = g;
data.ni = ni;

% Output regions and control function as text
file_reg = 1; %fopen("result_regions.txt", 'w');
for i = 1:1:n_reg
	fprintf(file_reg,'Region %d:\n',i);
	fprintf(file_reg,'d = %.4e i_L + %.4e v_C + %.4e i_o + %.4e v_in + %.4e\n',...
		F_c{i}(1),F_c{i}(2),F_c{i}(3),F_c{i}(4), g_c{i});
	fprintf(file_reg,'if\n');
	for j = 1:1:size(H_c{i},1)
		fprintf(file_reg,'%.4e i_L + %.4e v_C + %.4e i_o + %.4e v_in + %.4e <= 0\n',...
		H_c{i}(j,1),H_c{i}(j,2),H_c{i}(j,3),H_c{i}(j,4),-H_c{i}(j,5));
	end
	fprintf(file_reg,'\n');
end
fprintf(file_reg,'Separation function for saturated regions:\n');
fprintf(file_reg,'sigma = %.4e i_L + %.4e v_C + %.4e i_o + %.4e v_in + %.4e\n',...
	a_opt(1),a_opt(2),a_opt(3),a_opt(4),b_opt);
fprintf(file_reg,'If sigma > 0, d = %.4e; if sigma <= 0, d = %.4e\n', d_ub/ks_d, d_lb/ks_d);
fprintf(file_reg,'\n====================\n');

% save empc_params.mat a_opt b_opt F_c g_c H_c -mat

%% Estimator design
s= tf('s');
G1 = (C_nominal*(RL_nominal+Rc_nominal)*s+1)/(RL_nominal*(1+s*C_nominal*Rc_nominal));
Estim = [-G1, 1];

%% Simulation

N = 1600;

NUM_SIM = 1;

x0 = [0, 0]';

v_in = 10;
% Ceramic
% N_v_in_start = 1400;
% N_v_in_end = 1500;
% N_io_start = 1200;
% N_io_end = 1300;

% Elettrol
N_v_in_start = 1200;
N_v_in_end = 1300;
N_io_start = 1000;
N_io_end = 1100;

dt_sim = T_sw/100;

for idx=1:NUM_SIM
    RL = RL_min + 0.5*(RL_max-RL_min);
    io = (I_max - 5/RL);
    Rc = Rc_min +  0.5*(Rc_max-Rc_min);
    Co = C_range(1) + 0.5*(C_range(2)-C_range(1));
    L = L_range(1) + 0.5*(L_range(2)-L_range(1));

    if idx <= 250
        values_a{idx} = [RL, Rc, Co, L];
    else 
        values_b{idx-250} = [RL, Rc, Co, L];
    end
    fprintf(1, "Simulation %d\n", idx);
    %fprintf(1, "i_o step = %f\n", io);

    RLparRc = 1/(1/RL + 1/Rc);
    data.A = [-(RLparRc + R_Lf)/L, -RL/(L*(Rc+RL)); 
        RL/(Co*(Rc+RL)), -1/(Co*(Rc+RL))]; 
    data.B1 = [RLparRc/L; -RL/(Co*(Rc+RL))];
    data.B2 = [1/L; 0];
    data.C = [RLparRc, RL/(Rc+RL)];
    data.D1 = -RLparRc;
    data.T_sw = T_sw;
    if idx <= 250
        sim_out_a{idx} = sim('empc_yalmip_buck_scale_sl');
    else
        sim_out_b{idx-250} = sim('empc_yalmip_buck_scale_sl');
    end
end

%% Plot
% for idx = 1:length(sim_out)
%     figure(1), hold on;  plot(sim_out{idx}.d.time, sim_out{idx}.d.data);
%     grid on;
%     figure(2), hold on; plot(sim_out{idx}.y.time, sim_out{idx}.y.data);
%     grid on;
%     % figure(3), hold on; plot(sim_out{idx}.vc.time, sim_out{idx}.vc.data);
%     % plot(sim_out{idx}.vc_hat.time, sim_out{idx}.vc_hat.data);
%     % grid on;
% end

% cd G:\
% save sim_ceramic1.mat values_a sim_out_a -mat  
% save sim_ceramic2.mat values_b sim_out_b -mat  
% save sim_elettrol1.mat values_a sim_out_a -mat  
% save sim_elettrol2.mat values_b sim_out_b -mat  
save sim_elettrol_nom.mat sim_out_a -mat  
% simoutdata = whos('sim_out_a');
% fprintf(1,"Saved two files of size %.2f GB\n", simoutdata.bytes/(1024^3));