close all; clear; clc;

%% Loading acquired data
training_data_y = cell(1,2);
training_data_u = cell(1,2);

%training set 2: square sweep
load('../Data/Squaresweep 1 alpha.mat');   % loading alpha's
load('../Data/Squaresweep 1 theta.mat');   % loading theta's
load('../Data/Squaresweep 1 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 8000;
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end) - mean(theta(data_begin:data_end))];
uin = u(data_begin:data_end,2);   
uin = uin - mean(uin(data_begin:data_end));

training_data_y{2} = ymeas;
training_data_u{2} = uin;

%training set 1: sweep
load('../Data/Sweep 6 alpha.mat');   % loading alpha's
load('../Data/Sweep 6 theta.mat');   % loading theta's
load('../Data/Sweep 6 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 8000;
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end) - mean(theta(data_begin:data_end))];
uin = u(data_begin:data_end,2);     
uin = uin - mean(uin(data_begin:data_end));

training_data_y{1} = ymeas;
training_data_u{1} = uin;

dt = 0.01;
t = 0:dt:(data_end - data_begin)*dt;

%% Singular values
 s = 60;
 y_hank = hankel(ymeas(1:s),ymeas(s:end));
 [U,S,V] = svd(y_hank, 'econ');
 sing_vals = diag(S);
 figure
 semilogy(sing_vals,'o')
 grid on
 title('Singular values of the output signal')
 xlabel('Value index')
 ylabel('Singular value')

%% Run Algorithm

% Testing different sizes
n=4;
A0 = ones(n,n);
B0 = ones(n,1);
C0 = ones(2,n);
D0 = zeros(2,1);

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;

training_data = iddata(training_data_y,training_data_u,dt);
opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'zero';
%opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))]; %CHECK IF THIS WORKS BETTER
opt.N4Weight = 'MOESP';
opt.N4Horizon = [30 30 30]; % s = 60, symmetric Yo Yf
opt.EnforceStability = true;
opt.Advanced.DDC = 'on';
sys_init2 = n4sid(training_data, init_sys,opt);
sys_init2.Ts = 0;

opt.Regularization.Lambda = 1e-12;
opt.Regularization.Nominal = 'zero'; % prefer low entries in matrices for controller implementation and num. stability

sys = pem(training_data, sys_init2,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
Cbar = sys.C
D = sys.D
sys.x0 = zeros(n,1);

Config = RespConfig(Amplitude=0.01,Delay=2);
figure
impulse(sys, Config)
eig(sys.A)

%% Enforcing structure via state coordinate transformation (T = [C;CA])

T = [sys.C;sys.C*sys.A];
Ds = sys.D;
Cs = sys.C/T;
Bs = T*sys.B;
As = T*sys.A/T;

semi_struct_sys = ss(As,Bs,Cs,Ds);

%% Conditioning B to be [0;0;*;*]
corrected_B = semi_struct_sys.B;
corrected_B(1) = 0;
corrected_B(2) = 0;

sys_init3 = idss(semi_struct_sys.A,corrected_B,semi_struct_sys.C,semi_struct_sys.D);
sys_init3.Ts = 0;

sys_init3.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [1,1,1,1]; %A31 is a torsional factor induced by cable
                             [1,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
sys_init3.Structure.B.Free = [0;0;1;1];
sys_init3.Structure.C.Free = zeros(2,4);
sys_init3.Structure.D.Free = [0;0];


struct_sys = pem(training_data, sys_init3,opt);

impulse(struct_sys)
%% Validation (for every test do two runs)

load('../Data/Sweep 5 alpha.mat');   % loading alpha's
load('../Data/Sweep 5 theta.mat');   % loading theta's
load('../Data/Sweep 5 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 8000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end) - mean(theta(data_begin:data_end))];
uin = u(data_begin:data_end,2);     

validation_data = [ymeas,uin];
figure
compare(validation_data,sys, struct_sys)

%% Discretize for implementation
str_discr_sys = c2d(struct_sys,0.01,'zoh');
compare(validation_data,str_discr_sys,struct_sys);

%% Discrete Iterative Kalman Filter Design
close all

u_training = training_data_u{1};
y_training = training_data_y{1};
x_0 = [0;0;0;0];

y_est = lsim(str_discr_sys, u_training, t,x_0);

% estimate the variance
R1 = 6 * 1e-6 * eye(4); % this is the model uncertainty
R2 = 1e-10 * eye(2); % this is the measurement uncertainty
R12 = zeros(4,2);

aug_sys = ss(str_discr_sys.A,[str_discr_sys.B eye(n)],str_discr_sys.C,[str_discr_sys.D zeros(2,4)]);
[Kest,L,P] = kalman(aug_sys, R1, R2)

% redefine system

K_sys = Kest;

% test results
yest_k = lsim(K_sys, [u_training,y_training], t,x_0);
epsilon = y_training - yest_k(:,1:2);

var =  epsilon.' * epsilon

figure()
title('alpha')
plot(y_training(:,1), 'r:')
hold on
plot(y_est(:,1), 'b:')
plot(yest_k(:,1))
hold off
legend('y_training','yest','yest_k')

%% (DISCRETE VERSION) - Manual LQGR design
close all

Q = diag([1, 1, 1, 1]);
R = [0.1*1e4];
[P, cl_eig, K] = dare(str_discr_sys.A, str_discr_sys.B, Q, R)

lqsys = str_discr_sys % printing original matrices
lqsys.A = str_discr_sys.A - str_discr_sys.B*K;
DC_gain = dcgain(lqsys);
G = 1/DC_gain(1);
lqsys.B = G*lqsys.B;

figure()
pzplot(lqsys)


% adding kalman filter for simulation(aka augmenting state space)
aug_inputs = [u_training, y_training];
A_aug = lqsys.A; %  adding the kalman filter to the lqr (this makes it a lqgr)
B_aug = [lqsys.B L];
C_aug = [lqsys.C;eye(n)];
D_aug = [lqsys.D, zeros(2,2); zeros(4,3)];

lqgr_sys = ss(A_aug,B_aug,C_aug,D_aug,dt);

y_aug = step(lqgr_sys);
y_lqr = step(lqsys);

figure()
plot(y_aug(:,2,1))
title('step response theta lqr')
hold on 
plot(y_lqr(:,2))
legend('with kalman','without kalman')
grid on

stepinfo(lqsys).TransientTime
[~,~,x] = step(lqsys);

figure()
plot(K*x.');
grid on
u_max = max(abs(K*x.'));

%% (DISCRETE VERSION) - Manual LQR design
close all

Q = diag([1, 1, 1, 1]);
R = [1e3];
[P, cl_eig, K] = dare(str_discr_sys.A, str_discr_sys.B, Q, R)

lqsys = str_discr_sys % printing original matrices
lqsys.A = str_discr_sys.A - str_discr_sys.B*K;
DC_gain = dcgain(lqsys);
G = 1/DC_gain(1);
lqsys.B = G*lqsys.B;

figure()
pzplot(lqsys)

figure()
step(lqsys)
grid on

stepinfo(lqsys).TransientTime
[~,~,x] = step(lqsys);

figure()
plot(K*x.');
grid on
u_max = max(abs(K*x.'));

%% (Continuous VERSION) - Manual LQR design
close all

Q = diag([1, 1, 1, 1]);
R = [1e3];
[P, cl_eig, K] = care(struct_sys.A, struct_sys.B, Q, R)

clqsys = struct_sys % printing original matrices
clqsys.A = struct_sys.A - struct_sys.B*K;
DC_gain = dcgain(clqsys);
G = 1/DC_gain(1);
clqsys.B = G*clqsys.B;

figure()
pzplot(clqsys)

figure()
step(clqsys)
grid on

stepinfo(clqsys).TransientTime
[~,~,x] = step(clqsys);

figure()
plot(K*x.');
grid on
u_max = max(abs(K*x.'));

%% (Discrete VERSION) - Optimal LQR syntesis using genetic algorithm
close all

function [u_max, overshoot,settling_time, rise_time, f_fp] = results_d(system,K,dt)
% K is the resulting stabilizing feedbackgain from ARE
% system is the uncontrolled system of interest

system.A = system.A - system.B*K; % A-> Ak

DC_gain = dcgain(system); 
system.B = system.B / DC_gain(1); % xdot = Ak x + BG r 

info = stepinfo(system);
settling_time = [info.TransientTime].';
overshoot = [info.Overshoot].' + [info.Undershoot].';
rise_time = [info.RiseTime].';
[~,~,x] = step(system);
u_max = max(abs(1/DC_gain(1) - K*x.')); % step (1*G=1/DC_gain) + control inputs (Kx)
f_fp = max(abs(log(pole(system))/dt));
end

function f = evaluate_d(W, system,V,n,dt)
% W is a 1x3 vector of weights
% V = [q1 q2 ... r]
Q = diag(V(1:n));
R = 1;
[~, ~, K] = dare(system.A, system.B, Q, R);

[u_max, overshoot,settling_time, rise_time, f_fp] = results_d(system,K,dt);
%penalty = 1e6 / (1 + exp(-100*(u_max - 0.8)));
penalty = 1e6*max(0, u_max - 0.8)^4 + 1e6*max(0, f_fp - 18)^4; % restrict the input to 0.8 and the fastest pole to 18 rad/s

f = W(1)*u_max+[W(2), W(3)]*overshoot+[W(4) W(5)]*settling_time + W(6)*f_fp + penalty;
end

% LQR optimization
W = [1,0,1,0,0,1]; % {input | alpha overhoot | theta overshoot | alpha ts | theta ts | fastest pole coeficient | DC gain}

lb = zeros(1,n);
ub = [];     
V = optimvar('V',1,n,'LowerBound',0);
evaluateFcn = @(V) evaluate_d(W,str_discr_sys,V,n,dt);
est_max_V = 1e2;
initrange = [0,0,0,0;est_max_V,est_max_V,est_max_V,est_max_V];
options = optimoptions('ga','FunctionTolerance',1e-8, "UseParallel", ...
    true, "PopulationSize", 70, "EliteCount", 10, 'InitialPopulationRange',initrange,'PlotFcn', @gaplotbestf, ...
    MaxGenerations = 500, Display = "iter");
[sol,val] = ga(evaluateFcn, n, [], [], [], [], lb, ub, [], options);

Q = diag(sol(1:n));
R = 1;
[~,~,Kstar] = dare(str_discr_sys.A,str_discr_sys.B,Q,R)

str_discr_sys
olqsys = str_discr_sys;
olqsys.A = str_discr_sys.A - str_discr_sys.B*Kstar;
DC_gain = dcgain(olqsys);
olqsys.B = olqsys.B / DC_gain(1);

figure()
step(olqsys)
grid on

figure()
pzplot(olqsys)
grid on

DC_gain(1)

%% (DISCRETE VERSION) - pole placement design
%poles = [0.4,0.45,0.9898 + 0.0037i,0.9898 - 0.0037i,0.5]
poles = [0.96,0.97,0.98,0.9];
K = place(str_discr_sys.A, str_discr_sys.B, poles)

sfsys = str_discr_sys % printing original matrices
sfsys.A = sfsys.A - sfsys.B*K;
DC_gain = dcgain(sfsys);
G = 1/DC_gain(1);
sfsys.B = G*sfsys.B;
%pzplot(lqsys)
Config = RespConfig(Amplitude=1,Delay=0);
step(sfsys, Config)
grid on
stepinfo(sfsys).TransientTime
stepinfo(sfsys).Overshoot

%% compare results
figure()
title('Comparison of tuning methods')
step(olqsys)
hold on
step(lqsys)
step(sfsys)
hold off
legend('Optimal LQ','Manual LQ', 'Manual pole-placement')
grid on

% %% Transform system to CT for parameters
% 
% sys_c = d2c(sys,'zoh');
% impulse(sys_c)
% hold on
% impulse(sys)
% hold off
% legend()

%% Debugging
t_test = 0:0.01:120;
utest = 0.005 * chirp(t_test,0.1,120,5);
lsim(sys,utest,t_test)
%plot(utest)

%% T = [C;CA] from observability matrix!!

T = [sys.C;sys.C*sys.A];
Ds = sys.D;
Cs = sys.C/T
Bs = T*sys.B
As = T*sys.A/T

struct_sys = ss(As,Bs,Cs,Ds);

t_test = 0:0.01:120;
utest = 0.005 * chirp(t_test,0.1,120,5);
lsim(sys,utest,t_test)
legend('structured system response')