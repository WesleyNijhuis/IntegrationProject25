close all; clear; clc;

%% Loading acquired data
load('../Data/Sweep 6 alpha.mat');   % loading alpha's
load('../Data/Sweep 6 theta.mat');   % loading theta's
load('../Data/Sweep 6 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 8000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);            
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

plot(uin)
title('uin')
figure
plot(ymeas(:,1))
title('outputs')
hold on
plot(ymeas(:,2))
hold off
legend('alpha','theta')
%% Singular values
 s = 20;
 y_hank = hankel(ymeas(1:s),ymeas(s:end));
 [U,S,V] = svd(y_hank, 'econ');
 sing_vals = diag(S);
 figure
 semilogy(sing_vals,'o')
 grid on
 title('Singular values of the output signal')
 xlabel('Value index')
 ylabel('Singular value')

%% Initial guess can be entered in the function theta2matrices!
% function params = initialguess(mp,Lr,mp,Lp,Rm,kt,km,Jm)
% J1=Jr
% p1 = (g*mp^2*(1/4)*Lp^2*Lr)/()
% 
% theta_init(0.095,0.085,)
%% Splitting into training/validation set
%TODO: in our case the training set and validation cant be split as we're
%doing a frequency sweep
% f_t = 0.07;                 % fraction of data dedicated to training
% ut = u(1:round(end*f_t));            % input for train-set
% uv = u(round(end*f_t)+1:end);          % input for validation-set
% yt = ymeas(1:round2even(end*f_t));        % output for train-set
% yv = ymeas(round2even(end*f_t)+1:end);      % output for validation-set
%% Run Algorithm

% Testing different sizes
n=4;
A0 = ones(n,n);
B0 = ones(n,1);
C0 = ones(2,n);
D0 = zeros(2,1);


% BIG TODO: CHANGE THE SS TO DISCRETIZED, DIFFERENT STRUCTURE?

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;

training_data = iddata(ymeas, uin, dt);
opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'zero';
opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))]; %CHECK IF THIS WORKS BETTER
opt.InitializeMethod = 'n4sid';
opt.N4Weight = 'MOESP';
opt.EnforceStability = false;
opt.Advanced.DDC = 'on';
sys_init2 = n4sid(training_data, init_sys,opt);

opt.Regularization.Lambda = 1e-6;
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

%% Enforcing structure via state coordinate transformation

T = [sys.C;sys.C*sys.A];
Ds = sys.D;
Cs = sys.C/T
Bs = T*sys.B
As = T*sys.A/T

struct_sys = ss(As,Bs,Cs,Ds);

%% Validation (for every test do two runs)

load('../Data/Sweep 5 alpha.mat');   % loading alpha's
load('../Data/Sweep 5 theta.mat');   % loading theta's
load('../Data/Sweep 5 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 10000; %for debugging
data_begin = 1;
ymeas = [(alpha(data_begin:data_end) - mean(alpha(data_begin:data_end))), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);              
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

validation_data = [ymeas,uin];
figure
compare(validation_data,sys, struct_sys)

%% Optimal LQR syntesis using genetic algorithm

function [u_max, overshoot,settling_time, rise_time] = results(system,K)
% K is the resulting stabilizing feedbackgain from DARE
% system is the uncontrolled system of interest
system.A = system.A - system.B*K;
info = stepinfo(system);
settling_time = [info.TransientTime].';
overshoot = [info.Overshoot].' + [info.Undershoot].';
rise_time = [info.RiseTime].';
[~,~,x] = step(system);
u_max = max(abs(K*x.'));
end

function f = evaluate(W, system,V,n)
% W is a 1x3 vector of weights
% V = [q1 q2 ... r]
Q = diag(V(1:n));
R = 1;
[~, ~, K] = dare(system.A, system.B, Q, R);

[u_max, overshoot,settling_time, rise_time] = results(system,K);
penalty = 1e6 / (1 + exp(-100*(u_max - 0.8)));


f = W(1)*u_max+[W(2), W(3)]*overshoot+[W(4) W(5)]*settling_time + penalty;

end



% LQR optimization
W = [0,1,1,10,10]; % {input | alpha overhoot | theta overshoot | alpha ts | theta ts}
lb = zeros(1,n);
ub = [];     
V = optimvar('V',1,n,'LowerBound',0);
evaluateFcn = @(V) evaluate(W,struct_sys,V,n);
options = optimoptions('ga','FunctionTolerance',1e-8, "UseParallel", ...
    true, "PopulationSize", 100, "EliteCount", 10,'PlotFcn', @gaplotbestf, ...
    Display = "iter");
[sol,val] = ga(evaluateFcn, n, [], [], [], [], lb, ub, [], options)
Q = diag(sol(1:n));
R = 1;
[~,~,Kstar] = dare(struct_sys.A,struct_sys.B,Q,R);
olqsys = struct_sys
olqsys.A = struct_sys.A - struct_sys.B*Kstar;
DC_gain = dcgain(olqsys)
olqsys.B = olqsys.B / DC_gain(1);

step(olqsys)

%% (modal) Transformation
sys_m = canon(sys, 'modal)')

step(sys, Config)
hold on
step(sys_m, Config)
hold off
legend()

%% pole placement design
poles = [0.4,0.45,0.9898 + 0.0037i,0.9898 - 0.0037i,0.5]
poles = [0.95, 0.96,0.97,0.98,0.99];
K = place(sys_m.A, sys_m.B, poles)

sfsys = sys_m % printing original matrices
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
%[yout,tout]


%% LQR design
q1 = 3e3;
q2 = 1e1;
q3 = 1e1;
Q = diag([q1, q1, q2, q2, q3]);
R = [1];
[P, cl_eig, K] = dare(sys_m.A, sys_m.B, Q, R)

lqsys = sys_m % printing original matrices
lqsys.A = lqsys.A - lqsys.B*K;
DC_gain = dcgain(lqsys);
G = 1/DC_gain(1);
lqsys.B = G*lqsys.B;
%pzplot(lqsys)
step(lqsys)
grid on
stepinfo(lqsys).TransientTime
%[yout,tout]

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