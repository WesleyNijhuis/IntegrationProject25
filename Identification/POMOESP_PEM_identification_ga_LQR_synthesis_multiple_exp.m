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
ymeas = [alpha(data_begin:data_end) - mean(alpha), theta(data_begin:data_end) - mean(theta)];
uin = u(data_begin:data_end,2);     

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
ymeas = [alpha(data_begin:data_end) - mean(alpha), theta(data_begin:data_end) - mean(theta)];
uin = u(data_begin:data_end,2);     

training_data_y{1} = ymeas;
training_data_u{1} = uin;

dt = 0.01;

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

training_data = iddata(training_data_y,training_data_u,dt);
opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'zero';
%opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))]; %CHECK IF THIS WORKS BETTER
opt.EnforceStability = true;
opt.Advanced.DDC = 'on';

[Ap,Bp,Cp,Dp,x0p] = pomoesp(uin,ymeas,s,n); % only does in on the last training set loaded

sys_init2 = idss(Ap,Bp,Cp,Dp);
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

data_end = 10000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha), theta(data_begin:data_end) - mean(theta)];
uin = u(data_begin:data_end,2);              
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

validation_data = [ymeas,uin];
figure
compare(validation_data,sys, struct_sys)

% %% (CONTINUOUS VERSION) - Optimal LQR syntesis using genetic algorithm
% close all
% 
% function [u_max, overshoot,settling_time, rise_time, f_fp] = results(system,K)
% % K is the resulting stabilizing feedbackgain from ARE
% % system is the uncontrolled system of interest
% 
% system.A = system.A - system.B*K; % A-> Ak
% 
% DC_gain = dcgain(system); 
% system.B = system.B / DC_gain(1); % xdot = Ak x + BG r 
% 
% info = stepinfo(system);
% settling_time = [info.TransientTime].';
% overshoot = [info.Overshoot].' + [info.Undershoot].';
% rise_time = [info.RiseTime].';
% [~,~,x] = step(system);
% u_max = max(abs(1/DC_gain(1) - K*x.')); % step (1*G=1/DC_gain) + control inputs (Kx)
% f_fp = max(abs(pole(system)));
% end
% 
% function f = evaluate(W, system,V,n)
% % W is a 1x3 vector of weights
% % V = [q1 q2 ... r]
% Q = diag(V(1:n));
% R = 1;
% [~, ~, K] = care(system.A, system.B, Q, R);
% 
% [u_max, overshoot,settling_time, rise_time,~] = results(system,K);
% %penalty = 1e6 / (1 + exp(-100*(u_max - 0.8)));
% penalty = 1e6*max(0, u_max - 0.8)^4;
% 
% f = W(1)*u_max+[W(2), W(3)]*overshoot+[W(4) W(5)]*settling_time + W(6)*f_fp + penalty;
% end
% 
% % LQR optimization
% W = [1e-1,1e1,1,1e-1,1,1]; % {input | alpha overhoot | theta overshoot | alpha ts | theta ts | fastest pole}
% lb = zeros(1,n);
% ub = [];     
% V = optimvar('V',1,n,'LowerBound',0);
% evaluateFcn = @(V) evaluate(W,struct_sys,V,n);
% options = optimoptions('ga','FunctionTolerance',1e-8, "UseParallel", ...
%     true, "PopulationSize", 100, "EliteCount", 10,'PlotFcn', @gaplotbestf, ...
%     Display = "iter");
% [sol,val] = ga(evaluateFcn, n, [], [], [], [], lb, ub, [], options);
% 
% Q = diag(sol(1:n));
% R = 1;
% [~,~,Kstar] = care(struct_sys.A,struct_sys.B,Q,R)
% 
% struct_sys
% olqsys = struct_sys;
% olqsys.A = struct_sys.A - struct_sys.B*Kstar;
% DC_gain = dcgain(olqsys);
% olqsys.B = olqsys.B / DC_gain(1);
% 
% figure()
% step(olqsys)
% grid on
% 
% figure()
% pzplot(olqsys)
% grid on

%% Turning the pendulum upside down

%% Discretize for implementation
str_discr_sys = c2d(struct_sys,0.01,'zoh')

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

%% (DISCRETE VERSION) - Manual LQR design
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

figure()
step(lqsys)
grid on

stepinfo(lqsys).TransientTime
[~,~,x] = step(lqsys);

figure()
plot(K*x.');
grid on
u_max = max(abs(K*x.'));

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

%% PO MOESP ALGORITHM
function [A,B,C,D,x0]=pomoesp(u,y,s,n)

 % Instructions:
 % Implement your subspace ID method here.
 % Use the following function inputs and outputs.
 % Function INPUT 
% u         system input (matrix of size N x m)
 % y         system output (matrix of size N x l)
 % s         block size (scalar)
 % output
 % A         System matrix A (matrix of size n x n)
 % B         System matrix B (matrix of size n x m)
 % C         System matrix C (matrix of size l x n)
 % D         System matrix D (matrix of size l x m)
 % x0        Initial state (vector of size n x one)

N_u = length(u);
N_y = length(y);
U_p = zeros(s, N_u - 2 * s + 1); 
U_f = zeros(s, N_u - 2 * s + 1); 
Y_p = zeros(s, N_y - 2 * s + 1); 
Y_f = zeros(s, N_y - 2 * s + 1); 

for i = 1:(N_u - 2 * s + 1)
    U_p(:, i) = u(i:i+s-1);
    Y_p(:, i) = y(i:i+s-1);
end

% Future input/output
U_f(:, i) = u(i+s:i+2*s-1);
Y_f(:, i) = y(i+s:i+2*s-1);
Z = [U_p; Y_p]; 
UFYZ = [U_f; Z; Y_f]; 
[~, Q] = qr(UFYZ', 0); 
L_act = Q';
row_start = size(U_f, 1) + size(Z, 1) + 1; 
col_start = size(U_f, 1) + 1;              
L_32 = L_act(row_start:end, col_start:col_start + size(Z, 1) - 1);
singular_values = svd(L_32);

figure;
semilogy(singular_values, 'o-'); 
grid on;
xlabel('Index');
ylabel('Singular Value (Log Scale)');
title('Singular Values of L_{32}');
[U, S, V] = svd(L_32,'econ');
C = U(1, 1:n); 
A = U(1:end - 1, 1:n) \ U(2:end, 1:n);
[N, m] = size(u); 
l = 1; 
Phi = []; 
Y_vec = []; 
for k = 1:N
    CAk = C * (A^(k-1));
    sum_term = zeros(l, n*m); 
    for j = 0:(k-1)
        sum_term = sum_term + kron(u(j+1, :)', C * (A^(k-1-j)));
    end
    D_term = kron(u(k, :), eye(l));
    Phi_k = [CAk, sum_term, D_term];
    Phi((k-1)*l+1:k*l, :) = Phi_k;
end
Y_vec = [Y_vec; y(k, :)];
Phi_m = zeros(N * l, n + n * m + l * m);
Phi_m(:, 1:n) = Phi(:, 1:n);
Phi_m(:, n+n*m+1:n+n*m+m) = Phi(:, n+n*m+1:n+n*m+m);
Phi_m(l+1:N*l, n+1:n+n*m) = Phi(1:l*(N-1), n+1:n+n*m);
Phi_m(l, n+1:n+n*m) = zeros(l, n*m);
Phi = Phi_m;
theta = (Phi' *Phi)\Phi' * Y_vec; 
x0 = theta(1:n); 
B = reshape(theta(n+1:end-1), [n, m]);
D = reshape(theta(end), [l, m]); 
end