close all; clear; clc;

%% Loading acquired data
training_data_y = cell(1,2);
training_data_u = cell(1,2);

%training set 2: square sweep
load('../DataExp/prbs 1 alpha.mat');   % loading alpha's
load('../DataExp/prbs 1 theta.mat');   % loading theta's
load('../DataExp/prbs 1 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 20000;
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);   

training_data_y{2} = ymeas;
training_data_u{2} = uin;

%training set 1: sweep
load('../DataExp/doublesweep 8 epsilon00242 alpha.mat');   % loading alpha's
load('../DataExp/doublesweep 8 epsilon00242 theta.mat');   % loading theta's
load('../DataExp/doublesweep 8 epsilon00242 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 20000;
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);     

training_data_y{1} = ymeas;
training_data_u{1} = uin;

dt = 0.01;
t = 0:dt:(data_end - data_begin)*dt;

%% Singular values
 s = 100;
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
close all

% Testing different sizes
n=4;
A0 = ones(n,n);
B0 = ones(n,1);
C0 = ones(2,n);
D0 = zeros(2,1);

% check P.E.O.
%checkpoe(n, training_data_u{1},t,s)

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;

training_data = iddata(training_data_y{2},training_data_u{2},dt);
opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'estimate';
%opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))]; %CHECK IF THIS WORKS BETTER
opt.N4Weight = 'MOESP';
opt.N4Horizon = [s/2 s/2 s/2]; % s = 60, symmetric Yo Yf
opt.EnforceStability = true;
opt.Advanced.DDC = 'on';
sys_init2 = n4sid(training_data, init_sys,opt);
sys_init2.Ts = 0;

%opt.Regularization.Lambda = 1e-6;
%opt.Regularization.Nominal = 'zero'; % prefer low entries in matrices for controller implementation and num. stability

sys = pem(training_data, sys_init2,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
Cbar = sys.C
D = sys.D
sys.x0 = zeros(n,1);

Config = RespConfig(Amplitude=0.01,Delay=2);
figure()
impulse(sys, Config)

eig(sys.A)

figure()
ropt = residOptions;
ropt.MaxLag = 100;
resid(training_data, sys,ropt)
%% Enforcing structure via state coordinate transformation (T = [C;CA])

T = [sys.C;sys.C*sys.A];
Ds = sys.D;
Cs = sys.C/T;
Bs = T*sys.B;
As = T*sys.A/T;

semi_struct_sys = ss(As,Bs,Cs,Ds);

%% Conditioning B to be [0;0;*;*]
close all

corrected_B = semi_struct_sys.B;
corrected_B(1) = 0;
corrected_B(2) = 0;

Acorr = semi_struct_sys.A;

Acorr(4,1) = 0;
%Acorr(3,1) = 0;

sys_init3 = idss(Acorr,corrected_B,semi_struct_sys.C,semi_struct_sys.D);
sys_init3.Ts = 0;

sys_init3.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [1,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
sys_init3.Structure.B.Free = [0;0;1;1];
sys_init3.Structure.C.Free = zeros(2,4);
sys_init3.Structure.D.Free = [0;0];

%opt.Regularization.R = 1e16*eye(10);
%opt.Regularization.Nominal = 'zero';
struct_sys = pem(training_data, sys_init3,opt);

figure()
impulse(struct_sys)

figure()
resid(training_data, struct_sys,ropt)

%% Validation (for every test do two runs)

load('../Data/doublesweep 8 epsilon00242 alpha.mat');   % loading alpha's
load('../Data/doublesweep 8 epsilon00242 theta.mat');   % loading theta's
load('../Data/doublesweep 8 epsilon00242 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 20000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);     

validation_data = [ymeas,uin];

figure()
compare(validation_data,sys, struct_sys)

figure()
resid(validation_data, struct_sys,ropt)
%% Turning the pendulum upside down

As = struct_sys.A;
As(3,4) = -As(3,4);
As(4,2) = -As(4,2);
As(4,3) = -As(4,3);

struct_sys.A = As;

Bs = struct_sys.B;
Bs(4,1) = -Bs(4,1);

struct_sys.B = Bs;

%% Discretize for implementation
str_discr_sys = c2d(struct_sys,0.01,'zoh');

%% MPC synthesis
h=dt;

mpc_A = str_discr_sys.A(1:4,1:4);
mpc_B = str_discr_sys.B(1:4);
mpc_C = eye(4);
mpc_D = zeros(4,1);

mpc_model = ss(mpc_A,mpc_B,mpc_C,mpc_D,dt);


% creating mpcobject
MPC1 = mpc(mpc_model,dt);

% setting Q,R and N
MPC1.PredictionHorizon = 40;
MPC1.Weights.OutputVariables = [1e1, 1e-1, 1e-3,1e-3];
MPC1.Weights.ManipulatedVariables = 1e-1;

%omitting the build in state estimation
setEstimator(MPC1, "custom")
MPC1.Model.Plant = setmpcsignals(MPC1.Model.Plant, 'MO', []);
% Setting up terminal cost P
%P = ....
%Y.Weight = P;