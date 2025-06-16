close all; clear; clc;

%% Loading acquired data
training_data_y = cell(1,2);
training_data_u = cell(1,2);

%training set 2: square sweep
load('../Data/colored0.05 prbs 1 alpha.mat');   % loading alpha's
load('../Data/colored0.05 prbs 1 theta.mat');   % loading theta's
load('../Data/colored0.05 prbs 1 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 20000;
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
%ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2) - mean(u(data_begin:data_end,2));   

training_data_y{2} = ymeas;
training_data_u{2} = uin;

% %training set 1: sweep
% load('../Data/doublesweep 8 epsilon00242 alpha.mat');   % loading alpha's
% load('../Data/doublesweep 8 epsilon00242 theta.mat');   % loading theta's
% load('../Data/doublesweep 8 epsilon00242 input.mat');   % loading inputs
% alpha = alpha(:,2);
% theta = theta(:,2);
% 
% data_end = 20000;
% data_begin = 1;
% ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
% uin = u(data_begin:data_end,2);     
% 
% training_data_y{1} = ymeas;
% training_data_u{1} = uin;

dt = 0.01;
t = 0:dt:(data_end - data_begin)*dt;

%% plotting FFTs
close all;

figure()
plot((1/dt)/20001*(0:20001-1),abs(fft(u(:,2))))
% hold on
% plot((1/dt)/20001*(0:20001-1),abs(fft(theta)))
% xline(0.01,'k--')
% xline(12,'k--')
% hold off
grid on
legend('FFT of alpha','FFT of theta')
xlim([0,30]);
ylim([0,600]);
yscale log

%% Setting system parameters for simulink
h=0.01;
epsilon=0.00242;

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

init_sys = idss(A0, B0, C0, D0); %x0 is 0
init_sys.Ts = 0;

training_data = iddata(training_data_y{2},training_data_u{2},dt);
opt = ssestOptions('Display','off','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'estimate';
%opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))]; %CHECK IF THIS WORKS BETTER
opt.N4Weight = 'MOESP';
opt.N4Horizon = [s/2 s/2 s/2]; % This sets Moesp to PO-MOESP
opt.EnforceStability = true;
opt.Advanced.DDC = 'on';
opt.WeightingFilter = [1e-3,12];

sys_init2 = n4sid(training_data, init_sys,opt);
sys_init2.Ts = 0;

opt.Regularization.Lambda = 1e-1;
opt.Regularization.Nominal = 'zero'; % prefer low entries in matrices for controller implementation and num. stability

sys = pem(training_data, sys_init2,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
Cbar = sys.C
D = sys.D
sys.x0 = zeros(n,1);

%Config = RespConfig(Amplitude=0.01,Delay=2);
%figure()
%impulse(sys)

%eig(sys.A)

%figure()
ropt = residOptions;
ropt.MaxLag = 100;
%resid(training_data, sys,ropt)
%%% Enforcing structure via state coordinate transformation (T = [C;CA])

T = [sys.C;sys.C*sys.A];
Ds = sys.D;
Cs = sys.C/T;
Bs = T*sys.B;
As = T*sys.A/T;

semi_struct_sys = ss(As,Bs,Cs,Ds);

%%% Conditioning B to be [0;0;*;*]
%close all

corrected_B = semi_struct_sys.B;
corrected_B(1) = 0;
corrected_B(2) = 0;

Acorr = semi_struct_sys.A;

Acorr(4,1) = 0;
Acorr(3,1) = 0;

sys_init3 = idss(Acorr,corrected_B,semi_struct_sys.C,semi_struct_sys.D);
sys_init3.Ts = 0;

sys_init3.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [0,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
sys_init3.Structure.B.Free = [0;0;1;1];
sys_init3.Structure.C.Free = zeros(2,4);
sys_init3.Structure.D.Free = [0;0];

opt.Regularization.Lambda = 0;
%opt.Regularization.Nominal = 'zero';
%opt.WeightingFilter = [0,6];

struct_sys = pem(training_data, sys_init3,opt);

% figure()
% impulse(struct_sys)
% 
% figure()
% resid(training_data, struct_sys,ropt)

[struct_sys.A, struct_sys.B]

%%% Validation (for every test do two runs)

load('../Data/val4rad alpha.mat');   % loading alpha's
load('../Data/val4rad theta.mat');   % loading theta's
load('../Data/val4rad input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 2001; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end)];
%ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);     

validation_data = [ymeas,uin];

figure()
compare(validation_data,sys, struct_sys)

%figure()
%resid(validation_data, struct_sys,ropt)
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

%% Computing the terminal set with mpt3
close all;

mpc_A = str_discr_sys.A(1:4,1:4);
mpc_B = str_discr_sys.B(1:4);
mpc_C = eye(4);
mpc_D = zeros(4,1);

Q_mpc = diag([1e1, 1e-1, 1e-3,1e-3]);
R_mpc = 1e-1;

[P,K,~] = idare(str_discr_sys.A,str_discr_sys.B,Q_mpc,R_mpc);

% mptopt('qpsolver', 'quadprog');
% mptopt('lpsolver','LCP')

mpt3_model = ss(str_discr_sys.A,str_discr_sys.B,str_discr_sys.C,str_discr_sys.D,dt);
mpc_mpt3 = LTISystem(mpt3_model);
mpc_mpt3.x.min = [-pi/2; -inf; -inf; -inf];
mpc_mpt3.x.max = [pi/2; inf; inf; inf];
mpc_mpt3.u.min = -1;
mpc_mpt3.u.max = 1;
mpc_mpt3.x.with('reference');
mpc_mpt3.x.reference = 'free';

mpc_mpt3.x.penalty = QuadFunction(Q_mpc);
mpc_mpt3.u.penalty = QuadFunction(R_mpc);
Tset = mpc_mpt3.LQRSet();
PN = mpc_mpt3.LQRPenalty;
mpc_mpt3.x.with('terminalSet');
mpc_mpt3.x.terminalSet = Tset;
mpc_mpt3.x.with('terminalPenalty');
mpc_mpt3.x.terminalPenalty = QuadFunction(P);

mpt_horizon = 40;
ctrl = MPCController(mpc_mpt3,mpt_horizon);%.toExplicit();

reference = [1;0;0;0];

loop = ClosedLoop(ctrl, mpc_mpt3);
x0 = [-1; 0; 0; 0];
Nsim = round(mpt_horizon*2);
data = loop.simulate(x0, Nsim, 'x.reference', reference);

u_c = [1;-1];
y_c = [pi/2;-pi/2];

figure()
subplot(2,1,1)
plot(1:Nsim, data.Y);
hold on;
plot(1:Nsim, y_c*ones(1, Nsim), 'k--')
hold off;
title('outputs')
legend('alpha','theta')
grid on
subplot(2,1,2)
plot(1:Nsim, data.U);
hold on;
plot(1:Nsim, u_c*ones(1, Nsim), 'k--')
hold off
title('inputs')
grid on

H = Tset.A;
gamma = Tset.b;

% Export to a (faster) MPC controller
%ctrl.optimizer.toMatlab('mycontroller.m', 'primal', 'obj');

% precompile optimizer and using YALMIP directly for speed
optim = ctrl.toYALMIP();
constr = optim.constraints;
obj    = optim.objective;
vars   = optim.variables;
params = [vars.x(:,1); vars.filters.x.reference];
decision = vars.u;
fast_opt = optimizer(constr, obj, sdpsettings('solver','quadprog'), params, decision);



% function u0 = get_mu(ctrl,x,reference)
%     u0 = ctrl.evaluate(x,'x.reference', reference);
% end

% function u0 = get_mu(optim,x,reference)
%     u0 = optim{[x; reference]};
% end
%% TEST
% tic
% mu = get_mpc_u(loop,x0,mpt_horizon,reference)
% toc


% profile on;
% mu = get_mu(ctrl,x0,reference)
% profile off;
% profile viewer;

tic
for i=1:100
    mu = get_mu(fast_opt,x0,reference); % fast enough!
end
toc

figure()
plot(mu)
hold on
plot(data.U, '--')
hold off
legend('fast mpc','original mpc')
grid on

%% Testing fast mpc
close all;

x_his = zeros(4,Nsim);
y_his = zeros(2,Nsim);
u_his = zeros(1,Nsim);

x_prev=x0;
for i=1:Nsim
    mu = get_mu(fast_opt,x_prev,reference);
    x_next = mpt3_model.A * x_prev + mpt3_model.B*mu(1);
    x_his(:,i) = x_prev;
    u_his(i) = mu(1);
    y_his(:,i) = str_discr_sys.C * x_prev;

    x_prev = x_next;
end

lqr = mpt3_model;
lqr.A = mpt3_model.A - mpt3_model.B * K;
DC_gain = dcgain(lqr);
lqr.B = mpt3_model.B / DC_gain(1);

tsim = 0:dt:Nsim;
[y_lqr,~,x_lqr] = lsim(lqr,ones(Nsim/dt+1,1),tsim,x0);
u_lqr = -K*x_lqr.' + ones(Nsim/dt+1,1)/DC_gain(1);

figure()
plot(1:Nsim,y_his)
hold on
plot(1:Nsim, y_lqr(1:Nsim,:).', '--')
plot(1:Nsim, y_c*ones(1, Nsim), 'k--')
hold off
legend('fast mpc alpha','fast mpc theta','original mpc alpha','original mpc theta')
grid on

figure()
plot(1:Nsim,u_his)
hold on
plot(1:Nsim, u_lqr(:,1:Nsim),'r--')
plot(1:Nsim, u_c*ones(1, Nsim), 'k--')
hold off
legend('fast mpc control inputs','lqr control inputs')
grid on

%% Augmenting state space for disturbance rejection

A_aug = [str_discr_sys.A,str_discr_sys.B;
        zeros(1,n), 1];
B_aug = [str_discr_sys.B;0];
C_aug = [str_discr_sys.C, zeros(2,1)];
D_aug = str_discr_sys.D;

Extr_states = [eye(4),zeros(4,1)];
Kw = [0,0,0,0,1];

Qk = blkdiag((1e-2)*eye(4), 1e-7);
Rk = 1e-12;

x0 = [0;0;0;0;0];
%% Deactivate disturbance rejection
% A_aug = str_discr_sys.A;
% B_aug = str_discr_sys.B;
% C_aug = str_discr_sys.C;
% D_aug = str_discr_sys.D;
% 
% Extr_states = [eye(4)];
% Kw = [0,0,0,0];
% 
% Qk = (1e-1)*eye(4);
% Rk = 1e-12;
% 
% x0 = [0;0;0;0];

%%
% %% Creating the nonlinear multistage mpc object
% 
% mpc_model = mpt3_model;
% 
% 
% % creating mpcobject
% MPC1 = nlmpcMultistage(mpc_model,dt);
% 
% % setting Q,R and N
% MPC1.PredictionHorizon = 40;
% MPC1.Weights.OutputVariables = [1e1, 1e-1, 1e-3,1e-3];
% MPC1.Weights.ManipulatedVariables = 1e-1;
% 
% % Setting the terminal set
% setterminal(MPC1, struct('A', H, 'b', gamma));
% 
% % Setting the terminal cost
% P = idare(mpc_model.A, mpc_model.B,Q_mpc, R_mpc);
% MPC1.TerminalWeight = P;

% %% implementing the control
% '-----'
% Tsim = 10;
% 
% mdl = 'mpc_test';
% load_system(mdl);
% set_param(mdl, 'FastRestart', 'on');
% set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(h), 'StopTime', num2str(h));
% 
% xk = x0;   % initial state
% X = zeros(length(x0), Tsim);
% U = zeros(1, Tsim);
% 
% for k = 1:Tsim
%     tic
%     uk = -K*xk;   % your YALMIP function
% 
%     % Store
%     X(:,k) = xk;
%     U(:,k) = uk;
% 
%     simIn = Simulink.SimulationInput(mdl);
%     simIn = simIn.setVariable('u', uk);  % set Constant block value
%     %simIn = simIn.setInitialState(struct('x', xk));  % optional if plant is stateful
% 
%     out = sim(simIn);
% 
%     % Get next state from Outport
%     %xk = out.yout.getElement('x_out').Values.Data(end,:)';
%     toc
% end
% 
% set_param(mdl, 'FastRestart', 'off');

%% Functions

function u0 = get_mu(fast_opt,x,reference)
    mu = fast_opt{[x; reference]};
    u0=mu(1);
end
%% MPC() synthesis
% h=dt;
% 
% mpc_A = str_discr_sys.A(1:4,1:4);
% mpc_B = str_discr_sys.B(1:4);
% mpc_C = eye(4);
% mpc_D = zeros(4,1);
% 
% mpc_model = ss(mpc_A,mpc_B,mpc_C,mpc_D,dt);
% 
% 
% % creating mpcobject
% MPC1 = mpc(mpc_model,dt);
% 
% % setting Q,R and N
% MPC1.PredictionHorizon = 40;
% MPC1.Weights.OutputVariables = [1e1, 1e-1, 1e-3,1e-3];
% MPC1.Weights.ManipulatedVariables = 1e-1;
% 
% % Setting the terminal set
% setterminal(MPC1, struct('A', H, 'b', gamma));
% 
% % Setting the terminal cost
% P = idare(mpc_model.A, mpc_model.B,Q_mpc, R_mpc);
% MPC1.TerminalWeight = P;

%% Setting up terminal set AND terminal cost by expanding output vector

% Q_mpc = diag(MPC1.Weights.OutputVariables);
% R_mpc = MPC1.Weights.ManipulatedVariables;
% 
% P = idare(mpc_model.A, mpc_model.B,Q_mpc, R_mpc);
% 
% Qc = chol(P, 'lower'); %this will eventually enforce the terminal cost
% mpc_C_aug = [mpc_model.C;Qc];
% mpc_D_aug = [mpc_model.D,zeros(size(mpc_model.D,1), 1)];
% mpc_model_aug = ss(mpc_model.A,mpc_model.B,mpc_C_aug,mpc_D_aug,dt);
% 
% % creating mpcobject
% MPC_aug = mpc(mpc_model_aug,dt); %is matlab discretizing this again automatically?
% 
% % setting Q,R and N
% MPC1.PredictionHorizon = 40;
% MPC1.Weights.OutputVariables = [1e1, 1e-1, 1e-3,1e-3];
% MPC1.Weights.ManipulatedVariables = 1e-1;
% Y = struct('Weight', P); % Vf --> x^T P x 
% U = struct();
% 
% Y = struct('Weight', [zeros(1, 4), ones(1, 4)], 'Min', [], 'Max', []);
% setterminal(mpcobj, Y, []);

%% Setting up terminal cost P (x^T (beta*P) x < alpha) by expanding output vector
% beta = 1e2;
% 
% Q_mpc = diag(MPC1.Weights.OutputVariables);
% R_mpc = MPC1.Weights.ManipulatedVariables;
% 
% P = idare(mpc_model.A, mpc_model.B,Q_mpc, R_mpc);
% 
% Qc = chol(P, 'lower'); %this will eventually enforce the terminal cost
% mpc_C_aug = [mpc_model.C;Qc];
% mpc_D_aug = [mpc_model.D,zeros(size(mpc_model.D,1), 1)];
% mpc_model_aug = ss(mpc_model.A,mpc_model.B,mpc_C_aug,mpc_D_aug,dt);
% 
% % creating mpcobject
% MPC_aug = mpc(mpc_model_aug,dt); %is matlab discretizing this again automatically?
% 
% % setting Q,R and N
% MPC1.PredictionHorizon = 40;
% MPC1.Weights.OutputVariables = [1e1, 1e-1, 1e-3,1e-3];
% MPC1.Weights.ManipulatedVariables = 1e-1;
% Y = struct('Weight', beta*P); % Vf --> beta*Vf = x^T (beta * P) x 
% U = struct();
% 
% Y = struct('Weight', [zeros(1, 4), ones(1, 4)], 'Min', [], 'Max', []);
% setterminal(mpcobj, Y, []);
% omitting the build in state estimation
% setEstimator(MPC1, "custom")
% MPC1.Model.Plant = setmpcsignals(MPC1.Model.Plant, 'MO', []);
% Setting up terminal cost P
%P = ....
%Y.Weight = P;
