close all; clear; clc;

%% Loading acquired data
load('../Data/Sweep 5 alpha.mat');   % loading alpha's
load('../Data/Sweep 5 theta.mat');   % loading theta's
load('../Data/Sweep 5 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 12000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);            
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

% % Removing deadzone induced offset
% dz_right = 0.0028;
% dz_left = -0.0024;
% 
% for i=1:size(uin,1)
%     if uin(i,1) < dz_right
%         if uin(i,1) > 0
%             uin(i,1) = 0;
%         elseif uin(i,1) > dz_left
%             uin(i,1) = 0;
%         end
%     end
% end

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
 plot(sing_vals,'o')
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
theta_init =  [12.2;-13.6;-2.2;-1.4;-36.2;-1.4;-2.2;2.2;4.4;-7241;862];


% Testing different sizes
n=5;
A0 = ones(n,n);
B0 = ones(n,1);
C0 = ones(2,n);
D0 = zeros(2,1);


% BIG TODO: CHANGE THE SS TO DISCRETIZED, DIFFERENT STRUCTURE?

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
init_sys.Ts = 0.01;
%init_sys.x0 = x00;
init_sys.Ts = dt;


training_data = [ymeas,uin];
opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-12;
opt.InitialState = 'zero';
opt.InitializeMethod = 'n4sid';
opt.N4Weight = 'MOESP';
opt.EnforceStability = false;
sys_init2 = n4sid(training_data, init_sys,opt);

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

%% Validation (for every test do two runs)

load('../Data/Sweep 3v alpha.mat');   % loading alpha's
load('../Data/Sweep 3v theta.mat');   % loading theta's
load('../Data/Sweep 3v input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 10000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = u(data_begin:data_end,2);              
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

% % Removing deadzone induced offset
% dz_right = 0.0028;
% dz_left = -0.0024;
% 
% for i=1:size(uin,1)
%     if uin(i,1) < dz_right
%         if uin(i,1) > 0
%             uin(i,1) = 0;
%         elseif uin(i,1) > dz_left
%             uin(i,1) = 0;
%         end
%     end
% end

figure
plot(uin)
title('uin')
figure
plot(ymeas(:,1))
title('outputs')
hold on
plot(ymeas(:,2))
hold off
legend('alpha','theta')

validation_data = [ymeas,uin];
figure
compare(validation_data,sys)

%% (modal) Transformation
sys_m = canon(sys, 'modal)')

step(sys, Config)
hold on
step(sys_m, Config)
hold off
legend()

%% pole placement design
poles = [0.8,0.81,0.82,0.83,0.84]
K = place(sys_m.A, sys_m.B, poles)

sfsys = sys_m % printing original matrices
sfsys.A = sfsys.A - sfsys.B*K;
DC_gain = dcgain(sfsys);
G = 1/DC_gain(1);
sfsys.B = G*sfsys.B;
%pzplot(lqsys)
step(sfsys)
grid on
stepinfo(sfsys).TransientTime
%[yout,tout]


%% LQR design
q1 = 200;
q2 = 40;
q3 = 70;
Q = diag([q1, q1, q2, q2, q3]);
R = [0.001];
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

%% Transform system to CT for parameters

sys_c = d2c(sys,'zoh');
impulse(sys_c)
hold on
impulse(sys)
hold off
legend()


% %%%%%%%%%% Everything after this is not working correctly (yet) %%%%%%%%
% 
% %% Find a transformation matrix T to go from blackbox to structured parameter matrices
% O = obsv(sys_c.A, sys_c.C);
% Ostar = O(1:4,1:5);
% Ostar_t = [[1,0,0,0,0];
%           [0,1,0,0,0];
%           [0,0,1,0,0];
%           [0,0,0,1,0]];
% 
% H = kron(eye(5), Ostar_t);
% f = Ostar(:);
% 
% Aeq = zeros(4, 25);
% for i = 1:4
%     Aeq(i, (i-1)*5 + (1:5)) = sys_c.B';  % row i of T times B
% end
% 
% beq = zeros(4,1);
% 
% options = optimoptions('lsqlin', 'Display', 'iter', ...
%     'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-15);
% Tentries = lsqlin(H, f, [], [], Aeq, beq);
% T = reshape(Tentries,5,5);
% 
% %% Give structured parameter matrices
% A = T*Abar/T
% B = T*Bbar
% C = Cbar/T
% 
% sys2_c = ss(A,B,C,D)
% impulse(sys2_c)
% hold on
% impulse(sys)
% hold off
% legend()
% 
% % Note: this gives state matrices for the sampled system
% 
% 
% %% Estimating CT structured parameters
% 
% % B(1:4) = zeros(4,1);
% % A(:,1) = zeros(5,1);
% % A(5,2) = 0;
% % A(5,4) = 0;
% % 
% % init_sys2 = idss(A, B, C, D) 
% % init_sys2.Ts = 0;
% % init_sys2.Structure.A.Free = [[0,0,0,0,0];
% %                              [0,0,0,0,0];
% %                              [0,1,1,1,1];
% %                              [0,1,1,1,1];
% %                              [0,0,1,0,1]]; % only the parameter entries (1 or 'True') can be changed
% % init_sys2.Structure.B.Free = [0;0;0;0;1];
% % init_sys2.Structure.C.Free = zeros(2,5);
% % init_sys2.Structure.D.Free = 0;
% % 
% % opt2 = ssestOptions('Display','on','SearchMethod','auto');
% % opt2.SearchOptions.MaxIterations = 2000;
% % %opt2.InitialState = 'estimate';
% % sys2 = pem(training_data, init_sys2,opt2);
% % 
% % 
% % impulse(sys)
% % hold on
% % impulse(sys2)
% % hold off
% 
% %% Validation
% 
% % ypred = simsystem(Abar, Bbar, C, D, x0, uv);
% % ypred = ypred(:);
% % 
% % RMSE = rmse(ypred, yv);
% % VAF = 1 - var(yv - ypred)/var(yv);
% % 
% % figure('position', [0, 0, 800, 400])  % create new figure with specified size  
% % 
% % plot(ypred)
% % title(['Resuls for theta1, training set 1 with RMSE of ', num2str(RMSE), ' and VAF of ', num2str(VAF)])
% % 
% % hold on
% % plot(yv)
% % legend('ypred', 'yv')
% % hold off 
% %% Parameters to state matrices
function [Abar,Bbar,C,D,x0] = theta2matrices(theta)
%%
% Function INPUT
% theta Parameter vector (vector of size: according to the realization)
%
%
% Function OUTPUT
% Abar System matrix A (matrix of size n x n)
% Bbar System matrix B (matrix of size n x m)
% C System matrix C (matrix of size l x n)
% D System matrix D (matrix of size l x m)
% x0 Initial state (vector of size n x one)

Abar=[[0, 0, 1, 0, 0]; 
    [0, 0, 0, 1, 0]; 
    [0, theta(1), theta(3), theta(6), theta(8)]; 
    [0, theta(2), theta(4), theta(7), theta(9)]; 
    [0, 0, theta(5), 0, theta(10)]]; 

Bbar=[0 ; 0; 0; 0; theta(11)]; 

C=[[1, 0, 0, 0, 0];[0, 1, 0, 0, 0]];

D=[0;0];

x0=[0; 0; 0; 0; 0];
end

%% Simulating systems
function [y, x] = simsystem(A, B, C, D, x0, u)
% Instructions:
% Simulating a linear dynamic system given input u, matrices A,B,C,D ,and
% initial condition x(0)
%
% n = size(A, 1);
% m = size(B, 2);
% l = size(C, 1);
% N = size(u, 1);
%
% Function INPUT
% A system matrix (matrix of size n x n)
% B system matrix (matrix of size n x m)
% C system matrix (matrix of size l x n)
% D system matrix (matrix of size l x m)

% x0 initial state (vector of size n x one)
% u system input (matrix of size N x m)
%
% Function OUTPUT
% x state of system (matrix of size N x n)
% y system output (matrix of size N x l)
x(1,:) = x0';
%%%%%% YOUR CODE HERE %%%%%%
y = zeros(max(size(u,1), size(u,1)), 2);
for k=1:size(u,1)
 x(k+1, :) = x(k, :) * A.' + u(k) * B.'; % Transposed state equation
 y(k,:) = x(k, :) * C.' + u(k) * D.';
end
x=x(1:end -1, :); % removing the last entry to get Nxn x()
end

%% additional functions

function result=round2even(x)
    result = 2*round(x/2);
end