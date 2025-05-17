close all; clear; clc;

%% Loading acquired data
load('../Data/autonomous 1 alpha.mat');   % loading alpha's
load('../Data/autonomous 1 theta.mat');   % loading theta's

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 18000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
uin = zeros(size(alpha));            
dt = 0.01;
t = dt*(1:1:size(uin,1)).';

%Cutting the total dataset into subexperiments (here theta only starting from
%~+3.14)
Et_start = [524,1435,4212,6559,7740,10087,11180] + 500; % +550 to get rid of nonlinear domain
Et_end = [1275,2188,5089,7424,8620,10935,12051];

plot(uin)
title('uin')
figure
plot(ymeas(:,1))
title('outputs')
hold on
plot(ymeas(:,2))
hold off
legend('alpha','theta')
grid on

%% Initial guess can be entered in the function theta2matrices!
function [Ac, Bc, Cc, Dc] = initialguess(m1,m2,l1,l2,L1,J1,J2,Lm,Km,Rm,b1,b2,KT)
g = 9.81;

J1_hat = J1 + m1 * l1^2;
J2_hat = J2 + m2 * l2^2;

J0_hat = J1_hat + m2 * L1^2;

A41 = 0;
A42 = (g * m2 * l2 * J0_hat) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A43 = (-b1 * m2 * l2 * L1) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A44 = (-b2 * J0_hat) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);

A31 = -KT; %Torsion induced by cable
A32 = (g * (m2^2) * (l2^2) * L1) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A33 = (-b1 * J2_hat) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A34 = (-b2 * m2 * l2 * L1) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2); 

%switching signs for downwards position
A34 = -A34;
A42 = -A42;
A43 = -A43;

Ac = [[0,0,1,0];
     [0,0,0,1];
     [A31,A32,A33,A34];
     [A41,A42,A43,A44]];

Bc = [0;0;0;0];

Cc = [[1,0,0,0];
      [0,1,0,0]];

Dc = [0;0];
end
%% Splitting into training/validation set
%TODO: in our case the training set and validation cant be split as we're
%doing a frequency sweep
% f_t = 0.07;                 % fraction of data dedicated to training
% ut = u(1:round(end*f_t));            % input for train-set
% uv = u(round(end*f_t)+1:end);          % input for validation-set
% yt = ymeas(1:round2even(end*f_t));        % output for train-set
% yv = ymeas(round2even(end*f_t)+1:end);      % output for validation-set
%% Run Algorithm

% load training data
training_data_y = cell(1,size(Et_start,2));
training_data_u = cell(1,size(Et_start,2));

for i=1:size(Et_start,2)
    alpha_exp = ymeas(Et_start(i):Et_end(i),1);
    theta_exp = ymeas(Et_start(i):Et_end(i),2) - mean(ymeas(Et_start(i):Et_end(i),2)); % remove if unnecessary
    training_data_y{i} = [alpha_exp,theta_exp];
    training_data_u{i} = uin(Et_start(i):Et_end(i),1);
end

training_data = iddata(training_data_y,training_data_u,dt);
%% Initial guess Autonomous system
% initial guess (from data sheet and paper), still need to guess damping
b1 = 0.0015; 
b2 = 0.0005;

m1 = 0.095; %kg
m2 = 0.024; %kg
l1 = 0.085/2; %(m)
l2 = 0.129/2; %Lp/2 (m) 
L1 = 0.085; %(m)
L2 = 0.129 - 0.01; %(m)
J1 = 0.6*1e-6 + 4*1e-6 + m1*(L1^2)/12; %Jr+Jm (kg/m2) %Jp and Ja SHOULD BE AROUND THE CENTER!
J2 = m2*((L2)^2)/12; %(kg/m2)
Lm = 1.16 * 1e-3; %mH
Km = 0.042; %0.042(Nm/A)
Rm = 8.4; %(ohm)
KT = 1;

[A0,B0,C0,D0] = initialguess(m1,m2,l1,l2,L1,J1,J2,Lm,Km,Rm,b1,b2,KT);
% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;
init_sys.x0 = [-0.09;0.131;-1;3];

% testing the initial guess
T = dt*(Et_start:1:Et_end).';
ysim = lsim(init_sys,training_data_u{1}, T, init_sys.x0);

close all
figure('Position',[100 100 1000 400])
sgtitle('initial guess')
subplot(1,2,1)
plot(ysim(:,1));
hold on
plot(training_data_y{1}(:,1));
hold off
title('alpha')
legend('simulated','data')

subplot(1,2,2)
plot(ysim(:,2));
hold on
plot(training_data_y{1}(:,2));
hold off
title('theta')
legend('simulated','data')


%% PEM Autonomous system
% Setting some bounds for these parameters (we expect the initial guess of
% each parameter to be within CF*100% of initial guess (But if we're
% actually close, this shouldnt matter)

% CF = 0.5;
% init_sys.Structure.A.Minimum = A0 - abs(A0).*CF - 1000*[0,0,0,0,0;0,0,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0]; % relax cons for A31
% init_sys.Structure.A.Maximum = A0 + abs(A0).*CF + 1000*[0,0,0,0,0;0,0,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0]; % relax cons for A31
% init_sys.Structure.B.Minimum = B0 - abs(B0).*CF;
% init_sys.Structure.B.Maximum = B0 + abs(B0).*CF;

% Setting which entries are the parameters
init_sys.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [1,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
init_sys.Structure.B.Free = [0;0;0;0];
init_sys.Structure.C.Free = zeros(2,4);
init_sys.Structure.D.Free = [0;0];


opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-8;
opt.InitialState = 'estimate';
opt.EnforceStability = true;
%opt.OutputOffset = [mean(ymeas(:,1));0];
sys_aut = pem(training_data, init_sys,opt);

compare(training_data,init_sys,sys)

%% Loading acquired data full state space
load('../Data/Sweep 6 alpha.mat');   % loading alpha's
load('../Data/Sweep 6 theta.mat');   % loading theta's
load('../Data/Sweep 6 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 8000; %for debugging
data_begin = 1000;
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

%% Initial for full system using A0 from autonomous system identification
function [Ac, Bc, Cc, Dc] = initialguess_full(A0,m1,m2,l1,l2,L1,J1,J2,Km,Rm)
g = 9.81;

J1_hat = J1 + m1 * l1^2;
J2_hat = J2 + m2 * l2^2;

J0_hat = J1_hat + m2 * L1^2;

B31 = J2_hat / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B41 = (m2 * L1 * l2) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);

A41 = A0(4,1);
A42 = A0(4,2);
A43 = A0(4,3);
A44 = A0(4,4);

A31 = A0(3,1); %Torsion induced by cable
A32 = A0(3,2);
A33 = A0(3,3);
A34 = A0(3,4); 

%switching signs for downwards position
B41 = -B41;

% Adding extended state parameters
A33 = A33 + B31*Km^2/Rm;
A43 = A43 + B41*Km^2/Rm;

B31 = B31 * 5 * Km / Rm;
B41 = B41 * 5 * Km / Rm;

Ac = [[0,0,1,0];
     [0,0,0,1];
     [A31,A32,A33,A34];
     [A41,A42,A43,A44]];

Bc = [0;0;B31;B41];

Cc = [[1,0,0,0];
      [0,1,0,0]];

Dc = [0;0];
end
%% Splitting into training/validation set
%TODO: in our case the training set and validation cant be split as we're
%doing a frequency sweep
% f_t = 0.07;                 % fraction of data dedicated to training
% ut = u(1:round(end*f_t));            % input for train-set
% uv = u(round(end*f_t)+1:end);          % input for validation-set
% yt = ymeas(1:round2even(end*f_t));        % output for train-set
% yv = ymeas(round2even(end*f_t)+1:end);      % output for validation-set
%% Run Algorithm

% load training data
training_data = [ymeas,uin];

% additional parameters for B

Lm = 1.16 * 1e-3; %mH
Km = 0.042; %0.042(Nm/A)
Rm = 8.4; %(ohm)

% Using parameters acquired from autonomous experiments!
[Af0,Bf0,Cf0,Df0] = initialguess_full(sys_aut.A,m1,m2,l1,l2,L1,J1,J2,Km,Rm);

%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(Af0, Bf0, Cf0, Df0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;

% testing the initial guess
ysim = lsim(init_sys,uin,t);

close all
figure('Position',[100 100 1000 400])
sgtitle('initial guess')
subplot(1,2,1)
plot(ysim(:,1));
hold on
plot(ymeas(:,1));
hold off
title('alpha')
legend('simulated','data')

subplot(1,2,2)
plot(ysim(:,2));
hold on
plot(ymeas(:,2));
hold off
title('theta')
legend('simulated','data')


%% PEM
% Setting some bounds for these parameters (we expect the initial guess of
% each parameter to be within CF*100% of initial guess (But if we're
% actually close, this shouldnt matter)

% CF = 0.5;
% init_sys.Structure.A.Minimum = A0 - abs(A0).*CF - 1000*[0,0,0,0,0;0,0,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0]; % relax cons for A31
% init_sys.Structure.A.Maximum = A0 + abs(A0).*CF + 1000*[0,0,0,0,0;0,0,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0]; % relax cons for A31
% init_sys.Structure.B.Minimum = B0 - abs(B0).*CF;
% init_sys.Structure.B.Maximum = B0 + abs(B0).*CF;

% Setting which entries are the parameters
init_sys.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [0,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
init_sys.Structure.B.Free = [0;0;1;1];
init_sys.Structure.C.Free = zeros(2,4);
init_sys.Structure.D.Free = [0;0];


opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 4000;
opt.SearchOptions.Tolerance = 1e-8;
opt.InitialState = 'estimate';
opt.EnforceStability = false;
opt.OutputOffset = [mean(ymeas(:,1));0];
sys = pem(training_data, init_sys,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
C = sys.C
D = sys.D
x0 = sys.x0

compare(training_data,init_sys,sys)
% ypred = simsystem(Abar, Bbar, C, D, x0, uv);
% ypred = ypred(:);
% 
% RMSE = rmse(ypred, yv);
% VAF = 1 - var(yv - ypred)/var(yv);
% 
% figure('position', [0, 0, 800, 400])  % create new figure with specified size  
% 
% plot(ypred)
% title(['Resuls for theta1, training set 1 with RMSE of ', num2str(RMSE), ' and VAF of ', num2str(VAF)])
% 
% hold on
% plot(yv)
% legend('ypred', 'yv')
% hold off 
%% Parameters to state matrices
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