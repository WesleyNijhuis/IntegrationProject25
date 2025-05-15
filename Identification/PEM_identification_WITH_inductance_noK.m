close all; clear; clc;

%% Loading acquired data
load('../Data/Sweep 6 alpha.mat');   % loading alpha's
load('../Data/Sweep 6 theta.mat');   % loading theta's
load('../Data/Sweep 6 input.mat');   % loading inputs

alpha = alpha(:,2);
theta = theta(:,2);

data_end = 6000; %for debugging
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

B31 = J2_hat / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B41 = (m2 * L1 * l2) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B32 = (m2 * L1 * l2) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B42 = J0_hat / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);

A31 = -KT; %Torsion induced by cable
A32 = (g * (m2^2) * (l2^2) * L1) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A33 = (-b1 * J2_hat) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
A34 = (-b2 * m2 * l2 * L1) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2); 

%switching signs for downwards position
A34 = -A34;
A42 = -A42;
A43 = -A43;
B32 = -B32;
B41 = -B41;

% Adding extended state parameters
A35 = B31*Km;
A45 = B41*Km;

A53 = -Km/Lm;
A55 = -Rm/Lm;

B5 = 5*1/Lm; % Voltage sent is multiplied by 5 in Qube block!

% Theta is defined positive in CW direction?
A32 = -A32;
A34 = -A34;
A43 = -A43;
A45 = -A45;

Ac = [[0,0,1,0,0];
     [0,0,0,1,0];
     [A31,A32,A33,A34,A35];
     [A41,A42,A43,A44,A45];
     [0,0,A53,0,A55]];

Bc = [0;0;0;0;B5];

Cc = [[1,0,0,0,0];
      [0,1,0,0,0]];

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

% initial guess (from data sheet and paper), still need to guess damping
b1 = 0.0015; 
b2 = 0.0005;

m1 = 0.095; %kg
m2 = 0.024; %kg
l1 = 0.085/2; %(m)
l2 = 0.129/2; %Lp/2 (m) 
L1 = 0.085; %(m)
L2 = 0.129; %(m)
J1 = 0.6*1e-6 + 4*1e-6 + m1*(L1^2)/12; %Jr+Jm (kg/m2) %Jp and Ja SHOULD BE AROUND THE CENTER!
J2 = m2*(L2^2)/12; %(kg/m2)
Lm = 1.16 * 1e-3; %mH
Km = 0.042; %0.042(Nm/A)
Rm = 8.4; %(ohm)
KT = 1;

[A0,B0,C0,D0] = initialguess(m1,m2,l1,l2,L1,J1,J2,Lm,Km,Rm,b1,b2,KT);

% Set 1
theta_init =  [-38.26;-151.88;-2548.8*b1;-2519.2*b1;-36.21;-2519.2*b2;-10001*b2;107.05;105.8;-7241;862];
theta_init = rand(11,1)-5*rand(11,1);
[A02, B02, C02, D02, x00] = theta2matrices(theta_init);

% BIG TODO: CHANGE THE SS TO DISCRETIZED, DIFFERENT STRUCTURE?

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
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
init_sys.Structure.A.Free = [[0,0,0,0,0];
                             [0,0,0,0,0];
                             [1,1,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1,1];
                             [0,0,1,0,1]]; % only the parameter entries (1 or 'True') can be changed
init_sys.Structure.B.Free = [0;0;0;0;1];
init_sys.Structure.C.Free = zeros(2,5);
init_sys.Structure.D.Free = [0;0];


opt = ssestOptions('Display','on','SearchMethod','gna');
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-8;
opt.InitialState = 'zero';
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