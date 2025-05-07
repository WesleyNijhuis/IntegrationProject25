close all; clear; clc;

%% Loading acquired data
load('../Data/Sweep0.001to7 alpha.mat');   % loading alpha's
load('../Data/Sweep0.001to7 theta.mat');   % loading theta's
load('../Data/Sweep0.001to7 input.mat');   % loading inputs

data_end = 10000; %for debugging
data_begin = 1;
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
u = uin(data_begin:data_end).';            
dt = 0.01;
t = dt*(1:1:size(u,1)).';

plot(u)
figure
plot(ymeas(:,1))
figure
plot(ymeas(:,2))

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

% initial guess (from paper), still need to guess damping
b1 = 0.0015; 
b2 = 0.0005;
theta_init =  [-38.26;-151.88;-2548.8*b1;2519.2*b1;-36.21;-2519.2*b2;-10001*b2;107.05;105.8;-7241;862];

% Set 1
[A0, B0, C0, D0, x00] = theta2matrices(theta_init);


% BIG TODO: CHANGE THE SS TO DISCRETIZED, DIFFERENT STRUCTURE?

% PEM DT
%K = zeros(4,1); % no kalman observer for now: add later
%T_s = 0.01; % sampling time
init_sys = idss(A0, B0, C0, D0); %x0 is 0
%init_sys.x0 = x00;
init_sys.Ts = 0;
init_sys.Structure.A.Free = [[0,0,0,0,0];
                             [0,0,0,0,0];
                             [0,1,1,1,1];
                             [0,1,1,1,1];
                             [0,0,1,0,1]]; % only the parameter entries (1 or 'True') can be changed
init_sys.Structure.B.Free = [0;0;0;0;1];
init_sys.Structure.C.Free = zeros(2,5);
init_sys.Structure.D.Free = 0;

training_data = [ymeas,u];
opt = ssestOptions('Display','on','SearchMethod','lm');
opt.SearchOptions.MaxIterations = 100000;
opt.SearchOptions.Tolerance = 0.01;
%opt.InitialState = 'estimate';
sys = pem(training_data, init_sys,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
C = sys.C
D = sys.D
x0 = sys.x0

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