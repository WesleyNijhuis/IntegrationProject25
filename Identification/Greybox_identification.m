close all; clear; clc;

%% Loading acquired data
training_data_y = cell(1,2);
training_data_u = cell(1,2);

%training set 2: square sweep
load('../Data/prbs 1 alpha.mat');   % loading alpha's
load('../Data/prbs 1 theta.mat');   % loading theta's
load('../Data/prbs 1 input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 501;
data_begin = 1;
%ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end) - mean(theta(data_begin:data_end))];
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
%uin = u(data_begin:data_end,2) - mean(u(data_begin:data_end,2));   
uin = u(data_begin:data_end,2);   

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
%% Setting system parameters for simulink
h=0.01;
epsilon=0.00242;
%% graybox identification
close all;
% Furuta Pendulum Grey-box Model in Downward Position

% Dynamics based on linearized model from the paper:
% "On the Dynamics of the Furuta Pendulum" (Cazzolato & Prime)
% Section 6.2 - Downward linearized model (sign flipped entries)
L1 = 0.085;
L2 = 0.129;
m1 = 0.095;
m2 = 0.024;
b1 = 1e-4;
b2 = 1e-7;
Kt = 0.042;
Rm = 8.4;

J1 = 0.095*0.085^2/12+4.6e-6;
J2 = 0.024*0.129^2/12;
params = [J1,J2, m1, m2, L1/2, L2/2, L1, b1, b2,Kt,Rm];

init_sys = idgrey(@furuta_downward_model, params, 'c');
itit_sys.x0 = [0;0;0;0];
training_data = iddata(training_data_y{2},training_data_u{2},dt);
opt = greyestOptions;
opt.InitialState = 'estimate';
opt.WeightingFilter = [1,10];
opt.Display = 'off';
opt.SearchMethod = 'gna';
opt.SearchOptions.Tolerance = 1e-8;
opt.SearchOptions.MaxIterations = 4000;
opt.Regularization.Lambda = 1e7;
opt.Regularization.Nominal = 'model';
estimated_sys = greyest(training_data, init_sys,opt);
compare(training_data, estimated_sys,init_sys);

%  validation
load('../Data/val4rad alpha.mat');   % loading alpha's
load('../Data/val4rad theta.mat');   % loading theta's
load('../Data/val4rad input.mat');   % loading inputs
alpha = alpha(:,2);
theta = theta(:,2);

data_end = 2000;
data_begin = 1;
%ymeas = [alpha(data_begin:data_end) - mean(alpha(data_begin:data_end)), theta(data_begin:data_end) - mean(theta(data_begin:data_end))];
ymeas = [alpha(data_begin:data_end), theta(data_begin:data_end)];
%uin = u(data_begin:data_end,2) - mean(u(data_begin:data_end,2));   
uin = u(data_begin:data_end,2);   

validation_data_y = ymeas;
validation_data_u = uin;

validation_data = iddata(validation_data_y,validation_data_u,dt);
figure()
compare(validation_data, estimated_sys,init_sys);

estimated_sys

function [A,B,C,D] = furuta_downward_model(p,Ts)
    % Parameters vector:
    J1 = p(1);
    J2 = p(2);
    m1 = p(3);
    m2 = p(4);
    l1 = p(5);
    l2 = p(6);
    L1 = p(7);
    b1 = p(8);
    b2 = p(9);
    Km = p(10);
    Rm = p(11);
    g  = 9.81;  % gravitational constant [m/s^2]

    % Compute equivalent inertias using parallel axis theorem
    J0_hat = J1 + m1*l1^2 + m2*L1^2;  % Arm 1 + rotor + pendulum's arm contribution
    J2_hat = J2 + m2*l2^2;            % Pendulum's own moment

    % Common denominator
    denom = J0_hat * J2_hat - m2^2 * L1^2 * l2^2;

    % A matrix (state matrix)
    A = zeros(4);
    A(1,3) = 1;
    A(2,4) = 1;
    A(3,2) = g * m2^2 * l2^2 * L1 / denom;          
    A(3,3) = -b1 * J2_hat / denom;
    A(3,4) = -b2 * m2 * l2 * L1 / denom;               
    A(4,2) = g * m2 * l2 * J0_hat / denom;        
    A(4,3) = -b1 * m2 * l2 * L1 / denom;            
    A(4,4) = -b2 * J0_hat / denom;

    % B matrix (input matrix) - single input assumed
    B = zeros(4,1);
    B(3) = J2_hat / denom;
    B(4) = m2 * l2 * L1 / denom;        

    % Tau to u
    A(3,4) = A(3,4) - B(3)*(Km^2/Rm);
    A(4,4) = A(4,4) - B(4)*(Km^2/Rm);
    B(3) = B(3)*(5*Km/Rm);
    B(4) = B(4)*(5*Km/Rm);

    % To downwards position
    A(3,4) = -A(3,4);
    A(4,2) = -A(4,2);
    A(4,3) = -A(4,3);
    B(4) = -B(4);

    % C matrix (output matrix) - only angular positions are measured
    C = [1 0 0 0;
         0 1 0 0];

    % D matrix (direct feedthrough)
    D = zeros(2,1);
end
