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

%% Initial guess can be entered in the function theta2matrices!
function [Ac, Bc, Cc, Dc] = initialguess(m1,m2,l1,l2,L1,J1,J2,Lm,Km,Rm,b1,b2,KT)
g = 9.81;

J1_hat = J1 + m1 * l1^2;
J2_hat = J2 + m2 * l2^2;

J0_hat = J1_hat + m2 * L1^2;

B31 = J2_hat / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B41 = (m2 * L1 * l2) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B32 = (m2 * L1 * l2) / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);
B42 = J0_hat / (J0_hat * J2_hat - m2^2 * L1^2 * l2^2);

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
B32 = -B32;
B41 = -B41;

% Adding extended state parameters
A33 = A33 + B31*Km^2/Rm;
A43 = A43 + B41*Km^2/Rm;

Ac = [[0,0,1,0];
     [0,0,0,1];
     [A31,A32,A33,A34];
     [A41,A42,A43,A44]];

Bc = 5*[0;0;B31;B41];

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
KT = 3;

[A0,B0,C0,D0] = initialguess(m1,m2,l1,l2,L1,J1,J2,Lm,Km,Rm,b1,b2,KT);
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
init_sys.Structure.A.Free = [[0,0,0,0];
                             [0,0,0,0];
                             [1,1,1,1]; %A31 is a torsional factor induced by cable
                             [0,1,1,1]]; % only the parameter entries (1 or 'True') can be changed
init_sys.Structure.B.Free = [0;0;1;1];
init_sys.Structure.C.Free = zeros(2,4);
init_sys.Structure.D.Free = [0;0];

BIG = 1e6;
init_sys.Structure.A.Maximum = [0,0,1,0;0,0,0,1;0,BIG,0,BIG;0,0,0,0];
init_sys.Structure.A.Minimum = [0,0,1,0;0,0,0,1;-BIG,0,-BIG,0;0,-BIG,-BIG,-BIG];
init_sys.Structure.B.Maximum = [0,0,BIG,0];
init_sys.Structure.B.Minimum = [0,0,0,-BIG];

opt = ssestOptions('Display','on','SearchMethod','auto');
opt.SearchOptions.MaxIterations = 2000;
opt.SearchOptions.Tolerance = 1e-10;
opt.InitialState = 'zero';
opt.EnforceStability = false;
opt.Regularization.Lambda = 1e-10;
opt.Regularization.Nominal = 'model'; % restricting parameters to stay close to initial guess (using a priori knowledge)
opt.OutputOffset = [mean(ymeas(:,1));mean(ymeas(:,2))];
opt.Advanced.DDC = 'on';
sys = pem(training_data, init_sys,opt);

disp('Results theta 1, set 1')
Abar = sys.A
Bbar = sys.B
C = sys.C
D = sys.D
x0 = sys.x0

compare(training_data,init_sys,sys)
