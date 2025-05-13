clc;
clear all;
close all;
%% Set up parameters
br=0.001;
bp=0.00005;
% Resistance
Rm = 8.4;
%Rm = 12;
% Current-torque (N-m/A)
kt = 0.042;
% Back-emf constant (V-s/rad)
km = 0.042;

% % Rotary Arm
% Mass (kg)
mr = 0.095;
% Total length (m)
r = 0.085;

% Moment of inertia about pivot (kg-m^2)
Jr = mr*r^2/3;
%Jr = mr*r^2/12;
% Equivalent Viscous % % Pendulum Link
% Mass (kg)
mp = 0.024;
% Total length (m)
Lp = 0.129;
% Pendulum center of mass (m)
l = Lp/2;
% Moment of inertia about pivot (kg-m^2)
Jp = mp*Lp^2/3;
% Equivalent Viscous Damping Coefficient (N-m-s/rad)
% Gravity Constant
g = 9.81;
%   Copyright 2013 The MathWorks, Inc.
Jt = Jr*Jp - mp^2*r^2*l^2;

h=0.01;

% State Space Representation
A = [0 0 1 0;
    0 0 0 1;
    0 mp^2*l^2*r*g/Jt  -br*Jp/Jt  -mp*l*r*bp/Jt;
    0  mp*g*l*Jr/Jt   -mp*l*r*br/Jt -Jp*bp/Jt];

Bold = [0; 0; Jp/Jt; -mp*l*r/Jt];
C = eye(2,4);
D = zeros(2,1);
% Add actuator dynamics
gain=113;
% gain=113;
B = -km * Bold / Rm*gain;
B(4)=-B(4);
A(3,1) = A(3,1) - km*km/Rm*Bold(3);
A(4,1) = A(4,1) + km*km/Rm*Bold(4);
%% Swing up Control
u_max=10; %Input saturation
ke=0.1; %
%% pole-placement
T=ctrb(A,B);
rank(T)
% control specification
zeta=0.9;%damping ratio
wn=4; %natural frequency
% Location of dominant poles along real-axis
sigma=zeta*wn;
% Location of dominant poles along img axis
wd=wn*sqrt(1-zeta^2); %damped natural frequency
% Desired poles 
DP= [-sigma+1i*wd, -sigma-1i*wd, -20, -60];
% Find control gain using MATLAB pole-placement command
[K, prec]=place(A,B,DP)
eig(A-B*K)
% K=[0.33373 -1.4261 0.08819 -0.12198]
%% closed-loop identification (greyest)
alpha=alphadata.Data;
theta=thetadata.Data;

u=input.Data;
y=[alpha,theta];

data = iddata(y(300:end,:), u(300:end)*0.3,2*h, 'Name', 'Collected Measurement');
data.InputName = 'Voltage';
data.InputUnit = 'V';
data.OutputName = {'Alpha (rad)', 'Theta (rad)'};
data.OutputUnit = {'rad', 'rad'};
data.Tstart = 0;
data.TimeUnit = 's';
% ode1
odefun = 'DCMotorODE_up';
gain=140;
Ts=h;
parameters = {gain};
fcn_type = 'c';
sys = idgrey(odefun,parameters,fcn_type);
sys.Structure.Parameters(1).Minimum =140;
opt=greyestOptions('EnforceStability',true,'Initialstate',[0; pi; 0; 0]);
sys=greyest(data,sys, opt);

opt = compareOptions('InitialCondition',[0;0; 0; 0],'outputOffset',[-0.07;pi]);

compare(data,sys,Inf,opt)
legend('Collected measurement','Greybox model')
%% Observer plot
sys_closed = ss(A - B*K, B, C, D);
L=place(A',C', DP*4)'; %Observer Gain
gain=dcgain(sys_closed);
gain=gain(1);
At = [A   -B*K; L*C  A-(L*C)-(B*K)]; %Augmented State Matrices for Observer
Bt = [B*inv(gain); B*inv(gain)] ;  % 
Ct = [C zeros(2,size(C,2));zeros(2,size(C,2)) C];
sys__ob_closed = ss(At,Bt,Ct,[D;D]); %output feedback close-loop system

x0 = [-0.4*pi pi 0.02 0.02]; % Initial condition
t3 = 1:2*h:4;      % Simulation time 

u=input.Data(400:2*h:403);
[Y, Tsim, X] = lsim(sys__ob_closed,u,t3,[x0 0 0 0 0]);  % Simulate

figure(14)
stairs(Tsim,Y(:,1))  
hold on;
stairs(Tsim,Y(:,2))  
title('Step response for tracking the upright position')
grid on;
legend('State Feedback','State Feedback and Observer')
xlabel('Time(s)')
ylabel('Alpha(rad)')

figure(15)
subplot(2, 2, 1)
%estimation error of true state x1
plot(Tsim,X(:,1)-X(: ,5))
title('Estimation Error of State 1')
xlabel('Time(s)')
ylabel('Error')
grid on;
subplot(2, 2, 2)
%estimation error of true state x2
plot(Tsim,X(:,2)-X(: ,6))
title('Estimation Error of State 2')
xlabel('Time(s)')
ylabel('Error')
grid on;
subplot(2, 2, 3)
%estimation error of true state x3
plot(Tsim,X(:,3)-X(: ,7))
title('Estimation Error of State 3')
grid on;
xlabel('Time(s)')
ylabel('Error')
subplot(2, 2, 4)
%estimation error of true state x4
plot(Tsim,X(:,4)-X(: ,8))
title('Estimation Error of State 4')
grid on;
xlabel('Time(s)')
ylabel('Error')