% sim_furuta_sine.m
% Compare nonlinear and linearised Furuta pendulum responses
% to a sinusoidal motor‐voltage input u(t) = Amp·sin(2πf_sin·t)

clear; clc; close all

%% 1) Physical parameters
I_rzz = 5e-4;      % arm inertia about z [kg·m^2]
I_pzz = 2e-4;      % pendulum inertia about hinge [kg·m^2]
m_p   = 0.05;      % pendulum mass [kg]
L_r   = 0.1;       % arm length [m]
L_p   = 0.2;       % pendulum length [m]
mu_th = 1e-3;      % viscous friction on theta [N·m·s/rad]
mu_al = 1e-3;      % viscous friction on alpha [N·m·s/rad]
K     = 0.02;      % motor torque constant [N·m/A]
K_u   = 1;         % amplifier gain
R     = 2;         % motor resistance [Ω]
g     = 9.81;      % gravity [m/s^2]

%% 2) Sinusoidal input definition
t_final = 10;                       % total simulation time [s]
N       = 2000;                     % number of time‐steps
tspan   = linspace(0, t_final, N)'; % column‐vector of time instants

f_sin = 1;      % sine frequency [Hz]
Amp   = 0.0;   % input amplitude [V]

U = Amp * sin(2*pi*f_sin*tspan);    % u(t) = Amp·sin(2πf_sin·t)
%U=zeros(N,1);
%U(50:100) = Amp * 1;
%U=0.1*ones(N,1);


%% 3) Initial conditions near upright (α ≈ π)
theta_eq = 0;
alpha_eq = pi;
delta    = 0.10;                  
x0_nl    = [theta_eq; alpha_eq+delta; 0; 0];  % for nonlinear sim
x0_lin   = [0; delta; 0; 0];                  % deviation for linear

%% 4) Simulate full nonlinear model (ODE45 with u(t))
odefun = @(t,x) furutaEOM_sine(t, x, I_rzz,I_pzz,m_p,L_r,L_p,mu_th,mu_al,K,R,K_u,g, tspan, U);
[ t_nl, X_nl ] = ode45(odefun, tspan, x0_nl);

%% 5) Build linearised state‐space about α = π
J1 = m_p*L_p^2 + 4*I_pzz;
J2 = m_p*L_r^2 +   I_rzz;
Km = K^2 + R*mu_th;
D  = 4*I_rzz*m_p*L_p^2 + 16*I_pzz*m_p*L_r^2 + 16*I_pzz*I_rzz;

A = [ ...
    0, 0, 1, 0; ...
    0, 0, 0, 1; ...
    0, 4*L_p^2*L_r*g*m_p^2/D,   -4*Km*J1/(R*D),   8*L_p*L_r*m_p*mu_al/D; ...
    0, -8*L_p*g*m_p*J2/D,        8*L_p*L_r*m_p*Km/(R*D),  -16*mu_al*J2/D ...
    ];

B = [ ...
    0; ...
    0; ...
    4*K*K_u*J1/(R*D); ...
   -8*K*K_u*L_p*L_r*m_p/(R*D) ...
    ];

sys_lin = ss(A, B, eye(4), zeros(4,1));

%% 6) Simulate linearised response to the sine input
[Y_lin, t_lin] = lsim(sys_lin, U, tspan, x0_lin);

%% 7) Plot input and both responses
figure;

subplot(3,1,1);
plot(tspan, U, 'k', 'LineWidth',1.5);
ylabel('u(t) (V)');
title(sprintf('Sinusoidal Input: %.2f Hz, %.3f V', f_sin, Amp));

subplot(3,1,2);
plot(t_nl,  X_nl(:,1), 'b', ...
     t_lin, Y_lin(:,1), 'r--', 'LineWidth',1.5);
ylabel('\theta (rad)');
legend('nonlinear','linearised','Location','best');
title('Arm Angle \theta');

subplot(3,1,3);
plot(t_nl,        X_nl(:,2),      'b', ...
     t_lin, Y_lin(:,2)+pi, 'r--', 'LineWidth',1.5);
ylabel('\alpha (rad)');
xlabel('Time (s)');
legend('nonlinear','linearised','Location','best');
title('Pendulum Angle \alpha');

%% 8) 3-D animation of the nonlinear Furuta pendulum
figure('Name','Furuta Pendulum Animation','NumberTitle','off')
axis equal
R = L_r+L_p;
xlim([-R R]); ylim([-R R]); zlim([-L_p L_p]);
view(96,20)
grid on; hold on
xlabel('X'); ylabel('Y'); zlabel('Z');

% initialise arm and pendulum line objects
arm_line  = plot3([0 0],[0 0],[0 0],'LineWidth',3,'Color',[0 0.5 0]);
pend_line = plot3([0 0],[0 0],[0 0],'LineWidth',2,'Color',[0.8 0 0]);

fps=20;

for k = 1:5:length(t_nl)
    th = X_nl(k,1);
    al = X_nl(k,2);
    % arm tip
    A = [L_r*cos(th); L_r*sin(th); 0];
    % pendulum bob
    B = A + [ -L_p*sin(al)*sin(th);
               L_p*sin(al)*cos(th);
               L_p*cos(al) ];
    % update graphics
    set(arm_line,  'XData',[0 A(1)], 'YData',[0 A(2)], 'ZData',[0 A(3)]);
    set(pend_line, 'XData',[A(1) B(1)], 'YData',[A(2) B(2)], 'ZData',[A(3) B(3)]);
    title(sprintf('t = %.2f s', t_nl(k)))
    drawnow
    pause(1/fps)
end


%% Local function: nonlinear EOM with sine input u(t)
function dx = furutaEOM_sine(t, x, I_rzz,I_pzz,m_p,L_r,L_p,mu_th,mu_al,K,R,K_u,g, tvec, Uvec)
    % Unpack state
    theta = x(1);  alpha = x(2);
    dth   = x(3);  dal   = x(4);

    % Interpolate u at time t
    u_val = interp1(tvec, Uvec, t, 'linear', 0);
    tau_m = K*K_u/R * u_val;

    % Mass matrix
    M11 = I_rzz + m_p*L_r^2     + (m_p*L_p^2/4)*sin(alpha)^2;
    M12 = - (m_p*L_r*L_p/2)*cos(alpha);
    M22 = I_pzz + m_p*(L_p^2/4);
    M   = [M11, M12; M12, M22];

    % Coriolis / coupling
    C1 =  (m_p*L_p^2/2)*dal*dth*sin(alpha)*cos(alpha) ...
        + (m_p*L_r*L_p/2)*sin(alpha)*dal^2;
    C2 = -(m_p*L_p^2/4)*dth^2*sin(alpha)*cos(alpha);

    % Friction
    F1 = mu_th * dth;
    F2 = mu_al * dal;

    % Gravity (correct sign)
    G2 = 0.5 * m_p * g * L_p * sin(alpha);

    % RHS
    R1 = tau_m - F1 - C1;
    R2 =   C2 + G2 - F2;

    % Solve for accelerations
    acc = M \ [R1; R2];

    % Pack derivatives
    dx = [ dth; dal; acc(1); acc(2) ];
end
