close all; clear; clc;

syms b1 b2 m1 m2 l1 l2 L1 L2 Jhat1 Jhat2 Jhat0 Km Rm KT g J1 J2
syms alpha alphadot alphaddot theta thetadot thetaddot tau u x

Jhat1 = J1+m1*l1^2;
Jhat2 = J2+m2*l2^2;

Jhat0 = Jhat1 + m2*L1^2;

x = [alpha;theta;alphadot;thetadot];

D = Jhat0*Jhat2 + Jhat2^2 * sin(theta)^2 - m2^2*L1^2*l2^2*cos(theta)^2;

f1 = alphadot;
f2 = thetadot;
f3 = (-Jhat2*b1/D)*alphadot + (m2*L1*l2*cos(theta)*b2/D)*thetadot - (Jhat2^2*sin(2*theta)/D)*alphadot*thetadot -(Jhat2*m2*L1*l2*cos(theta)*sin(2*theta)/(2*D))*alphadot^2 + (Jhat2*m2*L1*l2*sin(theta)/D)*alpha + (m2^2*l2^2*L1*sin(2*theta)/(2*D)) * g;
f4 = (m2*L1*l2*cos(theta)*b1/D)*alphadot - (b2*(Jhat0+Jhat2*sin(theta)^2)/D)*thetadot + (m2*L1*l2*Jhat2*cos(theta)*sin(2*theta)/D)*alphadot*thetadot - (sin(2*theta)*(Jhat0*Jhat2+Jhat2^2*sin(theta)^2)/(2*D))*alphadot^2 - ((m2^2*L1^2*l2^2*sin(2*theta))/(2*D))*thetadot^2 - (m2*l2*sin(theta)*(Jhat0+Jhat2*sin(theta)^2)/D)*g;

fx = [f1;f2;f3;f4];

g1 = 0;
g2 = 0;
g3 = Jhat2/D;
g4 = -(m2*L1*l2*cos(theta)/D);

gx = [g1;g2;g3;g4];

hx = [1,0,0,0;0,1,0,0] * x;

f(alpha,theta,alphadot,thetadot) = fx;

% TODO change tau to -Km^2/Rm * alphadot + 5Km/Rm * u

%% input-output linearization

% Define performance output h 
n3 = g3*D;
n4 = g4*D;
h = theta;

% First Lie derivative
Lgh = jacobian(h, [alpha, theta, alphadot, thetadot]) * gx

% Second Lie derivatice
Lfh = jacobian(h, [alpha, theta, alphadot, thetadot]) * fx;
LgLfh = simplify(jacobian(Lfh, [alpha,theta,alphadot,thetadot]) * gx)

% relative degree is n/2 because of structure! also because of two outputs?

% building the diffeomorphism

psi_r = [hx;Lfh];

phi = [alpha;theta]; % again using the structure

etha = phi;
xi = psi_r;

T = [etha;xi]

%% Checking feedback linearizability using Frobenius Theorem

G1 = gx;
G2 = (jacobian(gx, [alpha,theta,alphadot,thetadot])*fx - jacobian(fx, [alpha,theta,alphadot,thetadot])*gx);
G3 = (jacobian(G2, [alpha,theta,alphadot,thetadot])*fx - jacobian(fx, [alpha,theta,alphadot,thetadot])*G2);
G4 = (jacobian(G3, [alpha,theta,alphadot,thetadot])*fx - jacobian(fx, [alpha,theta,alphadot,thetadot])*G3);

G = simplify([G1,G2,G3,G4])
%% check rank condition

rank(subs(G, [alpha, theta, alphadot, thetadot], [pi, pi, 0, 0])) % probably full rank for large D?
%% input-state linearization

% define a performance output that makes r=n
C1 = g4/g3/(cos(theta));
h = log(abs(sec(theta)+tan(theta))-C1*alpha);

g_f = gx;

% First Lie derivative
Lgh = jacobian(h, [alpha, theta, alphadot, thetadot]) * g_f

% Second Lie derivatice
Lfh = jacobian(h, [alpha, theta, alphadot, thetadot]) * fx
LgLfh = simplify(jacobian(Lfh, [alpha,theta,alphadot,thetadot]) * g_f)
test = simplify(jacobian(h,alpha) + jacobian(h,theta) * C1)

% Third Lie derivative



% building the diffeomorphism

psi_r = [h;Lfh];

phi = [alpha;theta]; % again using the structure

etha = phi;
xi = psi_r;

T = [etha;xi]
%% Finding and using a control Lyapunov function

Ep = g*m2*l2*(1-cos(theta)) + KT*alpha;
Ek = 0.5*alphadot^2*(m1*l1^2+J1) + 0.5*alphadot^2*(m2*L2^2+(m2*l2^2+J2)*sin(theta)^2 +J2*cos(theta)^2)+0.5*thetadot^2*(J2+m2*l2^2)+m2*L1*l2*cos(theta)*alphadot*thetadot;
H = Ep+Ek;
V = simplify(H - +m2*L1*l2*cos(theta)*alphadot*thetadot)
V = 1-cos(theta) + alpha^2 + 1- cos(alphadot) + 1 - cos(thetadot);

% Check if V is a Lyapunov (candidate) function

V0 = subs(V, [alpha, theta, alphadot, thetadot], [0, 0, 0, 0]) %Check: V(0) = 0

V % Check: V(x) >= 0, x in D/{0}

% todo : barbashin krasovski for alpha = 0

% Check if this could be a control Lyapunov function

DelV_f = jacobian(V, [alpha, theta, alphadot, thetadot]) * fx

DelV_g = jacobian(V, [alpha,theta,alphadot,thetadot])*gx

LHS = simplify(DelV_f + DelV_g * u) % becomes too complex
