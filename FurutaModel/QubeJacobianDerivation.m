% Define symbols
syms I_rzz I_pzz m_p L_r L_p mu_theta mu_alpha K R K_u g real
syms theta alpha dtheta dalpha u real

% Define second-derivative symbols
syms ddtheta ddalpha real

% EOM (left-hand sides minus right-hand sides = 0)
eq1 = (I_rzz + m_p*L_r^2 + (m_p*L_p^2/4)*sin(alpha)^2)*ddtheta ...
    + (m_p*L_p^2/2)*dot(alpha,dtheta)*sin(alpha)*cos(alpha) ...
    - (m_p/2)*L_r*L_p*( cos(alpha)*ddalpha - sin(alpha)*dalpha^2 ) ...
    + (mu_theta + K^2/R)*dtheta ...
    - (K*K_u/R)*u;

eq2 = (I_pzz + m_p*L_p^2/4)*ddalpha ...
    - (m_p/2)*L_r*L_p*ddtheta*cos(alpha) ...
    - (m_p*L_p^2/4)*dtheta^2*sin(alpha)*cos(alpha) ...
    - (m_p/2)*g*L_p*sin(alpha) ...
    + mu_alpha*dalpha;

% Solve for accelerations
sol = solve([eq1; eq2], [ddtheta; ddalpha]);
ddtheta = simplify(sol.ddtheta);
ddalpha = simplify(sol.ddalpha);

% Form state-vector and vector field
x = [theta; alpha; dtheta; dalpha];
f = [ dtheta;
      dalpha;
      ddtheta;
      ddalpha ];

% Compute Jacobian w.r.t. [theta; alpha; dtheta; dalpha; u]
J = jacobian(f, [theta; alpha; dtheta; dalpha; u]);

% Display results
disp('State-space vector f = '), disp(f)
disp('Jacobian J = '), disp(J)

J_stable = subs(J,{alpha,dalpha,dtheta}, {pi,0,0});
A_stable = J_stable(:,1:end-1)
B_stable = J_stable(:,end)

J_unstable = subs(J,{alpha,dalpha,dtheta}, {0,0,0});
A_unstable = J_unstable(:,1:end-1)
B_unstable = J_unstable(:,end)