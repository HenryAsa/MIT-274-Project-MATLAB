clear 
name = 'geometry'; 

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
% generalized coordinates + derivatives 
syms t th1 th2 dth1 dth2 ddth1 ddth2 g real

% linkage parameters 
syms m1 m2 l1 l2 real
syms l1c l2c real % COM 
syms Ir1 Ir2 real

% control parameters 
syms tau1 tau2 Fx Fy real

% Group them
q   = [th1  ; th2 ];      % generalized coordinates
dq  = [dth1 ; dth2];    % first time derivatives
ddq = [ddth1;ddth2];   % second time derivatives
% q   = [th1  ; th2 ; x; y];      % generalized coordinates
% dq  = [dth1 ; dth2; dx; dy];    % first time derivatives
% ddq = [ddth1;ddth2;ddx; ddy];   % second time derivatives
u   = [tau1 ; tau2];     % controls
F   = [Fx ; Fy];

p   = [m1 m2 Ir1 Ir2 l1 l2 l1c l2c g]';        % parameters

% Generate Vectors and Derivativess
ihat = [1; 0; 0];
jhat = [0; 1; 0];

khat = cross(ihat,jhat);
e1hat =  cos(th1)*ihat + sin(th1)*jhat;
e2hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

% end of linkage positions 
rA = [0 0];
rB = l1 * e1hat;
rC = rB + l2 * e2hat;

% COM positions 
rBcom = l1c * e1hat;
rCcom = rB +  l2c * e2hat;

% derviatives = velocities 
drB = ddt(rB);
drC = ddt(rC);

dr_Bcom = ddt(rBcom);
dr_Ccom = ddt(rCcom);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces
% F2Q = @(F,r) simplify(jacobian(r,q(1:2))'*(F));    % force contributions to generalized forces
% M2Q = @(M,w) simplify(jacobian(w,dq(1:2))'*(M));   % moment contributions to generalized forces

omega1 = dth1;
omega2 = dth1 + dth2;

% Kinetic Energy 
T1 = (1/2)*m1*dot(dr_Bcom, dr_Bcom) + (1/2)*Ir1*dot(omega1, omega1); 
T2 = (1/2)*m2*dot(dr_Ccom, dr_Ccom) + (1/2)*Ir2*dot(omega2, omega2); 

% Potential Energy 
V1 = m1*g*dot(rBcom, -jhat); 
V2 = m2*g*dot(rCcom, -jhat); 

% Sum the Energies 
T = simplify(T1 + T2);
V = V1 + V2;

% torques 
Q_tau1 = M2Q(tau1*khat,omega1*khat);
Q_tau2 = M2Q(tau2*khat,omega2*khat); 
Q_tau2R= M2Q(-tau2*khat,omega1*khat);


Q_tau = Q_tau1+Q_tau2 + Q_tau2R;

Q = Q_tau;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rB(1:2) rC(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;
% eom = ddt(jacobian(L,dq(1:2)).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = simplify(jacobian(eom,ddq));
b = A*ddq - eom;
% 
% A = simplify(jacobian(eom,ddq(1:2)));
% b = A*ddq(1:2) - eom;


% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(V, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);

% Compute hand jacobian
J = jacobian(rC,q);
% J = jacobian(rC,q(1:20);

% Compute wrist jacobian
JW = jacobian(rB,q);
% JW = jacobian(rB,q(1:2));

% Compute ddt( J )
dJ= reshape( ddt(J(:)) , size(J) );

% Write Energy Function and Equations of Motion
z  = [q ; dq];
% z  = [q(1:2) ; dq(1:2)];

rC = rC(1:2);
drC= drC(1:2);
J  = J(1:2,1:2);
dJ = dJ(1:2,1:2);

r_wrist_com = rCcom; 

y_dot = ddt(dot(rB, jhat)); % vertical velocity of B
x_dot = ddt(dot(rB, ihat)); % horizontal velocity of B

drE_q1 = diff(rB, th1);
drE_q2 = diff(rB, th2);
drE_q = [drE_q1 + drE_q2];

matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(rC,'file',['position_hand'],'vars',{z p});
matlabFunction(r_wrist_com,'file',['position_wrist_com'],'vars',{z p});
matlabFunction(rB,'file',['position_wrist'],'vars',{z p});
matlabFunction(drC,'file',['velocity_hand'],'vars',{z p});
matlabFunction(drB,'file',['velocity_wrist'],'vars',{z p});
matlabFunction(y_dot,'file',['velocity_wrist_vert'],'vars',{z p});
matlabFunction(x_dot,'file',['velocity_wrist_hor'],'vars',{z p});
matlabFunction(J ,'file',['jacobian_hand'],'vars',{z p});
matlabFunction(JW ,'file',['jacobian_wrist'],'vars',{z p});
matlabFunction(dJ ,'file',['jacobian_dot_hand'],'vars',{z p});
matlabFunction(dq ,'file',['dq_1'],'vars',{z p});
matlabFunction(drE_q,'file', 'drE_q','vars',{z p});

matlabFunction(Grav_Joint_Sp ,'file', ['Grav_arm'] ,'vars',{z p});
matlabFunction(Corr_Joint_Sp ,'file', ['Corr_arm']     ,'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});


