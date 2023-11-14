clear 
name = 'ball'; 

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 

% generalized coordinates + derivatives + mass
syms t x y dx dy ddx ddy g real
syms m r real

% control parameters
syms Fx Fy real

% Group them
q = [x; y];
dq = [dx; dy];
ddq = [ddx; ddy];

F = [Fx; Fy];

p = [m g r]';

%Generate Vectors and Derivatives
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

r = [x, y, 0];
ddt = @(r) jacobian(r,[q;dq])*[dq;ddq];

vr = ddt(r); 

% Calculate Kinetic Energy, Potential Energy
% T = (1/2)*m*(dot(dx, dx)+dot(dy, dy));
T = (1/2)*m*dot(vr, vr);
V = m*g*y;

%% Derive Energy Function & EOM
E = T+V;
L = T-V;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).';

% Rearrange EOM
A = simplify(jacobian(eom,ddq));
b = A*ddq - eom;

% EOM
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(V, q)');
Corr_Joint_Sp = simplify( eom - Grav_Joint_Sp - A*ddq);

z = [q; dq];

% Compute ball jacobian
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(q,'file',['position_ball'],'vars',{z p});
matlabFunction(dq, 'file', ['velocity_ball'], 'vars',{z p});






