function [POLAR, TARGET] = dynamics(dt, POLAR, TARGET)
% Propagates dynamics of POLAR satellite using RK4
state=[POLAR.r; POLAR.v; POLAR.q; POLAR.omega; POLAR.wheelspeeds];
k1 = SAT_dynamics(state, POLAR);
k2 = SAT_dynamics(state+k1*dt/2, POLAR);
k3 = SAT_dynamics(state+k2*dt/2, POLAR);
k4 = SAT_dynamics(state+k3*dt, POLAR);
state = state + (k1 + 2*k2 + 2*k3 + k4)*dt/6;

POLAR.r = state(1:3);
POLAR.v = state(4:6);                                           
POLAR.q = state(7:10)/norm(state(7:10))*sign(state(10));
POLAR.omega = state(11:13);
POLAR.wheelspeeds = state(14:end);
POLAR.wheelspeeds = max(POLAR.wheelspeeds, -POLAR.maxw_wheel);
POLAR.wheelspeeds = min(POLAR.wheelspeeds,  POLAR.maxw_wheel);
POLAR.bias = POLAR.bias + POLAR.sig_u*randn(3,1)*dt;


% Propagate target
state=[TARGET.r; TARGET.v; TARGET.q; TARGET.omega];
k1 = TAR_dynamics(state);
k2 = TAR_dynamics(state+k1*dt/2);
k3 = TAR_dynamics(state+k2*dt/2);
k4 = TAR_dynamics(state+k3*dt);
state = state + (k1 + 2*k2 + 2*k3 + k4)*dt/6;

TARGET.r = state(1:3);
TARGET.v = state(4:6);                                           
TARGET.q = state(7:10)/norm(state(7:10))*sign(state(10));
TARGET.omega = state(11:13);
end

function dx = SAT_dynamics(x, SAT)
mu = 398400;

% extract state for convenience
r = x(1:3);
v = x(4:6);
q = x(7:10);
w = x(11:13);
w_w = x(14:end);

% derivative of each state
dr = v;
dv = -mu/norm(r)^3*r;

% *** TO DO Add other perturbing forces, like J2, drag, sun, etc.
% for now will just be random forces
dv = dv + randn(3,1)*1e-9;

A = q2a(q); J = SAT.I;
dq = 1/2*Xi(q)*w;
h = SAT.wheelaxis*diag(SAT.wheel_J)*w_w;
dw = J\(-cross(w, J*w + h) + SAT.L);
% *** TO DO: add disturbance torques
dw = dw + randn(3,1)*1e-5;

dw_w = -SAT.wheeltorques./SAT.wheel_J;

% lump them together
dx = [dr; dv; dq; dw; dw_w];
end

function dx = TAR_dynamics(x)
mu = 398400;

% extract state for convenience
r = x(1:3);
v = x(4:6);
q = x(7:10);
w = x(11:13);

% derivative of each state
dr = v;
dv = -mu/norm(r)^3*r;

% *** TO DO Add other perturbing forces, like J2, drag, sun, etc.
% for now will just be random forces
dv = dv + randn(3,1)*1e-9;

dq = 1/2*Xi(q)*w; J = eye(3);
dw = -J\(cross(w, J*w));
% *** TO DO: add disturbance torques
dw = dw + randn(3,1)*1e-5;

% lump them together
dx = [dr; dv; dq; dw];
end
