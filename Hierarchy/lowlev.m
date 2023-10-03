function xdot = lowlev(t, x, x0m, params)
%LOWLEV Computes the closed loop dynamics under the low level controller
%   Inputs:
%   t:   time
%   x:   current position
%   x0m: initial position for midlevel
%   u0m: input of the mid-level
% Parse inputs
th = x(3);
v = x(4);
w = x(5);
xd = x0m(1:2) + x0m(3:4) * t + 0.5 * x0m(5:6) * t^2;
dxd = x0m(3:4) + x0m(5:6) * t;
ddxd = x0m(5:6);
dddxd = [0; 0];

%% e
e = x(1:2) - xd;

%% de
de = v * [cos(th); sin(th)] - dxd;

%% Eddbar
edd_bar = ddxd - params.Kp * e - params.Kd * de;

%% Accel
acc = [cos(th), sin(th)] * edd_bar;

%% dde
dde = - ddxd + [cos(th), -v * sin(th); sin(th), v * cos(th)] * [acc; w];

%% Omega_des
omega_des = [-1 / v * sin(th), 1 / v * cos(th)] * edd_bar;

%% Omegad_des
omegad_des = [acc / v^2 * sin(th) - w / v * cos(th), -acc / v^2 * cos(th) - w / v * sin(th)] * edd_bar + [-1 / v * sin(th), 1 / v * cos(th)] * (dddxd - params.Kp * de - params.Kd * dde);

%% Alph
alph = -params.sigma / 2 * (w - omega_des) + omegad_des - 2 * [e', de'] * params.P * [0; 0; -v * sin(th); v * cos(th)];

%% Torque
R = params.R;
L = params.L;
alpha = params.m + 2 * params.Iw / R.^2;
beta = params.I + 2 * L^.2 / R.^2 * params.Iw;
mcd = params.mc * params.d;
gamma = 1 / R / alpha;
delta = L / R / beta;

u = 1 / (2 * gamma * delta) * [delta gamma; delta -gamma] * ([acc; alph] - [mcd * w^2 / alpha; -mcd * w * v / beta]);

%% Dynamics
tr = u(1);
tl = u(2);
xdot = [
    v * cos(th);
    v * sin(th);
    w;
    (mcd * w.^2 + (tr + tl) / R) / alpha;
    (-mcd * w * v + (tr - tl) * L / R) / beta
];
end

