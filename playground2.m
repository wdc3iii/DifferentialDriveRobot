clear; clc; close all;
load("params.mat")

params.Kp = eye(2) * 100;
params.Kd = eye(2) * 15;
params.P = lyap([zeros(2, 2), eye(2); -params.Kp, -params.Kd]', eye(4));
params.sigma = 1;

%% Tracking
T = 0.1;
x0 = 0; y0 = 0; vx0 = 0.5; vy0 = 0; 
ax = 0.1; ay = 0.5;

fullX0 = [x0; y0; atan2(vy0, vx0); sqrt(vx0^2 + vy0^2); 0];

traj_t = @(t) [x0 + vx0 * t + ax / 2 * t .^2; y0 + vy0 * t + ay / 2 * t .^2];
dtraj_t = @(t) [vx0 + ax * t; vy0 + ay * t];
ddtraj_t = @(t) [ax; ay] * ones(size(t));
dddtraj_t = @(t) [0; 0];

[t, x] = ode45(@(t, x) dynamics(x, fblControl(x, [traj_t(t), dtraj_t(t), ddtraj_t(t), dddtraj_t(t)], params), params), linspace(0, T), fullX0);
xd = traj_t(t')';
dxd = dtraj_t(t')';
ddxd = ddtraj_t(t')';

figure()
hold on
plotTrajectory(xd, t);
plotTrajectory(x, t);
legend('RoM', 'FoM Tracking')
hold off

figure()
subplot(4, 1, 1)
hold on
plot(t, xd(:, 1))
plot(t, x(:, 1))
xlabel('x')
hold off

subplot(4, 1, 2)
hold on
plot(t, xd(:, 2))
plot(t, x(:, 2))
xlabel('x')
hold off

subplot(4, 1, 3)
hold on
plot(t, x(:, 3))
plot(t, atan2(dxd(:, 2), dxd(:, 1)))
xlabel('theta')
hold off

subplot(4, 1, 4)
hold on
plot(t, x(:, 4))
plot(t, vecnorm(dxd')')
xlabel('v')
hold off

% Compute V
e = x(:, 1:2) - xd;
ed = x(:, 4) .* [cos(x(:, 3)), sin(x(:, 3))] - dxd;

Ve = diag([e, ed] * params.P * [e, ed]');

edd_bar = ddxd' - params.Kp * e' - params.Kd * ed';
LgLfinv = 1 ./ x(:, 4) .* [-sin(x(:, 3)), cos(x(:, 3))];

accel = diag([cos(x(:,3)), sin(x(:, 3))] * edd_bar);
omega_des = diag(LgLfinv * edd_bar);

V = Ve + 0.5 * (omega_des - x(:, 5)).^2;

dVe = -diag([e, ed] * [e, ed]');
dV = dVe - params.sigma * (omega_des - x(:, 5)).^2;


figure()
subplot(2, 1, 1)
hold on
plot(t, V)
plot(t, Ve)
plot(t, 0.5 * (omega_des - x(:, 5)).^2)
xlabel('Time')
legend('V', 'Ve', 'omega')

subplot(2, 1, 2)
hold on
plot(t, dV)
plot(t, dVe)
plot(t, -params.sigma * (omega_des - x(:, 5)).^2)
plot(t(2:end), diff(V) ./ diff(t))
plot(t(2:end), diff(Ve) ./ diff(t))
plot(t(2:end), diff(0.5 * (omega_des - x(:, 5)).^2) ./ diff(t))
legend('expected dV', 'expected dVe', 'expected domega', 'dV', 'dVe', 'domega')
xlabel('Time')
