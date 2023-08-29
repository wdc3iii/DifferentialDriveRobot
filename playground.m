clear; clc; close all;
load("params.mat")

params.Kp = eye(2) * 100;
params.Kd = eye(2) * 15;
P = lyap([zeros(2, 2), eye(2); -params.Kp, -params.Kd]', eye(4));
params.sigma = 5;

%% Test FBL Controller

x_lin = [0; 0; 0; 1; 0];
% u_lin = [0; 0];
u_lin = [0.1; 0.2];
[A, B, C] = linearizeDynamics(x_lin, u_lin, params);
x0 = x_lin + [0.01; -0.01; 0.02; 0.01; 0.001];
T = 3;

[t, x] = ode45(@(t, x) dynamics(x, fblControl(x, linearizationSolution(A, B, C, x_lin, u_lin, t), params), params), linspace(0, T), x0);
[t2, x2] = ode45(@(t, x) A * x + B * u_lin + C, linspace(0, T), x_lin);
[t3, x3] = ode45(@(t, x) dynamics(x, u_lin, params), linspace(0, T), x_lin);
dx2 = (A * x2' + B * u_lin + C)';
figure()
hold on
plotTrajectory(x, t);
plotTrajectory(x2, t2);
plotTrajectory(x3, t3);
legend('Tracking Linearization', 'Linearization Plan', 'Nonlinear Solution')
hold off

figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(t2, x2(:, 1))
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(t2, x2(:, 2))
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 3))
plot(t2, x2(:, 3))
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4))
plot(t2, x2(:, 4))
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(t2, x2(:, 5))
hold off


figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(t2, x2(:, 1))
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(t2, x2(:, 2))
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 4) .* cos(x(:, 3)))
plot(t2, dx2(:, 1))
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4) .* sin(x(:, 3)))
plot(t2, dx2(:, 2))
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(t2, x2(:, 5))
hold off

figure()
plot(t, vecnorm([x(:, 1:2), x(:, 4) .* cos(x(:, 3)), x(:, 4) .* sin(x(:, 3))]' - [x2(:, 1:2), dx2(:, 1:2)]'))


