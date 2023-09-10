clear; clc; close all;
load("params.mat")

params.Kp = eye(2) * 100;
params.Kd = eye(2) * 15;
params.P = lyap([zeros(2, 2), eye(2); -params.Kp, -params.Kd]', eye(4));
params.sigma = 1;

%% Test FBL Controller 1
T = 14;
xl = [0; 0; -1; 1; 0];
ul = [0.1; 0.2];
[A, B, C] = linearizeDynamics(xl, ul, params);

x0 = xl + [0.01; -0.01; 0.02; 0.01; 0.001];

[t, x] = ode45(@(t, x) dynamics(x, fblControl(x, linearizationSolution(A, B, C, xl, ul, t), params), params), linspace(0, T), x0);
[tlinsol, xlinsol] = ode45(@(t, x) A * x + B * ul + C, linspace(0, T), xl);
[tdynsol, xdynsol] = ode45(@(t, x) dynamics(x, ul, params), linspace(0, T), x0);
dxlinsol = (A * xlinsol' + B * ul + C)';

figure()
hold on
plotTrajectory(x, t);
plotTrajectory(xlinsol, tlinsol);
plotTrajectory(xdynsol, tdynsol);
legend('Tracking Linearization', 'Linearization Plan', 'Nonlinear Solution')
hold off

figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(tlinsol, xlinsol(:, 1))
xlabel('x')
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(tlinsol, xlinsol(:, 2))
xlabel('y')
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 3))
plot(tlinsol, xlinsol(:, 3))
xlabel('theta')
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4))
plot(tlinsol, xlinsol(:, 4))
xlabel('v')
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(tlinsol, xlinsol(:, 5))
xlabel('omega')
hold off


figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(tlinsol, xlinsol(:, 1))
xlabel('x')
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(tlinsol, xlinsol(:, 2))
xlabel('y')
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 4) .* cos(x(:, 3)))
plot(tlinsol, dxlinsol(:, 1))
xlabel('xdot')
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4) .* sin(x(:, 3)))
plot(tlinsol, dxlinsol(:, 2))
xlabel('ydot')
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(tlinsol, xlinsol(:, 5))
xlabel('omega')
hold off

figure()
hold on
norm_err = vecnorm([x(:, 1:2), x(:, 4) .* cos(x(:, 3)), x(:, 4) .* sin(x(:, 3))]' - [xlinsol(:, 1:2), dxlinsol(:, 1:2)]');
plot(t, norm_err)
eg = eig(params.P);
k1 = min([min(eg), 0.5]);
k2 = max([max(eg), 0.5]);
k3 = min([1/max(eg), 2 * params.sigma]);

plot(t, k2 / k1 * exp(-k3 * t) * norm_err(1))

V = zeros(size(x, 1), 1);
dV = zeros(size(x, 1), 1);
for ii = 1:size(x, 1)
    V(ii) = Veta(x(ii, :)', linearizationSolution(A, B, C, xl, ul, t(ii)), params);
    dV(ii) = dotV(x(ii, :)', linearizationSolution(A, B, C, xl, ul, t(ii)), params);
end
figure()
hold on
plot(V)
plot(dV)
plot(diff(V) ./ diff(t))

%% Test FBL Controller
T = 0.2;
x1 = [0; 0; 0; 1; 0];
u1 = [0.1; 0.2];
[A1, B1, C1] = linearizeDynamics(x1, u1, params);
x2 = linearizationSolution(A1, B1, C1, x1, u1, T);
x2 = x2(:, 1);
u2 = [0.2; 0.1];
[A2, B2, C2] = linearizeDynamics(x2, u2, params);
x3 = linearizationSolution(A2, B2, C2, x2, u2, T);
x3 = x3(:, 1);
u3 = [1; 1];
[A3, B3, C3] = linearizeDynamics(x3, u3, params);
x4 = linearizationSolution(A3, B3, C3, x3, u3, T);
x4 = x4(:, 1);
u4 = [0; 0];
[A4, B4, C4] = linearizeDynamics(x4, u4, params);

A = {A1, A2, A3, A4};
B = {B1, B2, B3, B4};
C = {C1, C2, C3, C4};
xlin = {x1, x2, x3, x4};
ulin = {u1, u2, u3, u4};

x0 = x1; %+ [0.01; -0.01; 0.02; 0.01; 0.001];

[t, x] = ode45(@(t, x) dynamics(x, fblControl(x, linearizationSolution(A{floor(t / T) + 1}, B{floor(t / T) + 1}, C{floor(t / T) + 1}, xlin{floor(t / T) + 1}, ulin{floor(t / T) + 1}, t), params), params), linspace(0, 4*T - eps), x0);
[tlinsol, xlinsol] = ode45(@(t, x) A{floor(t / T) + 1} * x + B{floor(t / T) + 1} * ulin{floor(t / T) + 1} + C{floor(t / T) + 1}, linspace(0, 4 * T - eps), x1);
[tdynsol, xdynsol] = ode45(@(t, x) dynamics(x, ulin{floor(t / T) + 1}, params), linspace(0, 4 * T - eps), x0);

figure()
hold on
plotTrajectory(x, t);
plotTrajectory(xlinsol, tlinsol);
plotTrajectory(xdynsol, tdynsol);
legend('Tracking Linearization', 'Linearization Plan', 'Nonlinear Solution')
hold off

figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(tlinsol, xlinsol(:, 1))
xlabel('x')
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(tlinsol, xlinsol(:, 2))
xlabel('y')
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 3))
plot(tlinsol, xlinsol(:, 3))
xlabel('theta')
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4))
plot(tlinsol, xlinsol(:, 4))
xlabel('v')
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(tlinsol, xlinsol(:, 5))
xlabel('omega')
hold off


figure()
subplot(5, 1, 1)
hold on
plot(t, x(:, 1))
plot(tlinsol, xlinsol(:, 1))
xlabel('x')
hold off

subplot(5, 1, 2)
hold on
plot(t, x(:, 2))
plot(tlinsol, xlinsol(:, 2))
xlabel('y')
hold off

subplot(5, 1, 3)
hold on
plot(t, x(:, 4) .* cos(x(:, 3)))
% plot(tlinsol, dxlinsol(:, 1))
xlabel('xdot')
hold off

subplot(5, 1, 4)
hold on
plot(t, x(:, 4) .* sin(x(:, 3)))
% plot(tlinsol, dxlinsol(:, 2))
xlabel('ydot')
hold off

subplot(5, 1, 5)
hold on
plot(t, x(:, 5))
plot(tlinsol, xlinsol(:, 5))
xlabel('omega')
hold off

figure()
plot(t, vecnorm([x(:, 1:2), x(:, 4) .* cos(x(:, 3)), x(:, 4) .* sin(x(:, 3))]' - [xlinsol(:, 1:2), dx2(:, 1:2)]'))

V = zeros(size(x, 1), 1);
dV = zeros(size(x, 1), 1);
for ii = 1:size(x, 1)
    V(ii) = Veta(x(ii, :)', linearizationSolution(A, B, C, x_lin, u_lin, t(ii)), params);
    dV(ii) = dotV(x(ii, :)', linearizationSolution(A, B, C, x_lin, u_lin, t(ii)), params);
end
figure()
hold on
plot(V)
plot(dV)
plot(diff(V) ./ diff(t))

%% Troublehshooting
fblControl(x(end, :)', linearizationSolution(A, B, C, x_lin, u_lin, T), params);


function V = Veta(x, traj, params)
ybar = traj(1:2, 1);
ydbar = traj(1:2, 2);
yddbar = traj(1:2, 3);

e = x(1:2) - ybar;
ed = [x(4) * cos(x(3)); x(4) * sin(x(3))] - ydbar;

var1 = yddbar - params.Kp * e - params.Kd * ed;
var2 = [-1 / x(4) * sin(x(3)) 1 / x(4) * cos(x(3))];
omega_des = var2 * var1;

V = [e; ed]' * params.P * [e; ed] + 0.5 * (x(5) - omega_des)^2;
end

function dV = dotV(x, traj, params)
ybar = traj(1:2, 1);
ydbar = traj(1:2, 2);
yddbar = traj(1:2, 3);

e = x(1:2) - ybar;
ed = [x(4) * cos(x(3)); x(4) * sin(x(3))] - ydbar;

var1 = yddbar - params.Kp * e - params.Kd * ed;
var2 = [-1 / x(4) * sin(x(3)) 1 / x(4) * cos(x(3))];
omega_des = var2 * var1;

dV = norm([e; ed])^2 - params.sigma * (x(5) - omega_des)^2;
end