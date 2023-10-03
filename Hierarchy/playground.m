%% Work on predictive control stuff
clear; clc; close all;

T = 0.1;    % Time step 0.1s
N = 25;     % Horizon length
vmax = 1;   % maximum velocity, m/s
amax = 5;   % maximum acceleration, m/s^2
jmax = 1;   % maximum jerk (delta input), m/s^3

load("../inv_set.mat")

%% Define Dynamics

A = [
    1, 0, T, 0, 0.5*T^2, 0;
    0, 1, 0, T, 0, 0.5*T^2;
    0, 0, 1, 0, T, 0;
    0, 0, 0, 1, 0, T;
    zeros(2, 4), eye(2)
];
B = [zeros(4, 2); eye(2)];

% Hv x <= hv
Hx = [[eye(2); -eye(2)], zeros(4, 4)];
hx = [2; 1; 0; 0];

Hv = [zeros(4, 2), [eye(2); -eye(2)], zeros(4, 2)];
hv = [vmax; vmax; vmax; vmax];

Ha = [zeros(4, 2), zeros(4, 2), [eye(2); -eye(2)]];
ha = [amax; amax; amax; amax];

Hj = [eye(2); -eye(2)];
hj = [jmax; jmax; jmax; jmax];

X = Polyhedron('A', [Hx; Hv; Ha], 'b', [hx; hv; ha]);
U = Polyhedron('A', Hj, 'b', hj);
Xf =  Bctrlinv;

Q = diag([100, 100, 5, 5, 1, 1]);
R = .1 * eye(2);
Qf = diag([100, 100, 5, 5, 1, 1]);

x0 = [0.5; 0.25; 0; 0; 0; 0];
xf = [1.5; 0.75; 0; 0; 0; 0];
% x0 = [0.930000000000014;
% 0.560608666873694;
% 1;
% 0.544447626152260;
% 8.88178419700125e-16;
% -0.869135534234132];
% xf = [1.5; 0.75; 0; 0; 0; 0];

[x, u] = midlev(x0, xf, A, B, X, U, Xf, Q, Qf, R, N);

figure()
hold on
% axis("square")
Polyhedron('A', Hx(:, 1:2), 'b', hx).plot('wire', true, 'linewidth', 2)
Polyhedron('A', Hx(:, 1:2), 'b', [1,1,0,0]).plot('wire', true, 'linewidth', 2)
plot(x(1, :), x(2, :), '*-')


figure()
subplot(4, 1, 1)
t = linspace(0, N* T, N);
plot(t, x(1:2, :))
legend('x', 'y')

subplot(4, 1, 2)
plot(t, x(3:4, :))
legend('$\dot{x}$', '$\dot{y}$')

subplot(4, 1, 3)
plot(t, x(5:6, :))
legend('$u_x$', '$u_y$')

subplot(4, 1, 4)
plot(t(1:end-1), u)
legend('$\Delta u_x$', '$\Delta u_y$')

