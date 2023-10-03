%% Setup Simulation
load("params.mat")

% params.Kp = eye(2) * 100;
% params.Kd = eye(2) * 15;
rho = 0.01;
K = lqr([zeros(2), eye(2); zeros(2, 4)], [zeros(2); eye(2)], eye(4), rho * eye(2));
params.Kp = K(1:2, 1:2);
params.Kd = K(1:2, 3:4);
params.P = lyap([zeros(2, 2), eye(2); -params.Kp, -params.Kd]', eye(4));

params.sigma = 10;

lam = min([1 / max(eig(params.P)), params.sigma]);

T = 0.1;
x0 = 0; y0 = 0; vx0 = 0.5; vy0 = 0; ax0 = 0; ay0 = 0;
vx = [0.5, 0.5, 0.5, 0.7, 0.8, 0.9, 1, 0.9, 0.8, 0.7, 0.5, 0.5];
vy = [0, 0.05, 0.1, 0.15, 0.2, 0.2, 0.2, 0.15, 0.1, 0.05, 0];
ax = [0, 0, 0.5, 0.5, 0.5, 0, 0, 0, -0.5, -0.5, -0.5];
ay = [0.25, 0.25, 0.25, 0.25, 0, 0, 0, 0, -0.75, -0.75, -0.75];
jx = [0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, -1, -1, -1, 1, 1, 1];
jy = [1, 1, 1, -1, -1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 0, 0, 0];

input_order = 2; % 1 for vx, vy, 2 for ax,ay, then 3 for jx,jy

delta = [0.01; 0.01; 0; 0; 0];
% delta = [0; 0; 0; 0; 0];
fullX0 = [x0; y0; atan2(vy0, vx0); sqrt(vx0^2 + vy0^2); 0] + delta;

%% Process Simulation

x = squeeze(out.yout{1}.Values.Data);
dx = squeeze(out.yout{2}.Values.Data);
xd = squeeze(out.yout{3}.Values.Data);
dxd = squeeze(out.yout{4}.Values.Data);
ddxd = squeeze(out.yout{5}.Values.Data);
dddxd = squeeze(out.yout{6}.Values.Data);
e = squeeze(out.yout{7}.Values.Data);
de = squeeze(out.yout{8}.Values.Data);
dde = squeeze(out.yout{9}.Values.Data);
dde_bar = squeeze(out.yout{10}.Values.Data);
omega_des = squeeze(out.yout{11}.Values.Data);
domega_des = squeeze(out.yout{12}.Values.Data);
accel = squeeze(out.yout{13}.Values.Data);
alpha = squeeze(out.yout{14}.Values.Data);
V = squeeze(out.yout{15}.Values.Data);
dV_obs = squeeze(out.yout{16}.Values.Data);
dV_des = squeeze(out.yout{17}.Values.Data);
dV_exp = squeeze(out.yout{18}.Values.Data);
t = squeeze(out.yout{19}.Values.Data);
de_obs = squeeze(out.yout{20}.Values.Data);
dde_obs = squeeze(out.yout{21}.Values.Data);
domega_obs = squeeze(out.yout{22}.Values.Data);
u = squeeze(out.yout{23}.Values.Data);

% Visualize
close all;
% Tracking
figure()
subplot(5, 1, 1)
hold on
plot(t, x(1, :))
plot(t, xd(1, :))
hold off
xlabel('Time (s)')
ylabel('X (m)')
legend('Actual', 'Desired')

subplot(5, 1, 2)
hold on
plot(t, x(2, :))
plot(t, xd(2, :))
hold off
xlabel('Time (s)')
ylabel('Y (m)')

subplot(5, 1, 3)
hold on
plot(t, dx(1, :))
plot(t, dxd(1, :))
hold off
xlabel('Time (s)')
ylabel('Xdot (m/s)')

subplot(5, 1, 4)
hold on
plot(t, dx(2, :))
plot(t, dxd(2, :))
hold off
xlabel('Time (s)')
ylabel('Ydot (m/s)')

subplot(5, 1, 5)
hold on
plot(t, x(5, :))
plot(t, omega_des)
hold off
xlabel('Time (s)')
ylabel('omega (rad/s)')

% Tracking Error
figure()
subplot(3, 1, 1)
hold on
plot(t, e(1, :))
plot(t, e(2, :))
hold off
xlabel('Time (s)')
ylabel('Position error (m)')
legend('X', 'Y')

subplot(3, 1, 2)
hold on
plot(t, de(1, :))
plot(t, de(2, :))
plot(t, de_obs(1, :), '--')
plot(t, de_obs(2, :), '--')
hold off
xlabel('Time (s)')
ylabel('Velocity error (m/s)')
legend('X', 'Y', 'X obs', 'Y obs')


subplot(3, 1, 3)
hold on
plot(t, dde(1, :))
plot(t, dde(2, :))
plot(t, dde_obs(1, :), '--')
plot(t, dde_obs(2, :), '--')
% plot(t, dde_bar(1, :))
% plot(t, dde_bar(2, :))
xlabel('Time (s)')
ylabel('Accel error (m/s)')
legend('X', 'Y', 'X obs', 'Y obs')

% Lyapunov
figure()
subplot(2, 1, 1)
hold on
plot(t, V)
plot(t, exp(-lam * t) * V(1))
xlabel('Time(s)')
ylabel('V')
legend('V', 'V bound')

subplot(2, 1, 2)
hold on
plot(t, dV_obs)
plot(t, dV_des, '.-')
plot(t, dV_exp, '--')
hold off
legend('Observed Vdot', 'Desired Vdot under feedback', 'Expected Vdot under feedback')
xlabel('Time (s)')
ylabel('Vdot')

% Omega
figure()
subplot(2, 1, 1)
hold on
plot(t, omega_des)
plot(t, x(5, :))
hold off
legend('omega_des', 'omega')
xlabel('Time (s)')
ylabel('OmegaDes')

subplot(2, 1, 2)
hold on
plot(t, domega_des)
plot(t, domega_obs, '--')
hold off
legend('domega des', 'domega obs')
xlabel('Time (s)')
ylabel('dOmegaDes')

% Inputs
figure()
subplot(2, 1, 1)
hold on
plot(t, accel)
plot(t, alpha)
hold off
legend('accel', 'alpha')
xlabel('Time (s)')
ylabel('Input')

subplot(2, 1, 2)
hold on
plot(t, u(1, :))
plot(t, u(2, :))
hold off
legend('tau1', 'tau2')
xlabel('Time (s)')
ylabel('Input')

%% Bounds?
delta = 1; % V(0) exp(-lambda * t)
M = [eye(2) zeros(2); zeros(2, 4)];
Phalfinv = params.P ^ -0.5;
Q = Phalfinv * M * Phalfinv * delta;
eig_max = eigs(Q, 1);
disp(eig_max)
