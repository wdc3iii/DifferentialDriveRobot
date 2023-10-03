%% Define An Interesting Problem...
clear; clc; close all;
load("../params.mat");

%% High-level constraints
% 10 x 10 grid, 1m square boxes
% Obstacles at (2, 9), (3, 5), (3, 7), (7, 6), (8, 3), (9, 2), (10, 7)

% From verification, we know that adjacent 1x1 boxes are reachable.
% Adjacency matrix is therefore
dis = 100;
adj = eye(dis);
adj = adj + diag(ones(1,dis-1),1) + diag(ones(1,dis-1), -1); % Encodes all columnwise shifts
adj = adj + diag(ones(1,dis-10), 10) + diag(ones(1,dis-10), -10); % Encodes all rowwise shifts

obst = [2, 3, 3, 7, 8, 9, 10; 9, 5, 7, 6, 3, 2, 7];

for ii = 1:size(obst, 2)
    r = obst(2, ii);
    c = obst(1, ii);
    loc = c + 10 * (r - 1);
    adj(loc, loc) = 0;
    if r < 10
        adj(loc + 10, loc) = 0;
    end
    if r > 1
        adj(loc - 10, loc) = 0;
    end
    if c < 10
        adj(loc, loc + 1) = 0;
    end
    if c > 1
        adj(loc, loc - 1) = 0;
    end
end

%% Mid-level Constraints
load("../inv_set.mat")
T = 0.1;    % Time step 0.1s
N = 30;     % Horizon length
vmax = 1;   % maximum velocity, m/s
amax = 1;   % maximum acceleration, m/s^2
jmax = 1;   % maximum jerk (delta input), m/s^3

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

%% Mid-level Dynamics
A = [
    1, 0, T, 0, 0.5*T^2, 0;
    0, 1, 0, T, 0, 0.5*T^2;
    0, 0, 1, 0, T, 0;
    0, 0, 0, 1, 0, T;
    zeros(2, 6)
    ];
B = [zeros(4, 2); eye(2)];

%% Cost matrices
Q = eye(6);
R = eye(2);
Qf = eye(6);

%% Control Invariant Sets
ctrlInvs = cell(10, 10);
for r = 1:10
    for c = 1:10
        off = [c-1; r - 1; 0; 0; 0; 0];
        ctrlInvs{c, r} = Polyhedron('A', Actrlinv.A, 'b', Actrlinv.b + Actrlinv.A * off);
    end
end

%% Low Level Control
rho = 0.01;
K = lqr([zeros(2), eye(2); zeros(2, 4)], [zeros(2); eye(2)], eye(4), rho * eye(2));
params.Kp = K(1:2, 1:2);
params.Kd = K(1:2, 3:4);
params.P = lyap([zeros(2, 2), eye(2); -params.Kp, -params.Kd]', eye(4));

params.sigma = 10;
