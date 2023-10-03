%% Work on predictive control stuff
clear; clc; close all;

T = 0.1;    % Time step 0.1s
vmax = 1;   % maximum velocity, m/s
amax = 5;   % maximum acceleration, m/s^2
jmax = 1;  % maximum jerk (delta input), m/s^3

height = 3; % environment size
width = 3;  % environment size

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
Hv = [zeros(4, 2), [eye(2); -eye(2)], zeros(4, 2)];
hv = [vmax; vmax; vmax; vmax];

Ha = [zeros(4, 2), zeros(4, 2), [eye(2); -eye(2)]];
ha = [amax; amax; amax; amax];


Hj = [eye(2); -eye(2)];
hj = [jmax; jmax; jmax; jmax];

%% Define a decomposition of space
% y^ 
%  | G, H, I
%  | D, E, F
%  | A, B, C
%  ---------> x

% A1 x <= a1

box = [[eye(2); -eye(2)], zeros(4, 4)];
box_lim = [height; width; 0; 0];

A1 = [[eye(2); -eye(2)], zeros(4, 4)];
a1 = [1; 1; 0; 0];
B1 = [[eye(2); -eye(2)], zeros(4, 4)];
b1 = [2; 1; -1; 0];
C1 = [[eye(2); -eye(2)], zeros(4, 4)];
c1 = [3; 1; -2; 0];
D1 = [[eye(2); -eye(2)], zeros(4, 4)];
d1 = [1; 2; 0; -1];
E1 = [[eye(2); -eye(2)], zeros(4, 4)];
e1 = [2; 2; -1; -1];
F1 = [[eye(2); -eye(2)], zeros(4, 4)];
f1 = [3; 2; -2; -1];
G1 = [[eye(2); -eye(2)], zeros(4, 4)];
g1 = [1; 3; 0; -2];
H1 = [[eye(2); -eye(2)], zeros(4, 4)];
h1 = [2; 3; -1; -2];
I1 = [[eye(2); -eye(2)], zeros(4, 4)];
i1 = [3; 3; -2; -2];

A1uB1 = [[eye(2); -eye(2)], zeros(4, 4)];
a1ub1 = [2; 1; 0; 0];

%% Now, lets compute the maximal control invariant set inside a set
% For instance, E
X = [A1; Hv; Ha];
x = [a1; hv; ha];


Omega = Polyhedron('A', X, 'b', x);
Omega.minHRep();
Hprev = zeros(size(Omega.H));
k = 0;
tol = 1e-4;
while size(Hprev, 1) ~= size(Omega.H, 1) || norm(sortrows(Hprev) - sortrows(Omega.H), "fro") > tol
    H = Omega.H;
    Hx = H(:, 1:end-1);
    hx = H(:, end);

    % Compute Pre(Omega_k)
    extendedOmega = Polyhedron('A', [Hx * A, Hx * B; zeros(4,6), Hj], 'b', [hx; hj]);
    Omega_k1 = extendedOmega.projection(1:6);
    % Omega_k+1 = Pre(omega_k) intersected Omega_k
    Omega = Polyhedron('H', [Omega.H; Omega_k1.H]);
    Omega.minHRep();
    k = k + 1;
    Hprev = H;
    fprintf("Iteration %d\n", k)
end

Actrlinv = Omega;
b_off = [1; 0; 0; 0; 0; 0];
Bctrlinv = Polyhedron('A', Omega.A, 'b', Omega.b + Omega.A * b_off);
save('inv_set.mat', 'Actrlinv', 'Bctrlinv');

% This algorithm appears to be identical
sys = LTISystem('A', A, 'B', B);
sys.x.min = [0; 0; -1; -1; -5; -5];
sys.x.max = [1; 1; 1; 1; 5; 5];
sys.u.min = [-1; -1];
sys.u.max = [1; 1];

InvSet = sys.invariantSet();

%% Check how big inv set is if projected onto position
figure()
pos = Omega.projection([1:2, 3]);
pos.plot()
% Note no volume is lost

%% Compute Vertex Representation of A
Actrlinv.computeVRep();

%% Now, see if B is reachable from all points in A
X = [A1uB1; Hv; Ha];
x = [a1ub1; hv; ha];

X = Polyhedron('A', X, 'b', x);
X.minHRep();

Kappa = Bctrlinv;
maxHorizon = 30;
for ii = 1:maxHorizon
    H = Kappa.H;
    Hx = H(:, 1:end-1);
    hx = H(:, end);

    % Compute Pre(Omega_k)
    extendedKappa = Polyhedron('A', [Hx * A, Hx * B; zeros(4,6), Hj], 'b', [hx; hj]);
    Kappa_k1 = extendedKappa.projection(1:6);
    % Kappa_k+1 = Pre(Kappa_k) intersected X
    Kappa = Polyhedron('H', [Kappa_k1.H; X.H]);
    Kappa.minHRep();
    contains_extreme_points = Kappa.A * Actrlinv.V' <= Kappa.b + tol;
    reachable = all(contains_extreme_points, 1);
    fprintf("Iteration %d, Reachable: %0.4f\n", ii, sum(reachable) / length(reachable))
end

%% Check how big inv set is if projected onto position
figure()
pos = Kappa.projection([1:2, 3]);
pos.plot()
% Note no volume is lost

save('inv_set.mat', 'Actrlinv', 'Bctrlinv');