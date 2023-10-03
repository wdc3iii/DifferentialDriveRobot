function [x, u] = midlev(x0, xf, A, B, X, U, Xf, Q, Qf, R, N, Nf)
%MIDLEV Solves an MPC problem, returning solution
% Inputs:
%   x0: vector, initial state
%   X:  polyhedron, valid state region
%   U:  polyhedron, valid input region
%   Xf: polyhedron, terminal constraint (assumed Xf subseteq X)
%   Q:  State cost matrix
%   xQ: desired final position
%   R:  Input cost matrix
if nargin < 12
    Nf = N;
end

% State bounds + terminal constraint
Astate = blkdiag(kron(eye(Nf-1), [X.A, zeros(size(X.A, 1), size(R, 2))]), kron(eye(N-Nf), [Xf.A, zeros(size(Xf.A, 1), size(R, 2))]), Xf.A);
bstate = [kron(ones(Nf-1, 1), X.b); kron(ones(N-Nf, 1), Xf.b); Xf.b];

% Input bounds
Ainput = [kron(eye(N-1), [zeros(size(U.A, 1), size(Q, 2)), U.A]), zeros((N-1) * size(U.A, 1), size(Q, 2))];
binput = kron(ones(N-1, 1), U.b);

% Dynamics Constraints
AB = kron(eye(N-1), [A, B]);
AC = kron(eye(N-1), [zeros(size(B)), -eye(length(x0))]);
Adyn = [AB, zeros(size(AB, 1), length(x0))] + ...
    [zeros(size(AB, 1), length(x0)), AC];
bdyn = zeros((N-1) * length(x0), 1);

% Initial Condition
Aic = eye(length(x0), size(Astate, 2));
bic = x0;

% Construct constraints
Ac = [Astate; Ainput; Adyn; -Adyn; Aic; -Aic];
bc = [bstate; binput; bdyn; -bdyn; bic + 1e-8; -bic + 1e-8];

% Construct H
H = blkdiag(kron(eye(N-1), blkdiag(Q, R)), Qf);

% Construct f
f = -[kron(ones(N-1, 1), [Q; zeros(size(R, 1), size(Q, 2))]); Qf] * xf;

% Solve Quadratic Program
out = quadprog(H, f, Ac, bc); % Solves 0.5x'Hx + f'x s.t. Ax <= b

% Extract x, u
x = blkdiag(kron(eye(N-1), [eye(length(x0)), zeros(length(x0), size(R, 2))]), eye(length(x0))) * out;
M = kron(eye(N-1), [zeros(size(R, 2), length(x0)), eye(size(R, 1))]);
u = [M, zeros(size(M, 1), length(x0))] * out;
x = reshape(x, [length(x0), N]);
u = reshape(u, [size(R, 1), N - 1]);
end

