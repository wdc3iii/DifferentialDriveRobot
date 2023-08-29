function xdot = dynamics(x, u, params)
%DYNAMICS compute the time rate of change of state given the input,
%assuming no constraint violation
%   Detailed explanation goes here
th = x(3);
v = x(4);
w = x(5);
R = params.R;
L = params.L;
alpha = params.m + 2 * params.Iw / R.^2;
beta = params.I + 2 * L^.2 / R.^2 * params.Iw;
mcd = params.mc * params.d;

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

