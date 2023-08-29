function [A, B, C] = linearizeDynamics(x0, u0, params)
%LINEARIZE Summary of this function goes here
%   Detailed explanation goes here

th = x0(3);
v = x0(4);
w = x0(5);
R = params.R;
L = params.L;
alpha = params.m + 2 * params.Iw / R.^2;
beta = params.I + 2 * L^.2 / R.^2 * params.Iw;
mcd = params.mc * params.d;

A = [
    0, 0, -v * sin(th), cos(th), 0;
    0, 0, v * cos(th), sin(th), 0;
    0, 0, 0, 0, 1;
    0, 0, 0, 0, 2 * mcd / alpha * w;
    0, 0, 0, -mcd * w / beta, mcd * v / beta
];
B = [
    0 0;
    0 0;
    0 0;
    1 / (R * alpha), 1 / (R * alpha);
    L / (R * beta), -L / (R * beta)
];
C = dynamics(x0, u0, params) - A * x0 - B*u0;

end

