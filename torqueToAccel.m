function a = torqueToAccel(x, u, params)
%ACCELTOTORQUE Summary of this function goes here
%   Detailed explanation goes here
v = x(4);
w = x(5);
R = params.R;
L = params.L;
alpha = params.m + 2 * params.Iw / R.^2;
beta = params.I + 2 * L^.2 / R.^2 * params.Iw;
mcd = params.mc * params.d;

gamma = 1 / R / alpha;
delta = L / R / beta;
a  = [mcd * w^2 / alpha; -mcd * w * v / beta] + [gamma gamma; delta -delta] * u;
end


