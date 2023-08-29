function traj = linearizationSolution(A, B, C, x0, u, t)
%LINEARATIONSOLUTION Summary of this function goes here
%   Detailed explanation goes here
if t == 0
    x_t = x0;
else
    [~, x_t] = ode45(@(t, x) A * x + B * u + C, [0, t], x0);
    x_t = x_t(end, :)';
end

dx_t = A * x_t + B * u + C;
ddx_t = A * dx_t;
dddx_t = A * ddx_t;
traj = [x_t, dx_t, ddx_t, dddx_t];
end

