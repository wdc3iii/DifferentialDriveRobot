function u = fblControl(x, traj, params)
ybar = traj(1:2, 1);
ydbar = traj(1:2, 2);
yddbar = traj(1:2, 3);
ydddbar = traj(1:2, 4);

e = x(1:2) - ybar;
ed = [x(4) * cos(x(3)); x(4) * sin(x(3))] - ydbar;

var1 = yddbar - params.Kp * e - params.Kd * ed;
var2 = [-1 / x(4) * sin(x(3)) 1 / x(4) * cos(x(3))];

accel = [cos(x(3)) sin(x(3))] * var1;
edd = - yddbar + [cos(x(3)), -x(4) * sin(x(3)); sin(x(3)), x(4) * cos(x(3))] * [accel; x(5)];
omega_des = var2 * var1;
omegad_des = [1 / x(4)^2 * sin(x(3)) - x(5) / x(4) * cos(x(3)), -1 / x(4)^2 * cos(x(3)) - x(5) / x(4) * sin(x(3))] * var1 + var2 * (ydddbar - params.Kp * ed - params.Kd * edd);
alpha = -params.sigma * (x(5) - omega_des) + omegad_des - 2 * [e', ed'] * params.P * [0; 0; -x(4) * cos(x(3)); x(4) * sin(x(3))];
u = accelToTorque(x, [accel; alpha], params);
end