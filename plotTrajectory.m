function plotTrajectory(x, t)
%PLOTTRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
plot3(x(:, 1), x(:, 2), t)
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Time (sec)")
end

