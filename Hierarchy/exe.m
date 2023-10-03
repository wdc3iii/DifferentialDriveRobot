%% Solve the problem
problem_setup

%% Solve High Level Problem
%TODO:

sol = [
    2, 3, 4, 5, 6, 6, 6, 6, 7, 8, 9, 9, 9
    8, 8, 8, 8, 8, 7, 6, 5, 5, 5, 5, 4, 3
];
% sol = [1, 2, 3; 1, 1, 1];

x0h = [2, 8];
x0m = [1.25; 7.25; 0.01; 0; 0; 0]; 
x0l = [1.245; 7.247; 0.01; 0.01; 0];

% x0h = [1, 1];
% x0m = [0.5; 0.5; 0.01; 0; 0; 0]; 
% x0l = [0.5; 0.5; 0; 0.01; 0];

%% Solve Mid Level Problem
x_m = x0m;
x_m_hor = {};
x_l = x0l;

cur_set = sol(:, 1);
next_set = sol(:, 2);
xcurr_m = x0m;
xcurr_l = x0l;
% Outer loop through high-level events
high_loops = 1;
midLoops = 0;
while true
    % Construct mid-level constraints
    hx = [max(cur_set(1), next_set(1)); max(cur_set(2), next_set(2)); -min(cur_set(1) - 1, next_set(1) - 1); -min(cur_set(2) - 1, next_set(2) - 1)];
    X = Polyhedron('A', [Hx; Hv; Ha], 'b', [hx; hv; ha]);
    Xf = ctrlInvs{next_set(1), next_set(2)};
    xf = [next_set(1) - 0.5; next_set(2) - 0.5; 0; 0; 0; 0];
    % Middle loop through MPC iterations
    Nf = N;
    while any(Xf.A * xcurr_m > Xf.b)
        disp(midLoops)
        [x_m_T, u_m_T] = midlev(xcurr_m, xf, A, B, X, U, Xf, Q, Qf, R, N, Nf);
        % Nf = Nf - 1;
        
        %% Solve Low Level Problem
        [~, x_l_T] = ode45(@(t, x) lowlev(t, x, xcurr_m, params), [0, T], xcurr_l);

        %% Record
        x_l = [x_l, x_l_T(2:end, :)'];
        xcurr_l = x_l_T(end, :)';
        x_m_hor = [x_m_hor, x_m_T];
        xcurr_m = x_m_T(:, 2);
        x_m = [x_m, xcurr_m];
        midLoops = midLoops + 1;
    end
    high_loops = high_loops + 1;
    cur_set = sol(:, high_loops);
    if high_loops + 1 <= size(sol, 2)
        next_set = sol(:, high_loops + 1);
    else
        break;
    end
end

%% Visualize
figure()
hold on

for r = 0:10
    plot([0, 10], [r, r], 'k')
    plot([r, r], [0, 10], 'k')
end
% axis square
for ii = 1:size(obst, 2)
    r = obst(2, ii);
    c = obst(1, ii);
    rectangle('Position', [c-1, r-1, 1, 1], 'FaceColor', 'k')
end
plot(x_m(1, :), x_m(2, :), '.')
plot(x_l(1, :), x_l(2, :))

figure()
subplot(4, 1, 1)
t = linspace(0, T * size(x_m, 2), size(x_m, 2));
plot(t, x_m(1:2, :))
legend('x', 'y')

subplot(4, 1, 2)
plot(t, x_m(3:4, :))
legend('$\dot{x}$', '$\dot{y}$')

subplot(4, 1, 3)
plot(t, x_m(5:6, :))
legend('$u_x$', '$u_y$')

% subplot(4, 1, 4)
% plot(t(1:end-1), u)
% legend('$\Delta u_x$', '$\Delta u_y$')


