% clear all
close all

pendulum = true;
cartpole = false;
drones = false;

if pendulum
    xn = 2;
    un = 1;
elseif cartpole
    xn = 4;
    un = 1;
elseif drones
    xn = 12;
    un = 6;
end

b = 0.3;
k = 1;
l = 1;

Q = eye(xn); %*100; % kind of works with * 1
R = eye(un);
Qf = eye(xn); %*100; % kind of works * 1

q = eye(xn);
r = eye(un);

num_steps = 100;
dt = 0.1;

% compute desired trajectory using direct collocation

if pendulum
    x0 = [0, 0];
    u0 = [0];
    xd = [pi, 0];
    xd2 = -xd;
elseif cartpole
    x0 = [0, 0, 0, 0];
    u0 = [0];
    xd = [0, pi, 0, 0];
    xd2 = -xd;
elseif drones
    % x0 = [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
    x0 = [-2, -2, 1, 1, 4, 4, 0, 0, 0, 0, 0, 0];
    u0 = zeros(1, un);
    xd = [0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
end

[x_d, u_d] = direct_collocation(x0, xd, u0, q, r, num_steps+1, dt, b, k, l, pendulum, cartpole, drones);
% x_d = get_xd(x0, xd, dt, num_steps);
% u_d = zeros(un, num_steps);
u_d = u_d(1:num_steps);

disp('Computed desired trajectory');

% compute nominal trajectory using direct collocation

[x0, u0] = direct_collocation(x0, xd2, u0, q, r, num_steps+1, dt, b, k, l, pendulum, cartpole, drones);
% x0 = zeros(xn, num_steps+1);
% u0 = zeros(un, num_steps);
u0 = u0(1:num_steps);

disp('Computed nominal trajectory');
Q = dt*Q;
R = dt*R;
Qf = dt*Qf;
[x_traj, u_traj] = demo_pendulum(x0, u0, x_d, u_d, Q, R, Qf, num_steps, dt);

T = num_steps * dt;
time = linspace(0, T, num_steps);

if pendulum

figure(2)
hold on
% for i = 1:num_iter+1
% plot(x{i}(1, :), x{i}(2, :));
% end
% plot(x{num_iter+1}(1, :), x{num_iter+1}(2, :));
plot(x_traj(1, :), x_traj(2, :))
plot(x_d(1, :), x_d(2, :));
xlabel('theta (rad)');
ylabel('thetadot (rad/s)');

figure(3)
hold on
% for i = 1:num_iter+1
% plot(time, u{i}(1, :));
% end
% plot(time, u{num_iter+1}(1, :));
plot(time, u_traj)
plot(time, u_d(1, :));
xlabel('time (s)');
ylabel('control');

elseif cartpole

figure(1)
hold on
% for i = 1:num_iter+1
% plot(x{i}(1, :), x{i}(3, :));
% end
plot(x{num_iter+1}(1, :), x{num_iter+1}(3, :));
plot(x_d(1, :), x_d(3, :));
xlabel('position (m)');
ylabel('velocity (m/s)');

figure(2)
hold on
for i = 1:num_iter+1
plot(x{i}(2, :), x{i}(4, :));
end
% plot(x{num_iter+1}(2, :), x{num_iter+1}(4, :));
plot(x_d(2, :), x_d(4, :));
xlabel('theta (rad)');
ylabel('thetadot (rad/s)');

figure(3)
hold on
for i = 1:num_iter+1
plot(time, u{i}(1, :));
end
% plot(time, u{num_iter+1}(1, :));
plot(time, u_d(1, :));
xlabel('time (s)');
ylabel('control');

elseif drones

% figure(1)
% hold off
% % D1 = rectangle('Position', [x{num_iter+1}(1, 1), x{num_iter+1}(2, 1), 1, 1], 'Curvature', [1, 1]);
% % D2 = rectangle('Position', [x{num_iter+1}(3, 1), x{num_iter+1}(4, 1), 1, 1], 'Curvature', [1, 1]);
% % D3 = rectangle('Position', [x{num_iter+1}(5, 1), x{num_iter+1}(6, 1), 1, 1], 'Curvature', [1, 1]);
% % axis([-10, 10, -10, 10])
% % for i = 1:num_steps
% %     set(D1, 'Position', [x{num_iter+1}(1, i), x{num_iter+1}(2, i), 1, 1]);
% %     set(D2, 'Position', [x{num_iter+1}(3, i), x{num_iter+1}(4, i), 1, 1]);
% %     set(D3, 'Position', [x{num_iter+1}(5, i), x{num_iter+1}(6, i), 1, 1]);
% %     S12 = line([x{num_iter+1}(1, i), x{num_iter+1}(2, i)], [x{num_iter+1}(3, i), x{num_iter+1}(4, i)]);
% %     S23 = line([x{num_iter+1}(3, i), x{num_iter+1}(4, i)], [x{num_iter+1}(5, i), x{num_iter+1}(6, i)]);
% %     F(i) = getframe;
% % end
% D1 = rectangle('Position', [x_d(1, 1), x_d(2, 1), 1, 1], 'Curvature', [1, 1]);
% D2 = rectangle('Position', [x_d(3, 1), x_d(4, 1), 1, 1], 'Curvature', [1, 1]);
% D3 = rectangle('Position', [x_d(5, 1), x_d(6, 1), 1, 1], 'Curvature', [1, 1]);
% axis([-10, 10, -10, 10])
% for i = 1:num_steps
%     F(i) = getframe;
%     set(D1, 'Position', [x_d(1, i), x_d(2, i), 1, 1]);
%     set(D2, 'Position', [x_d(3, i), x_d(4, i), 1, 1]);
%     set(D3, 'Position', [x_d(5, i), x_d(6, i), 1, 1]);
%     S12 = line([x_d(1, i), x_d(2, i)], [x_d(3, i), x_d(4, i)]);
%     S23 = line([x_d(3, i), x_d(4, i)], [x_d(5, i), x_d(6, i)]);
% end
% F(end+1) = getframe;
% movie(F, 10)

figure(2)
hold on
plot(x_d(1, :), x_d(2, :))
plot(x_d(3, :), x_d(4, :))
plot(x_d(5, :), x_d(6, :))
plot(x{num_iter+1}(1, :), x{num_iter+1}(2, :))
plot(x{num_iter+1}(3, :), x{num_iter+1}(4, :))
plot(x{num_iter+1}(5, :), x{num_iter+1}(6, :))
legend('d1d', 'd2d', 'd3d', 'd1', 'd2', 'd3')

% figure(3)
% hold off
% D1 = rectangle('Position', [x{num_iter+1}(1, 1), x{num_iter+1}(2, 1), 1, 1], 'Curvature', [1, 1]);
% D2 = rectangle('Position', [x{num_iter+1}(3, 1), x{num_iter+1}(4, 1), 1, 1], 'Curvature', [1, 1]);
% D3 = rectangle('Position', [x{num_iter+1}(5, 1), x{num_iter+1}(6, 1), 1, 1], 'Curvature', [1, 1]);
% axis([-10, 10, -10, 10])
% for i = 1:num_steps
%     F(i) = getframe;
%     set(D1, 'Position', [x{num_iter+1}(1, i), x{num_iter+1}(2, i), 1, 1]);
%     set(D2, 'Position', [x{num_iter+1}(3, i), x{num_iter+1}(4, i), 1, 1]);
%     set(D3, 'Position', [x{num_iter+1}(5, i), x{num_iter+1}(6, i), 1, 1]);
%     S12 = line([x{num_iter+1}(1, i), x{num_iter+1}(2, i)], [x{num_iter+1}(3, i), x{num_iter+1}(4, i)]);
%     S23 = line([x{num_iter+1}(3, i), x{num_iter+1}(4, i)], [x{num_iter+1}(5, i), x{num_iter+1}(6, i)]);
% end
% F(end+1) = getframe;
% movie(F, 10)

end
