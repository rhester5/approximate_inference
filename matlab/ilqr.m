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
num_iter = 10;

Q = 100*eye(xn); %*100; % kind of works with * 1
R = 10*eye(un);
Qf = 100*eye(xn); %*100; % kind of works * 1

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

[x_d, u_d] = direct_collocation(x0, xd, u0, q, r, num_steps, dt, b, k, l, pendulum, cartpole, drones);
% x_d = get_xd(x0, xd, dt, num_steps);
% u_d = zeros(un, 100);

disp('Computed desired trajectory');

% compute nominal trajectory using direct collocation

% [x0, u0] = direct_collocation(x0, xd2, u0, q, r, num_steps, dt, b, k, l, pendulum, cartpole, drones);
x0 = zeros(xn, 100);
u0 = zeros(un, 100);

disp('Computed nominal trajectory');

x = cell([1, num_iter+1]);
u = cell([1, num_iter+1]);

for it=1:num_iter+1
   x{it} = x0;
   u{it} = u0; 
end

alpha = 10.^linspace(0, -3, 11);
% alpha = 10.^linspace(1, -5, 101);

tStart = cputime;
tIter = tStart;

for it = 1:num_iter

msg = ['iLQR iteration #', num2str(it)];
disp(msg);

% compute improved nominal control sequence by linearizing dynamics about
% nominal trajectory and performing Ricatti recursion

dx = x_d - x{it};
du = u_d - u{it};
S2 = zeros(xn, xn*num_steps+1);
S2(:, xn*num_steps+1:xn*num_steps+xn) = Qf;
S1 = zeros(xn, num_steps+1);
S1(:, num_steps+1) = -2 * Qf * dx(:, end);
S0 = zeros(1, num_steps+1);x0 = [0, 0];
S0(num_steps+1) = dx(:, end)' * Qf * dx(:, end);

for i = 1:num_steps
    t = num_steps-i+1;
    if pendulum
        [At, Bt] = linear_pendulum_dynamics(x{it}(:, t), b);
    elseif cartpole
        At = A_cartpole(x{it}(:, t), u{it}(:, t));
        Bt = B_cartpole(x{it}(:, t), u{it}(:, t));
    elseif drones
        At = A_n_drone(x{it}(:, t), u{it}(:, t), k, l);
        Bt = B_n_drone(x{it}(:, t), u{it}(:, t));
    end
    S2t = S2(:, xn*t+1:xn*t+xn);
    S1t = S1(:, t+1);
    S0t = S0(t+1);
    S2dot = Q - S2t * Bt * inv(R) * Bt' * S2t + S2t * At + At' * S2t;
    S1dot = -2 * Q * dx(:, t) + (At' - S2t * Bt * inv(R) * Bt') * S1t + 2 * S2t * Bt * du(:, t);
    S0dot = dx(:, t)' * Q * dx(:, t) - 0.25 * S1t' * Bt * inv(R) * Bt' * S1t + S1t' * Bt * du(:, t);
    S2(:, xn*t-xn+1:xn*t) = S2t + dt * S2dot;
    S1(:, t) = S1t + dt * S1dot;
    S0(t) = S0t + dt * S0dot;
end

disp('Completed backward pass');

% compute improved nominal trajectory by forward simulating improved
% nominal control sequence

best_cost = cost(x{it}, u{it}, x_d, u_d, Q, R, Qf, num_steps);
no_change = true;
for a = alpha
% for a = [1 1 1 1 1]
    new_x = zeros(xn, num_steps);
    new_u = zeros(un, num_steps);
    % new_x = x{it};
    % new_u = u{it};
    for t = 1:num_steps-1
        if pendulum
            [At, Bt] = linear_pendulum_dynamics(x{it}(:, t), b);
        elseif cartpole
            At = A_cartpole(x{it}(:, t), u{it}(:, t));
            Bt = B_cartpole(x{it}(:, t), u{it}(:, t));
        elseif drones
            At = A_n_drone(x{it}(:, t), u{it}(:, t), k, l)
            Bt = B_n_drone(x{it}(:, t), u{it}(:, t));
        end
        % u{it+1}(:, t) = u_d(:, t) - inv(R) * Bt' * (S2(:, xn*t-xn+1:xn*t) * (x{it+1}(:, t) - x{it}(:, t)) + 0.5 * S1(:, t));
        new_u(:, t) = u_d(:, t) - inv(R) * Bt' * (S2(:, xn*t-xn+1:xn*t) * (new_x(:, t) - x{it}(:, t)) + a * 0.5 * S1(:, t));
        % if pendulum
        %     xdot = pendulum_dynamics(x{it+1}(:, t), u{it+1}(:, t));
        % elseif cartpole
        %     xdot = cartpole_dynamics(x{it+1}(:, t), u{it+1}(:, t));
        % elseif drones
        %     xdot = drone_dynamics(x{it+1}(:, t), u{it+1}(:, t), k, l);
        % end
        % x{it+1}(:, t+1) = x{it+1}(:, t) + dt * xdot;
        if pendulum
            xdot = pendulum_dynamics(new_x(:, t), new_u(:, t));
        elseif cartpole
            xdot = cartpole_dynamics(new_x(:, t), new_u(:, t));
        elseif drones
            xdot = drone_dynamics(new_x(:, t), new_u(:, t), k, l);
        end
        new_x(:, t+1) = new_x(:, t) + dt * xdot;
    end
    if cost(new_x, new_u, x_d, u_d, Q, R, Qf, num_steps) < best_cost
        best_cost = cost(new_x, new_u, x_d, u_d, Q, R, Qf, num_steps);
        x{it+1} = new_x;
        u{it+1} = new_u;
        no_change = false;
    else
        x{it+1} = new_x;
        u{it+1} = new_u;
    end
end

if no_change
    x{it+1} = x{it};
    u{it+1} = u{it};
    % x{it+1} = new_x;
    % u{it+1} = new_u;
end

disp('Completed forward pass');

if pendulum
    msg = ['Cost:', num2str(x_cost(x{it+1}, u{it+1}, x_d, u_d, Q, R, Qf, num_steps))];
    disp(msg)
    % tIter = cputime - tIter;
    % disp(tIter)
    % tIter = cputime;
    % msg = ['Theta Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    % disp(msg);
    % msg = ['Final Theta Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    % disp(msg)
    % msg = ['Omega Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    % disp(msg);
    % msg = ['Final Omega Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    % disp(msg)

    % msg = ['Control Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    % disp(msg);
    % msg = ['Final Control Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    % disp(msg);

elseif cartpole
    msg = ['Position Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Position Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    disp(msg)
    msg = ['Theta Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    disp(msg);
    msg = ['Final Theta Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    disp(msg)
    msg = ['Velocity Error:', num2str(sum(abs(x_d(3, :)-x{it+1}(3, :))))];
    disp(msg);
    msg = ['Final Velocity Error:', num2str(x_d(3, end)-x{it+1}(3, end))];
    disp(msg)
    msg = ['Omega Error:', num2str(sum(abs(x_d(4, :)-x{it+1}(4, :))))];
    disp(msg);
    msg = ['Final Omega Error:', num2str(x_d(4, end)-x{it+1}(4, end))];
    disp(msg)

    msg = ['Control Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Control Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    disp(msg)

elseif drones
    msg = ['x1 Error:', num2str(sum(abs(x_d(1, :)-x{it+1}(1, :))))];
    disp(msg);
    msg = ['Final x1 Error:', num2str(x_d(1, end)-x{it+1}(1, end))];
    disp(msg)
    msg = ['y1 Error:', num2str(sum(abs(x_d(2, :)-x{it+1}(2, :))))];
    disp(msg);
    msg = ['Final y1 Error:', num2str(x_d(2, end)-x{it+1}(2, end))];
    disp(msg)
    msg = ['x2 Error:', num2str(sum(abs(x_d(3, :)-x{it+1}(3, :))))];
    disp(msg);
    msg = ['Final x2 Error:', num2str(x_d(3, end)-x{it+1}(3, end))];
    disp(msg)
    msg = ['y2 Error:', num2str(sum(abs(x_d(4, :)-x{it+1}(4, :))))];
    disp(msg);
    msg = ['Final y2 Error:', num2str(x_d(4, end)-x{it+1}(4, end))];
    disp(msg)
    msg = ['x3 Error:', num2str(sum(abs(x_d(5, :)-x{it+1}(5, :))))];
    disp(msg);
    msg = ['Final x3 Error:', num2str(x_d(5, end)-x{it+1}(5, end))];
    disp(msg)
    msg = ['y3 Error:', num2str(sum(abs(x_d(6, :)-x{it+1}(6, :))))];
    disp(msg);
    msg = ['Final y3 Error:', num2str(x_d(6, end)-x{it+1}(6, end))];
    disp(msg)

    msg = ['Control x1 Error:', num2str(sum(abs(u_d(1, :)-u{it+1}(1, :))))];
    disp(msg);
    msg = ['Final Control x1 Error:', num2str(u_d(1, end) - u{it+1}(1, end))];
    disp(msg)
    msg = ['Control y1 Error:', num2str(sum(abs(u_d(2, :)-u{it+1}(2, :))))];
    disp(msg);
    msg = ['Final Control y1 Error:', num2str(u_d(2, end) - u{it+1}(2, end))];
    disp(msg)
    msg = ['Control x2 Error:', num2str(sum(abs(u_d(3, :)-u{it+1}(3, :))))];
    disp(msg);
    msg = ['Final Control x2 Error:', num2str(u_d(3, end) - u{it+1}(3, end))];
    disp(msg)
    msg = ['Control y2 Error:', num2str(sum(abs(u_d(4, :)-u{it+1}(4, :))))];
    disp(msg);
    msg = ['Final Control y2 Error:', num2str(u_d(4, end) - u{it+1}(4, end))];
    disp(msg)
    msg = ['Control x3 Error:', num2str(sum(abs(u_d(5, :)-u{it+1}(5, :))))];
    disp(msg);
    msg = ['Final Control x3 Error:', num2str(u_d(5, end) - u{it+1}(5, end))];
    disp(msg)
    msg = ['Control y3 Error:', num2str(sum(abs(u_d(6, :)-u{it+1}(6, :))))];
    disp(msg);
    msg = ['Final Control y3 Error:', num2str(u_d(6, end) - u{it+1}(6, end))];
    disp(msg)

end

end

tTotal = cputime - tStart;
disp(tTotal)
disp(tTotal/num_iter)

T = num_steps * dt;
time = linspace(0, T, num_steps);

if pendulum

figure(1)
hold on
for i = 1:num_iter+1
plot(x{i}(1, :), x{i}(2, :));
end
plot(x{num_iter+1}(1, :), x{num_iter+1}(2, :), 'LineWidth', 5);
plot(x_d(1, :), x_d(2, :), 'LineWidth', 2.5);
xlabel('theta (rad)');
ylabel('thetadot (rad/s)');

figure(2)
hold on
for i = 1:num_iter+1
plot(time, u{i}(1, :));
end
plot(time, u{num_iter+1}(1, :), 'LineWidth', 5);
plot(time, u_d(1, :), 'LineWidth', 2.5);
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
