function [x, u] = direct_collocation(x0, xd, u0, Q, R, num_steps, dt, damping, spring_constant, rest_length, pendulum, cartpole, drones)

% parameter set x, u
% x = [theta, thetadot]
% u = [f]
xn = size(x0);
xn = xn(2);

un = size(u0);
un = un(2);

% no linear equality or inequality constraints
A = [];
b = [];
Aeq = [];
beq = [];

X0 = zeros(1, (xn+un) * num_steps); % initial guess is zero state, zero control

lb = ones(1, (xn+un) * num_steps); % lower bound on state, control
if pendulum
	lb(1:2*num_steps) = -inf * lb(1:2*num_steps);
	lb(2*num_steps+1:end) = -5 * lb(2*num_steps+1:end);
elseif cartpole
	lb(1:4*num_steps) = -inf * lb(1:4*num_steps);
	lb(4*num_steps+1:end) = -5 * lb(4*num_steps+1:end);
elseif drones
	lb(1:xn*num_steps) = -inf * lb(1:xn*num_steps);
	lb(xn*num_steps+1:end) = -5 * lb(xn*num_steps+1:end);
end

ub = ones(1, (xn+un) * num_steps);  % upper bound on state, control
if pendulum
	ub(1:2*num_steps) = inf * ub(1:2*num_steps);
	ub(2*num_steps+1:end) = 5 * ub(2*num_steps+1:end);
elseif cartpole
	ub(1:4*num_steps) = inf * ub(1:4*num_steps);
	ub(4*num_steps+1:end) = 5 * ub(4*num_steps+1:end);
elseif drones
	ub(1:xn*num_steps) = inf * ub(1:xn*num_steps);
	ub(xn*num_steps+1:end) = 5 * ub(xn*num_steps+1:end);
end

options = optimoptions(@fmincon, 'TolFun', 0.00000001, 'MaxIter', 10000, ...
            'MaxFunEvals', 100000, 'Display', 'iter', ...
            'DiffMinChange', 0.001, 'Algorithm', 'sqp');

fun = @(X)cost_function(X, xn, un, Q, R, num_steps);
nonlcon = @(X)general_collocation_constraints(X, x0 ,xd, u0, num_steps, dt, damping, spring_constant, rest_length, pendulum, cartpole, drones);
        
X = fmincon(fun, X0, A, b, Aeq, beq, lb, ub, nonlcon, options);

x = zeros(xn, num_steps);
u = zeros(un, num_steps);

for i = 0:xn-1
x(i+1, :) = X(i*num_steps+1:(i+1)*num_steps);
end

for i = 1:un
u(i, :) = X((xn+i-1)*num_steps+1:(xn+i)*num_steps);
end

end