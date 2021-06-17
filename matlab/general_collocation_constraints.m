function [c, ceq] = general_collocation_constraints(X, x0 ,xd, u0, num_steps, dt, damping, spring_constant, rest_length, pendulum, cartpole, drones)

xn = size(x0);
xn = xn(2);

un = size(u0);
un = un(2);

c = [];
ceq = zeros(1, xn*num_steps+xn);

x = zeros(xn, num_steps);
u = zeros(un, num_steps);

for i = 0:xn-1
x(i+1, :) = X(i*num_steps+1:(i+1)*num_steps);
end

for i = 1:un
u(i, :) = X((xn+i-1)*num_steps+1:(xn+i)*num_steps);
end

% constrain final state
for i = 1:xn
ceq(i) = x(i, end) - xd(i);
end

% constrain initial state
for i = 1:xn
ceq(xn+i) = x(i, 1) - x0(i);
end

% defect constraints
i = 2*xn+1;
for t = 1:num_steps-1
    xt = x(:, t);
    xt1 = x(:, t+1);
    ut = u(:, t);
    ut1 = u(:, t+1);
    if pendulum
        ft = pendulum_dynamics(xt, ut);
        ft1 = pendulum_dynamics(xt1, ut1);
    elseif cartpole
        ft = cartpole_dynamics(xt, ut);
        ft1 = cartpole_dynamics(xt1, ut1);
    elseif drones
        ft = drone_dynamics(xt, ut, spring_constant, rest_length);
        ft1 = drone_dynamics(xt1, ut1, spring_constant, rest_length);
    end
    xtc  = 0.5*(xt+xt1) + dt/8 * (ft - ft1);
    utc = 0.5*(ut+ut1);
    if pendulum
        ftc = pendulum_dynamics(xtc, utc);
    elseif cartpole
        ftc = cartpole_dynamics(xtc, utc);
    elseif drones
        ftc = drone_dynamics(xtc, utc, spring_constant, rest_length);
    end
    defect = (xt-xt1) + dt/6 * (ft + 4*ftc + ft1);
    for j = 0:xn-1
    ceq(i+j) = defect(j+1);
    end
    i = i + xn;
end
end