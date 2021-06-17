function [c, ceq] = collocation_constraints(X, x0 ,xd, num_steps, dt, b, simple)

c = [];
ceq = zeros(1, 2*num_steps+2);

theta = X(1:num_steps);
thetadot = X(num_steps+1:2*num_steps);
u = X(2*num_steps+1:end);

% constrain final state
ceq(1) = theta(end) - xd(1);
ceq(2) = thetadot(end) - xd(2);

% constrain initial state
ceq(3) = theta(1) - x0(1);
ceq(4) = thetadot(1) - x0(2);

% defect constraints
i = 5;
xvec = [theta; thetadot];
for t = 1:num_steps-1
    xt = xvec(:, t);
    xt1 = xvec(:, t+1);
    ut = u(t);
    ut1 = u(t+1);
    if simple
        ft = simple_pendulum_dynamics(xt, ut, b);
        ft1 = simple_pendulum_dynamics(xt1, ut1, b);
    else
        ft = pendulum_dynamics(xt, ut);
        ft1 = pendulum_dynamics(xt1, ut1);
    end
    xtc  = 0.5*(xt+xt1) + dt/8 * (ft - ft1);
    utc = 0.5*(ut+ut1);
    if simple
        ftc = simple_pendulum_dynamics(xtc, utc, b);
    else
        ftc = pendulum_dynamics(xtc, utc);
    end
    defect = (xt-xt1) + dt/6 * (ft + 4*ftc + ft1);
    ceq(i) = defect(1);
    ceq(i+1) = defect(2);
    i = i + 2;
end
