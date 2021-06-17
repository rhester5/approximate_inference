function [c, ceq] = nonlinear_constraints(u, theta0, thetadot0, theta_d, thetadot_d, num_iter, dt)

[theta, thetadot] = simulate_pendulum(u, theta0, thetadot0, num_iter, dt);

theta_final = theta(end);
thetadot_final = thetadot(end);

c = [];
ceq = [theta_final - theta_d; thetadot_final - thetadot_d];

end