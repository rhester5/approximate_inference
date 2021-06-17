function [theta, thetadot] = simulate_pendulum(u, theta0, thetadot0, num_iter, dt)

theta = zeros(1, num_iter+1);
thetadot = zeros(1, num_iter+1);

theta(1) = theta0;
thetadot(1) = thetadot0;

    for t = 1:num_iter
        thetaddot = u(t) - thetadot(t) - sin(theta(t));
        theta(t+1) = theta(t) + dt * thetadot(t);
        thetadot(t+1) = thetadot(t) + dt * thetaddot;
    end

end