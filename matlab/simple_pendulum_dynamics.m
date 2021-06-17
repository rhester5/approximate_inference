function xdot = simple_pendulum_dynamics(x, u, b)

xdot = zeros(2, 1);
xdot(1) = x(2);
xdot(2) = u - b * x(2) - sin(x(1));
end