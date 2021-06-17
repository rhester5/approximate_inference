function xdot = pendulum_dynamics(x, u)

xdot = zeros(2, 1);
xdot(1) = x(2);
xdot(2) = u - (x(1)^2 - 1) * x(2) - sin(x(1));
end