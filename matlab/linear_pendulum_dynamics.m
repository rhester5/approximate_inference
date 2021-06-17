function [A, B] = linear_pendulum_dynamics(x, b)

A = [0 1; -cos(x(1)) -b];
% A = [0, 1; -2*x(1)*x(2) - cos(x(1)), -x(1)^2 + 1];
B = [0; 1];

end

