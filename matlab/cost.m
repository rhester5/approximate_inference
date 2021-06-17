function cost = cost(x, u, xd, ud, Q, R, Qf, num_steps)

cost = 0;
for i =1:num_steps-1
    cost = cost + (xd(:, i) - x(:, i))' * Q * (xd(:, i) - x(:, i)) + (ud(:, i) - u(:, i))' * R * (ud(:, i) - u(:, i));
end
cost = cost + (xd(:, end) - x(:, end))' * Qf * (xd(:, end) - x(:, end));

end