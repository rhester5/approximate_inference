function cost = cost_function(X, xn, un, Q, R, num_steps)

x = zeros(xn, num_steps);
u = zeros(un, num_steps);

for i = 0:xn-1
x(i+1, :) = X(i*num_steps+1:(i+1)*num_steps);
end

for i = 1:un
u(i, :) = X((xn+i-1)*num_steps+1:(xn+i)*num_steps);
end

cost = 0;
for i =1:num_steps
    cost = cost + x(:, i)' * Q * x(:, i) + u(:, i)' * R * u(:, i);
end

end