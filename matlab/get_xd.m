function traj = get_xd(x0, xd, dt, num_steps)
	xn = size(x0);
	xn = xn(2);
	step = (xd-x0)/num_steps;
	traj = zeros(xn, num_steps);
	traj(:, 1) = x0;
	for i =2:num_steps
		traj(1:xn/2, i) = traj(1:xn/2, i-1)' + step(1:xn/2);
		traj(xn/2+1:end, i) = step(1:xn/2)/dt;
end