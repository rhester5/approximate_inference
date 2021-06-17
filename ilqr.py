import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from scipy.stats import multivariate_normal
import copy

from quadrotor import discrete_quadrotor_dynamics, quadrotor_dynamics, linear_quadrotor_dynamics, single_quadrotor_dynamics, linear_single_quadrotor_dynamics

def linear_pendulum_dynamics(x, u):#, dt):
	b = 0.3
	A = np.array([[0, 1], [-np.cos(x[0]), -b]])
	B = np.array([[0], [1]])
	# A = np.array([[1, 0.1], [-np.cos(x[0]), -b]])
	# B = np.array([[0], [0.1]])
	# A = np.array([[1, 0.1], [-np.cos(x[0]), 1]])
	# B = np.array([[0], [0.1]])
	return A, B
	# return np.eye(A.shape[0]) + dt*A, dt*B

def pendulum_dynamics(x, u):
	xdot = np.array([x[1], u[0] - (x[0]**2 - 1) * x[1] - np.sin(x[0])])
	return xdot

def linear_cartpole_dynamics(x, u):
	# TODO is this right?

	M = np.array([[2, np.cos(x[1])], [np.cos(x[1]), 1]])
	C = np.array([[0, -np.sin(x[1])], [0, 0]])
	taug = np.array([[0], [-np.sin(x[1])]])
	dtaug_dq = np.array([[0, 0], [0, -np.cos(x[1])]])
	B = np.array([[1], [0]])

	Minv = np.linalg.inv(M)
	Minv_dtaug_dq = np.dot(Minv, dtaug_dq)
	Minv_B = np.dot(Minv, B)

	A = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [Minv_dtaug_dq[0, 0], Minv_dtaug_dq[0, 1], 0, 0], [Minv_dtaug_dq[1, 0], Minv_dtaug_dq[1, 1], 0, 0]])
	B = np.array([[0], [0], [Minv_B[0, 0]], [Minv_B[1, 0]]])

	return A, B

def cartpole_dynamics(x, u):
	xdot = np.zeros(x.shape)
	xdot[0] = x[2]
	xdot[1] = x[3]
	xdot[2] = (u[0] + np.sin(x[1]) * (x[3]**2 + np.cos(x[1]))) / (1 + np.sin(x[1])**2)
	xdot[3] = (-u[0] * np.cos(x[1]) - x[3]**2 * np.cos(x[1]) * np.sin(x[1]) - 2 * np.sin(x[1])) / (1 + np.sin(x[1])**2)
	return xdot

def get_pendulum_position(x):
	l = 5
	pendulum_x = l * np.sin(x[0])
	pendulum_y = -l * np.cos(x[0])
	return pendulum_x, pendulum_y

def get_cartpole_pendulum_position(x):
	l = 5
	pendulum_x = x[0] + l * np.sin(x[1])
	pendulum_y = -l * np.cos(x[1])
	return pendulum_x, pendulum_y

def aico(dynamics, linear_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf):
	# x0 is initial trajectory
	# alpha is convergence rate
	# A is linear state dynamics
	# a is Gaussian noise
	# B is linear control dynamics
	# Q is covariance
	# R is Q in LQR
	# r is R as a vector instead of a diagonal matrix
	# H is R in LQR

	q = np.diag(Q)
	qf = np.diag(Qf)
	c = np.zeros((x0.shape[0],))
	C = np.eye(x0.shape[0])*0.1 # TODO covariance should be an input
	alpha = 0.9

	# x_traj0 = np.zeros((x0.shape[0], num_steps))
	# x_traj0[:, 0] = x0
	# x_traj0[:, -1] = xd
	# for i in range(1, num_steps):
	# 	x_traj0[:, i] = xd

	# u_traj0 = np.zeros((u0.shape[0], num_steps))
	# u_traj0[:, 0] = u0
	# u_traj0[:, -1] = ud
	# for i in range(1, num_steps):
	# 	u_traj0[:, i] = ud

	x_traj0 = xd
	u_traj0 = ud

	x_trajectories = [x_traj0]
	u_trajectories = [u_traj0]

	error = np.inf
	it = 1

	# a = []
	# for i in range(num_steps):
	# 	a.append(np.random.multivariate_normal(c, C))

	# s = num_steps * [x0]
	# Sinv = num_steps * [np.eye(x0.shape[0]) * 1e10]
	# S = num_steps * [np.eye(x0.shape[0]) * 1e-10]
	# v = np.zeros((x0.shape[0], num_steps))
	# Vinv = num_steps * [np.zeros((x0.shape[0], x0.shape[0]))]
	# V = num_steps * [np.zeros((x0.shape[0], x0.shape[0]))]
	# b = num_steps * [np.inf]

	# a = []
	s = []
	Sinv = []
	S = []
	v = []
	Vinv = []
	V = []
	b = []
	for i in range(num_steps+1):
		# a.append(np.random.multivariate_normal(c, C))
		# a.append(x0)
		s.append(x0)
		Sinv.append(np.eye(x0.shape[0]) * 1e10)
		S.append(np.eye(x0.shape[0]) * 1e-10)
		v.append(np.zeros((x0.shape[0],)))
		Vinv.append(np.zeros((x0.shape[0], x0.shape[0])))
		V.append(np.zeros((x0.shape[0], x0.shape[0])))
		b.append(np.ones((x0.shape[0],)))

	while error > 1e-3 and it < max_num_iter:
		# print(it)
		# x_traj = np.zeros(x_traj0.shape)
		# u_traj = np.zeros(u_traj0.shape)
		x_traj = copy.deepcopy(x_trajectories[it-1])
		u_traj = copy.deepcopy(u_trajectories[it-1])

		a = []
		for i in range(num_steps):
			a.append(np.random.multivariate_normal(c, C))

		# forward sweep
		for t in range(1, num_steps):
			# print(t)
			count = 0
			while np.linalg.norm(x_traj[:, t] - b[t]) > 1e0:
				# print(count, np.linalg.norm(x_traj[:, t] - b[t]), b[t])
				Atm1, Btm1 = linear_dynamics(x_traj[:, t-1], u_traj[:, t-1], dt)
				s[t] = a[t-1] + Atm1 @ np.linalg.inv((np.linalg.inv(S[t-1]) + Q)) @ (Sinv[t-1] @ s[t-1] + q)
				S[t] = C + Btm1 @ np.linalg.inv(R) @ Btm1.T + Atm1 @ np.linalg.inv(Sinv[t-1] + Q) @ Atm1.T
				Sinv[t] = np.linalg.inv(S[t])
				msg_m = s[t]
				# msg_m = np.random.multivariate_normal(s[t], S[t])
				# print(msg_m)
				if it == 1:
					x_traj[:, t] = s[t]
				else:
					x_traj[:, t] = (1 - alpha) * x_traj[:, t] + alpha * b[t]
				At, Bt = linear_dynamics(x_traj[:, t], u_traj[:, t], dt)
				msg = q
				# msg = np.random.multivariate_normal(q, Q)
				v[t] = -np.linalg.inv(At) @ a[t] + np.linalg.inv(At) @ np.linalg.inv(Vinv[t+1] + Q) @ (Vinv[t+1] @ v[t+1] + q)
				V[t] = np.linalg.inv(At) @ (C + Bt @ np.linalg.inv(R) @ Bt.T + np.linalg.inv(Vinv[t+1] + Q)) @ At.T
				Vinv[t] = np.linalg.inv(V[t])
				# TODO is it inv At.T or just At.T
				msg_p =  v[t]
				# msg_p = np.random.multivariate_normal(v[t], V[t])
				# b[t] = msg_m * msg_p * msg
				b[t] = np.multiply(np.multiply(msg_m, msg_p), msg)
				a[t] = np.random.multivariate_normal(c, C)
				count+=1
				# u_traj[:, t] = ud[:, t] - np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ x_traj[:, t] + a[t]) - v[t+1])

		# for t in range(num_steps-1):
			# x_traj[:, t+1] = x_traj[:, t] + dynamics(x_traj[:, t], u_traj[:, t])
		
		# backward sweep
		for i in range(num_steps-1):
			t = num_steps-i-1
			# print(t)
			count = 0
			while np.linalg.norm(x_traj[:, t] - b[t]) > 1e0:
				# print(count, np.linalg.norm(x_traj[:, t] - b[t]), b[t])
				Atm1, Btm1 = linear_dynamics(x_traj[:, t-1], u_traj[:, t-1], dt)
				s[t] = a[t-1] + Atm1 @ np.linalg.inv((np.linalg.inv(S[t-1]) + Q)) @ (Sinv[t-1] @ s[t-1] + q)
				S[t] = C + Btm1 @ np.linalg.inv(R) @ Btm1.T + Atm1 @ np.linalg.inv(Sinv[t-1] + Q) @ Atm1.T
				Sinv[t] = np.linalg.inv(S[t])
				msg_m = s[t]
				# msg_m = np.random.multivariate_normal(s[t], S[t])
				# print(msg_m)
				if it == 1:
					x_traj[:, t] = s[t]
				else:
					x_traj[:, t] = (1 - alpha) * x_traj[:, t] + alpha * b[t]
				At, Bt = linear_dynamics(x_traj[:, t], u_traj[:, t], dt)
				msg = q
				# msg = np.random.multivariate_normal(q, Q)
				v[t] = -np.linalg.inv(At) @ a[t] + np.linalg.inv(At) @ np.linalg.inv(Vinv[t+1] + Q) @ (Vinv[t+1] @ v[t+1] + q)
				V[t] = np.linalg.inv(At) @ (C + Bt @ np.linalg.inv(R) @ Bt.T + np.linalg.inv(Vinv[t+1] + Q)) @ At.T
				Vinv[t] = np.linalg.inv(V[t])
				# TODO is it inv At.T or just At.T
				msg_p =  v[t]
				# msg_p = np.random.multivariate_normal(v[t], V[t])
				# b[t] = msg_m * msg_p * msg
				b[t] = np.multiply(np.multiply(msg_m, msg_p), msg)
				a[t] = np.random.multivariate_normal(c, C)
				count+=1
				# u_traj[:, t] = ud[:, t] - np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ x_traj[:, t] + a[t]) - v[t+1])
				# print(np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ x[:, t] + a[t]) - v[t+1]))
		
		# print(x_traj)
		# for t in range(num_steps-1):
		# 	x_traj[:, t+1] = x_traj[:, t] + dynamics(x_traj[:, t], u_traj[:, t])
		# print(x_traj)

		x_trajectories.append(x_traj)
		u_trajectories.append(u_traj)

		error = np.linalg.norm(x_traj[:, t] - x_trajectories[it-1][:, t])
		print(error)
		it += 1

	return x_trajectories, u_trajectories

def ilqg(dynamics, linear_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf):
	# x0 is initial trajectory
	# alpha is convergence rate
	# A is linear state dynamics
	# a is Gaussian noise
	# B is linear control dynamics
	# Q is covariance
	# R is Q in LQR
	# r is R as a vector instead of a diagonal matrix
	# H is R in LQR

	q = np.diag(Q)
	qf = np.diag(Qf)
	c = np.zeros((x0.shape[0],))
	C = np.eye(x0.shape[0])*0.01 # TODO covariance should be an input
	# C[0, 0] = 0.01
	# C[1, 1] = 0.1
	alpha = 0.9

	x_traj0 = np.zeros((x0.shape[0], num_steps))
	x_traj0[:, 0] = x0

	u_traj0 = np.zeros((u0.shape[0], num_steps))
	u_traj0[:, 0] = u0

	x_trajectories = [x_traj0]
	u_trajectories = [u_traj0]

	error = np.inf
	it = 1

	# a = []
	# for i in range(num_steps):
	# 	a.append(np.random.multivariate_normal(c, C))

	while error > 1e-3 and it < max_num_iter:

	# compute improved nominal control sequence by linearizing dynamics about nominal trajectory and performing Ricatti recursion

		x_traj = np.zeros(x_traj0.shape)
		u_traj = np.zeros(u_traj0.shape)

		V = [Qf]
		v = [qf]
		# v = [2 * Qf @ (x_traj[:, -1] - xd)]

		a = []
		for i in range(num_steps):
			a.append(np.random.multivariate_normal(c, C))
		
		for i in range(num_steps-1):
			t = num_steps-i-1
			# q = 2 * Q @ (x_traj[:, t] - xd)
			At, Bt = linear_dynamics(x_trajectories[it-1][:, t], u_trajectories[it-1][:, t], dt)
			# TODO should this be xd-x, ud-u?
			Vtp1 = V[i]
			vtp1 = v[i]
			# a = np.random.multivariate_normal(c, C)
			# K = At.T @ Vtp1.T @ np.linalg.inv(Vtp1 + Bt.T @ R @ np.linalg.inv(Bt))
			# K = At.T @ Vtp1.T @ np.linalg.inv(Vtp1 + Bt @ R @ Bt.T)
			K = At.T @ Vtp1.T @ np.linalg.inv(Vtp1 + Bt @ np.linalg.inv(R) @ Bt.T)
			# TODO thinking it might be Rinv instead of R because that's how it is later in the paper and what's in the paper is completely wrong
			Vt = Q + At.T @ Vtp1 @ At - K @ Vtp1 @ At
			# vt = q + At.T @(vtp1 - Vtp1 @ a) - K @ (vtp1 - Vtp1 @ a)
			vt = q + At.T @(vtp1 - Vtp1 @ a[t]) - K @ (vtp1 - Vtp1 @ a[t])
			# print(Vt)
			# print(vt)
			V.append(Vt)
			v.append(vt)

		V = V[::-1]
		v = v[::-1]

		# compute improved nominal trajectory by forward simulating improved nominal control sequence

		for t in range(num_steps-1):
			# At, Bt = linear_dynamics(x_traj[:, t], u_traj[:, t])
			At, Bt = linear_dynamics(x_trajectories[it-1][:, t], u_trajectories[it-1][:, t], dt)
			# TODO should this actually be the same At, Bt as before?
			# a = np.random.multivariate_normal(c, C)
			# TODO should these a_t be precomputed? are they actually Gaussian noise or am I missing something?
			# u_traj[:, t] = -np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ x_traj[:, t] + a[t]) - v[t+1])
			# u_traj[:, t] = -np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ (xd - x_traj[:, t]) + a) - v[t+1])
			# u_traj[:, t] = ud[:, t]-np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ (xd[:, t] - x_traj[:, t]) + a[t]) - v[t+1])
			u_traj[:, t] = ud[:, t]-np.linalg.inv(R + Bt.T @ V[t+1] @ Bt) @ Bt.T @ (V[t+1] @ (At @ x_traj[:, t] + a[t]) - v[t+1])
			# print(u_traj[:, t])
			# TODO should this be xd-x? and is x_traj[:, t] the right shape?
			x_traj[:, t+1] = (1-alpha) * x_trajectories[it-1][:, t+1] + alpha * (At @ x_traj[:, t] + a[t] + Bt @ u_traj[:, t])
			# x_traj[:, t+1] = (1-alpha) * x_trajectories[it-1][:, t+1] + alpha * (At @ (xd - x_traj[:, t]) + a + Bt @ u_traj[:, t])
			# x_traj[:, t+1] = (1 - alpha) * x_trajectories[it-1][:, t+1] + alpha * (x_traj[:, t] + dt * (At @ (xd - x_traj[:, t]) + a + Bt @ u_traj[:, t]))
			# x_traj[:, t+1] = (1 - alpha) * x_trajectories[it-1][:, t+1] + alpha * (x_traj[:, t] + dt * (At @ (x_traj[:, t] - x_trajectories[it-1][:, t]) + a + Bt @ u_traj[:, t]))
			# x_traj[:, t+1] = (1 - alpha) * x_trajectories[it-1][:, t+1] + alpha * (x_traj[:, t] + dt * (At @ (x_traj[:, t] - x_trajectories[it-1][:, t]) + a[t] + Bt @ u_traj[:, t]))
			# x_traj[:, t+1] = (1-alpha) * x_trajectories[it-1][:, t+1] + alpha * (At @ (x_traj[:, t] - x_trajectories[it-1][:, t]) + a[t] + Bt @ u_traj[:, t])
			# x_traj[:, t+1] = x_traj[:, t] + dt * dynamics(x_traj[:, t], u_traj[:, t])

		x_trajectories.append(x_traj)
		u_trajectories.append(u_traj)

		error = np.linalg.norm(x_traj[:, t] - x_trajectories[it-1][:, t])
		print(error)
		it += 1

	return x_trajectories, u_trajectories

def ilqr(dynamics, linear_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf):

	x_traj0 = np.zeros((x0.shape[0], num_steps))
	x_traj0[:, 0] = x0

	u_traj0 = np.zeros((u0.shape[0], num_steps))
	u_traj0[:, 0] = u0

	x_trajectories = [x_traj0]
	u_trajectories = [u_traj0]

	error = np.inf
	it = 1

	while error > 1e-3 and it < max_num_iter:
		print(it-1)

		# compute improved nominal control sequence by linearizing dynamics about nominal trajectory and performing Ricatti recursion

		x_traj = np.zeros(x_traj0.shape)
		u_traj = np.zeros(u_traj0.shape)
		# x_traj = x_trajectories[it-1]
		# u_traj = u_trajectories[it-1]
		# x_traj = copy.deepcopy(x_trajectories[it-1])
		# u_traj = copy.deepcopy(u_trajectories[it-1])

		dx = xd[:, -1]-x0
		# dx = -dx
		du = ud[:, -1]-u0
		# du = -du

		S2 = [Qf]
		S1 = [-2 * Qf @ dx]
		S0 = [dx.T @ Qf @ dx]
		
		# TODO am I not getting all the way there because I should be doing range(num_steps-1)?
		for i in range(num_steps-1): # TODO or just num_steps?
		# for i in range(num_steps):
			t = num_steps-i-1
			# At, Bt = linear_dynamics(x_traj[:, t], u_traj[:, t])
			# At, Bt = linear_dynamics(x_traj[:, t-1], u_traj[:, t-1])
			At, Bt = linear_dynamics(x_trajectories[it-1][:, t], u_trajectories[it-1][:, t])
			# At, Bt = linear_dynamics(x_trajectories[it-1][:, t-1], u_trajectories[it-1][:, t-1])
			# dx = xd - x_trajectories[it-1][:, t-1]
			# du = ud - u_trajectories[it-1][:, t-1]
			dx = xd[:, t] - x_trajectories[it-1][:, t]
			# dx = -dx
			du = ud[:, t] - u_trajectories[it-1][:, t]
			# du = -du
			# dx = xd - x_traj[:, t]
			# du = ud - u_traj[:, t]
			# dx = xd - x_traj[:, t-1]
			# du = ud - u_traj[:, t-1]
			S2t = S2[i]
			S1t = S1[i]
			S0t = S0[i]
			S2dot = Q - S2t @ Bt @ np.linalg.inv(R) @ Bt.T @ S2t + S2t @ At + At.T @ S2t
			# S1dot = (At.T - S2t @ Bt @ np.linalg.inv(R) @ Bt.T) @ S1t
			S1dot = -2 * Q @ dx + (At.T - S2t @ Bt @ np.linalg.inv(R) @ Bt.T) @ S1t + 2 * S2t @ Bt @ du
			# S0dot = -0.25 * S1t.T @ Bt @ np.linalg.inv(R) @ Bt.T @ S1t
			S0dot = dx.T @ Q @ dx - 0.25 * S1t.T @ Bt @ np.linalg.inv(R) @ Bt.T @ S1t + S1t.T @ Bt @ du

			S2.append(S2t + dt * S2dot)
			S1.append(S1t + dt * S1dot)
			S0.append(S0t + dt * S0dot)

		S2 = S2[::-1]
		S1 = S1[::-1]
		S0 = S0[::-1]

		# compute improved nominal trajectory by forward simulating improved nominal control sequence

		for t in range(num_steps-1):
			# At, Bt = linear_dynamics(x_traj[:, t], u_traj[:, t])
			# At, Bt = linear_dynamics(x_traj[:, t-1], u_traj[:, t-1])
			At, Bt = linear_dynamics(x_trajectories[it-1][:, t], u_trajectories[it-1][:, t])
			# At, Bt = linear_dynamics(x_trajectories[it-1][:, t-1], u_trajectories[it-1][:, t-1])
			u_traj[:, t] = ud[:, t] - np.linalg.inv(R) @ Bt.T @ (S2[t] @ (x_traj[:, t] - x_trajectories[it-1][:, t]) + 0.5 * S1[t])
			# u_traj[:, t] = -np.linalg.inv(R) @ Bt.T @ (S2[t] @ (xd - x_traj[:, t]) + 0.5 * S1[t])
			xdot = dynamics(x_traj[:, t], u_traj[:, t])
			# A, B = linear_dynamics(x_traj[:, t], u_traj[:, t])
			# xdot = A @ x_traj[:, t] + B @ u_traj[:, t]
			# xdot = A @ (xd - x_traj[:, t]) + B @ (ud - u_traj[:, t])
			# TODO it shouldn't be xd - x but x - x_nom
			x_traj[:, t+1] = x_traj[:, t] + dt * xdot
			# x_traj[:, t+1] = x_traj[:, t] + dt * (xd - xdot)

		x_trajectories.append(x_traj)
		u_trajectories.append(u_traj)

		# error = np.linalg.norm(x_traj[:, :] - x_trajectories[it-1][:, :])
		error = np.linalg.norm(u_traj[:, :] - u_trajectories[it-1][:, :])
		# error = np.linalg.norm(xd - x_traj[:, -1])
		print(error)
		it += 1

	return x_trajectories, u_trajectories
	
if __name__ == '__main__':

	# system = 'pendulum'
	# system = 'cartpole'
	# system = 'arm' # TODO
	system = 'quadrotor'
	# system = 'drone rope'
	# algorithm = 'ilqr'
	# algorithm = 'ilqg'
	algorithm = 'aico'

	if system == 'pendulum':

		max_num_iter = 100
		num_steps = 100
		dt = 0.1

		T = num_steps * dt
		time = np.linspace(0, T-dt, num_steps)

		Q = 100*np.eye(2)
		R = 10*np.eye(1)
		Qf = 100*np.eye(2)

		x0 = np.array([0, 0])
		xd = np.array([np.pi, 0])
		u0 = np.array([0])
		ud = np.array([0])

		if algorithm == 'ilqr':
			x_trajectories, u_trajectories = ilqr(pendulum_dynamics, linear_pendulum_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'ilqg':
			x_trajectories, u_trajectories = ilqg(pendulum_dynamics, linear_pendulum_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'aico':
			x_trajectories, u_trajectories = aico(pendulum_dynamics, linear_pendulum_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)


		plt.figure(1)
		legend1 = []
		for i, x_traj in enumerate(x_trajectories):
			plt.plot(x_traj[0, :], x_traj[1, :])
			legend1.append('iter ' + str(i))
		plt.xlabel('theta (rad)')
		plt.ylabel('thetadot (rad/s)')
		# plt.legend(legend1)

		plt.figure(2)
		legend2 = []
		for i, u_traj in enumerate(u_trajectories):
			plt.plot(time, u_traj[0, :])
			legend2.append('iter ' + str(i))
		plt.xlabel('time (s)')
		plt.ylabel('control')
		# plt.legend(legend2)

		x_traj = x_trajectories[-1]

		fig = plt.figure(3)
		fig.set_dpi(100)
		fig.set_size_inches(7, 6.5)

		ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
		base = plt.Rectangle((-1, -1), 2, 2, fc='r')
		mass = plt.Circle((0, -5), 1, fc='b')
		beam = plt.Line2D((0, 0), (0, -5), lw=2.5)

		def init():
			ax.add_patch(base)
			ax.add_patch(mass)
			ax.add_line(beam)
			return base, mass, beam

		def animate(i):
			x, y = get_pendulum_position(x_traj[:, i])
			mass.center = (x, y)
			beam.set_xdata([0, x])
			beam.set_ydata([0, y])
			return base, mass, beam

		anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=num_steps, 
                               interval=60,
                               blit=True)

		plt.show()

	elif system == 'cartpole':

		max_num_iter = 100
		num_steps = 100
		dt = 0.1

		T = num_steps * dt
		time = np.linspace(0, T-dt, num_steps)

		Q = np.eye(4)
		R = np.eye(1)
		Qf = np.eye(4)

		x0 = np.array([0, 0, 0, 0])
		xd = np.array([0, np.pi, 0, 0])
		u0 = np.array([0])
		ud = np.array([0])

		if algorithm == 'ilqr':
			x_trajectories, u_trajectories = ilqr(cartpole_dynamics, linear_cartpole_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'ilqg':
			x_trajectories, u_trajectories = ilqg(cartpole_dynamics, linear_cartpole_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'aico':
			x_trajectories, u_trajectories = aico(cartpole_dynamics, linear_cartpole_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)


		plt.figure(1)
		legend1 = []
		for i, x_traj in enumerate(x_trajectories):
			plt.plot(x_traj[0, :], x_traj[2, :])
			legend1.append('iter ' + str(i))
		plt.xlabel('x (m)')
		plt.ylabel('xdot (m/s)')
		# plt.legend(legend1)

		plt.figure(2)
		legend2 = []
		for i, x_traj in enumerate(x_trajectories):
			plt.plot(x_traj[1, :], x_traj[3, :])
			legend2.append('iter ' + str(i))
		plt.xlabel('theta (rad)')
		plt.ylabel('thetadot (rad/s)')
		# plt.legend(legend2)

		plt.figure(3)
		legend3 = []
		for i, u_traj in enumerate(u_trajectories):
			plt.plot(time, u_traj[0, :])
			legend3.append('iter ' + str(i))
		plt.xlabel('time (s)')
		plt.ylabel('control')
		# plt.legend(legend3)

		x_traj = x_trajectories[-1]

		fig = plt.figure(4)
		fig.set_dpi(100)
		fig.set_size_inches(7, 6.5)

		ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
		base = plt.Rectangle((-1, -1), 2, 2, fc='r')
		mass = plt.Circle((0, -5), 1, fc='b')
		beam = plt.Line2D((0, 0), (0, -5), lw=2.5)

		def init():
			ax.add_patch(base)
			ax.add_patch(mass)
			ax.add_line(beam)
			return base, mass, beam

		def animate(i):
			x, y = get_cartpole_pendulum_position(x_traj[:, i])
			mass.center = (x, y)
			base.set_x(x_traj[0, i] - 1)
			beam.set_xdata([x_traj[0, i], x])
			beam.set_ydata([0, y])
			return base, mass, beam

		anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=num_steps, 
                               interval=60,
                               blit=True)

		plt.show()

	elif system == 'quadrotor':

		max_num_iter = 100
		num_steps = 100
		dt = 0.1

		T = num_steps * dt
		time = np.linspace(0, T-dt, num_steps)

		Q = np.eye(4)
		R = np.eye(2)
		Qf = np.eye(4)

		x0 = np.array([0, 0, 0, 0])
		x = np.zeros((x0.shape[0], num_steps))
		u = np.ones((2, num_steps))
		u[:, :int(num_steps/2)]*=0.2
		u[:, int(num_steps/2):]*=-0.2
		for i in range(num_steps-1):
			x[:, i+1] = x[:, i] + dt * single_quadrotor_dynamics(x[:, i], u[:, i])
		# print(x[:, -1])

		x0 = np.array([0, 0, 0, 0])
		# xd = np.array([5, 5, 0, 0])
		xd = x
		u0 = np.array([0, 0])
		# ud = np.array([0, 0])
		ud = u

		if algorithm == 'ilqr':
			x_trajectories, u_trajectories = ilqr(single_quadrotor_dynamics, linear_single_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'ilqg':
			x_trajectories, u_trajectories = ilqg(single_quadrotor_dynamics, discrete_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'aico':
			x_trajectories, u_trajectories = aico(single_quadrotor_dynamics, discrete_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)

		x = x_trajectories[-1]

		plt.figure(1)
		legend1 = []
		for i, x_traj in enumerate(x_trajectories):
			plt.plot(x_traj[0, :], x_traj[1, :])
			legend1.append('iter ' + str(i))
		plt.xlabel('x (m)')
		plt.ylabel('y (m)')
		# plt.legend(legend1)

		fig = plt.figure(2)
		fig.set_dpi(100)
		fig.set_size_inches(7, 6.5)

		ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
		quad1 = plt.Circle((x[0, 0], x[1, 0]), 1, fc='b')

		def init():
			plt.plot(x[0, :], x[1, :])
			ax.add_patch(quad1)
			return quad1,

		def animate(i):
			quad1.center = (x[0, i], x[1, i])
			return quad1,

		anim = animation.FuncAnimation(fig, animate, 
	                           init_func=init, 
	                           frames=num_steps, 
	                           interval=60,
	                           blit=True)

		plt.show()

	elif system == 'drone rope':

		max_num_iter = 100
		num_steps = 100
		dt = 0.1

		T = num_steps * dt
		time = np.linspace(0, T-dt, num_steps)

		Q = np.eye(12)
		R = np.eye(6)
		Qf = np.eye(12)

		x0 = np.array([-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])
		xd = np.array([0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
		u0 = np.array([0, 0, 0, 0, 0, 0])
		ud = np.array([0, 0, 0, 0, 0, 0])

		if algorithm == 'ilqr':
			x_trajectories, u_trajectories = ilqr(quadrotor_dynamics, linear_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'ilqg':
			x_trajectories, u_trajectories = ilqg(quadrotor_dynamics, linear_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)
		elif algorithm == 'aico':
			x_trajectories, u_trajectories = aico(quadrotor_dynamics, linear_quadrotor_dynamics, max_num_iter, num_steps, dt, x0, xd, u0, ud, Q, R, Qf)

		x = x_trajectories[-1]

		fig = plt.figure(1)
		fig.set_dpi(100)
		fig.set_size_inches(7, 6.5)

		ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
		quad1 = plt.Circle((x[0, 0], x[1, 0]), 1, fc='b')
		spring1 = plt.Line2D((x[0, 0], x[2, 0]), (x[1, 0], x[3, 0]), lw=2.5)
		quad2 = plt.Circle((x[2, 0], x[3, 0]), 1, fc='r')
		spring2 = plt.Line2D((x[2, 0], x[4, 0]), (x[3, 0], x[5, 0]), lw=2.5)
		quad3 = plt.Circle((x[4, 0], x[5, 0]), 1, fc='g')

		def init():
			ax.add_patch(quad1)
			ax.add_patch(quad2)
			ax.add_patch(quad3)
			ax.add_line(spring1)
			ax.add_line(spring2)
			return quad1, spring1, quad2, spring2, quad3

		def animate(i):
			quad1.center = (x[0, i], x[1, i])
			quad2.center = (x[2, i], x[3, i])
			quad3.center = (x[4, i], x[5, i])
			spring1.set_xdata((x[0, i], x[2, i]))
			spring1.set_ydata((x[1, i], x[3, i]))
			spring2.set_xdata((x[2, i], x[4, i]))
			spring2.set_ydata((x[3, i], x[5, i]))
			return quad1, spring1, quad2, spring2, quad3

		anim = animation.FuncAnimation(fig, animate, 
	                           init_func=init, 
	                           frames=num_steps, 
	                           interval=60,
	                           blit=True)

		plt.show()

	anim.save(system + "_" + algorithm + ".gif")
