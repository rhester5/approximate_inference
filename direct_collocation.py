import numpy as np
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
from matplotlib import pyplot as plt
from matplotlib import animation

from ilqr import get_pendulum_position

def pendulum_dynamics(x, u):
	xdot = np.array([x[1], u - (x[0]**2 - 1) * x[1] - np.sin(x[0])])
	return xdot

# def direct_collocation(dynamics, x0, xd, u0, ud, Q, R, num_steps, dt):
# 	X0 = np.zeros((1, (x0.shape[0] + u0.shape[0]) * num_steps))

# 	def cost_function(X):
# 		x = np.zeros((x0.shape[0], num_steps))
# 		for i in range(x0.shape[0]):
# 			x[i, :] = X[i*num_steps:(i+1)*num_steps]
# 		u = np.zeros((u0.shape[0], num_steps))
# 		for i in range(x0.shape[0], x0.shape[0]+u0.shape[0]):
# 			u[i-x0.shape[0], :] = X[i*num_steps:(i+1)*num_steps]
# 		cost = 0
# 		for i in range(num_steps):
# 			cost += x[:, i].T @ Q @ x[:, i] + u[:, i].T @ R @ u[:, i]
# 		return cost

# 	cons = []



# 	X = minimize(cost_function, X0, constraints=cons)

# Constrain x[0] < sin(x[1]) + 1.9
# con = lambda x: x[0] - np.sin(x[1])
# nlc = NonlinearConstraint(con, -np.inf, 1.9)

# bnds = ((0, None), (0, None)) # bounds=bnds in solver


'''
The constraint has the general inequality form:

lb <= fun(x) <= ub
Here the vector of independent variables x is passed as ndarray of shape 
(n,) and fun returns a vector with m components.

It is possible to use equal bounds to represent an equality constraint or 
infinite bounds to represent a one-sided constraint.
'''

def shooting(dynamics, x0, xd, u0, ud, num_steps, dt):
	X0 = np.zeros((u0.shape[0]*num_steps,))

	def cost_function(X):
		u = np.zeros((u0.shape[0], num_steps))
		for i in range(u0.shape[0]):
			u[i, :] = X[i*num_steps:(i+1)*num_steps]
		cost = 0.0
		for i in range(num_steps):
			cost += u[:, i].T @ u[:, i]
		return cost

	cons = []
	for i in range(x0.shape[0]):
		def shoot(X):
			x = x0
			for j in range(num_steps):
				x += dt * dynamics(x, X[j])
			return x[i] - xd[i]
		cons.append(NonlinearConstraint(shoot, 0.0, 0.0))

	X = minimize(cost_function, X0, constraints=cons, tol=0.00000001, options={'maxiter':10000, 'disp':True})

	return X

if __name__ == '__main__':
	num_steps = 100
	dt = 0.1

	T = num_steps * dt
	time = np.linspace(0, T-dt, num_steps)

	x0 = np.array([0.0, 0.0])
	xd = np.array([np.pi, 0.0])
	u0 = np.array([0.0])
	ud = np.array([0.0])

	X = shooting(pendulum_dynamics, x0, xd, u0, ud, num_steps, dt)
	X = X.x
	u_traj = np.zeros((u0.shape[0], num_steps))
	for i in range(u0.shape[0]):
		u_traj[i, :] = X[i*num_steps:(i+1)*num_steps]
	x_traj = np.zeros((x0.shape[0], num_steps))
	x_traj[:, 0] = x0
	for i in range(num_steps-1):
		x_traj[:, i+1] = x_traj[:, i] + dt * pendulum_dynamics(x_traj[:, i], u_traj[:, i])

	plt.figure(1)
	plt.plot(x_traj[0, :], x_traj[1, :])
	plt.xlabel('theta (rad)')
	plt.ylabel('thetadot (rad/s)')

	plt.figure(2)
	plt.plot(time, u_traj[0, :])
	plt.xlabel('time (s)')
	plt.ylabel('control')

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