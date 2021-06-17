import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def get_distance(x1, y1, x2, y2):
	return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def get_angle(x1, y1, x2, y2):
	return np.arctan2(y2-y1, x2-x1)

def x_derivs(a, c, b, d, k, l):

	# a_deriv = (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2) - (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2)
	# b_deriv = (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2)
	# c_deriv = -(k * (c - d)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) + (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) + (k * (c - d)**2)/((a - b)**2 + (c - d)**2)
	# d_deriv = (k * (c - d)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) - (k * (c - d)**2)/((a - b)**2 + (c - d)**2)

	a_deriv = -(k * (a - b)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) + (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) + (k * (a - b)**2)/((a - b)**2 + (c - d)**2)
	b_deriv = (k * (a - b)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) - (k * (a - b)**2)/((a - b)**2 + (c - d)**2)
	c_deriv = (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2) - (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2)
	d_deriv = (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2)

	return [-a_deriv, -c_deriv, -b_deriv, -d_deriv]
	# return [a_deriv, c_deriv, d_deriv, d_deriv]
	# return [-a_deriv, -b_deriv, -c_deriv, -d_deriv]
	# return [a_deriv, b_deriv, c_deriv, d_deriv]

def y_derivs(a, c, b, d, k, l):

	# a_deriv = -(k * (a - b)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) + (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) + (k * (a - b)**2)/((a - b)**2 + (c - d)**2)
	# b_deriv = (k * (a - b)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) - (k * (a - b)**2)/((a - b)**2 + (c - d)**2)
	# c_deriv = (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2) - (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2)
	# d_deriv = (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2)

	a_deriv = (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2) - (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2)
	b_deriv = (k * (a - b) * (c - d) * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (a - b) * (c - d))/((a - b)**2 + (c - d)**2)
	c_deriv = -(k * (c - d)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) + (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) + (k * (c - d)**2)/((a - b)**2 + (c - d)**2)
	d_deriv = (k * (c - d)**2 * (np.sqrt((a - b)**2 + (c - d)**2) - l))/((a - b)**2 + (c - d)**2)**(3/2) - (k * (np.sqrt((a - b)**2 + (c - d)**2) - l))/np.sqrt((a - b)**2 + (c - d)**2) - (k * (c - d)**2)/((a - b)**2 + (c - d)**2)

	return [-a_deriv, -c_deriv, -b_deriv, -d_deriv]
	# return [a_deriv, c_deriv, d_deriv, d_deriv]
	# return [-a_deriv, -b_deriv, -c_deriv, -d_deriv]
	# return [a_deriv, b_deriv, c_deriv, d_deriv]

def single_quadrotor_dynamics(x, u):
	xdot = np.zeros((4,))
	x1, y1, x1d, y1d = x
	x1dd, y1dd = u
	xdot[0] = x1d
	xdot[1] = y1d
	xdot[2] = x1dd
	xdot[3] = y1dd
	return xdot

def linear_single_quadrotor_dynamics(x, u):
	A = np.zeros((x.shape[0], x.shape[0]))
	B = np.zeros((x.shape[0], u.shape[0]))
	A[0:2, 2:] = np.eye(2)
	B[2:, :] = np.eye(2)
	return A, B

def discrete_quadrotor_dynamics(x, u, dt):
	A = np.zeros((x.shape[0], x.shape[0]))
	B = np.zeros((x.shape[0], u.shape[0]))
	A[0:2, 2:] = np.eye(2)
	B[2:, :] = np.eye(2)
	return np.eye(A.shape[0]) + dt*A, dt*B


def quadrotor_dynamics(x, u):
	k = 1 # spring constant
	rl = 1 # rest_length
	xdot = np.zeros((12,))
	x1, y1, x2, y2, x3, y3, x1d, y1d, x2d, y2d, x3d, y3d = x
	x1dd, y1dd, x2dd, y2dd, x3dd, y3dd = u
	xdot[0] = x1d
	xdot[1] = y1d
	xdot[2] = x2d
	xdot[3] = y2d
	xdot[4] = x3d
	xdot[5] = y3d
	dist21 = get_distance(x2, y2, x1, y1)
	dist12 = get_distance(x1, y1, x2, y2)
	dist32 = get_distance(x3, y3, x2, y2)
	dist23 = get_distance(x2, y2, x3, y3)
	angle21 = get_angle(x2, y2, x1, y1)
	angle12 = get_angle(x1, y1, x2, y2)
	angle32 = get_angle(x3, y3, x2, y2)
	angle23 = get_angle(x2, y2, x3, y3)
	xdot[6] = x1dd - k * (dist21 - rl) * np.cos(angle21)
	xdot[7] = y1dd - k * (dist21 - rl) * np.sin(angle21)
	xdot[8] = x2dd - k * (dist12 - rl) * np.cos(angle12) - k * (dist32 - rl) * np.cos(angle32)
	xdot[9] = y2dd - k * (dist12 - rl) * np.sin(angle12) - k * (dist32 - rl) * np.sin(angle32)
	xdot[10] = x3dd - k * (dist23 - rl) * np.cos(angle23)
	xdot[11] = y3dd - k * (dist23 - rl) * np.sin(angle23)
	return xdot

def linear_quadrotor_dynamics(x, u):
	# TODO fml this does not work
	k = 1
	rl = 1
	A = np.zeros((x.shape[0], x.shape[0]))
	B = np.zeros((x.shape[0], u.shape[0]))
	A[0:6, 6:] = np.eye(6)
	x1, y1, x2, y2, x3, y3, x1d, y1d, x2d, y2d, x3d, y3d = x
	x1dd, y1dd, x2dd, y2dd, x3dd, y3dd = u
	x21_derivs = x_derivs(x2, y2, x1, y1, k, rl)
	y21_derivs = y_derivs(x2, y2, x1, y1, k, rl)
	x12_derivs = x_derivs(x1, y1, x2, y2, k, rl)
	y12_derivs = y_derivs(x1, y1, x2, y2, k, rl)
	x32_derivs = x_derivs(x3, y3, x2, y2, k, rl)
	y32_derivs = y_derivs(x3, y3, x2, y2, k, rl)
	x23_derivs = x_derivs(x2, y2, x3, y3, k, rl)
	y23_derivs = y_derivs(x2, y2, x3, y3, k, rl)
	A[6, 0:4] = x21_derivs
	A[7, 0:4] = y21_derivs
	A[8, 0:4] = x12_derivs
	A[9, 0:4] = y12_derivs
	A[8, 2:6] += x32_derivs
	A[9, 2:6] += y32_derivs
	A[10, 2:6] = x23_derivs
	A[11, 2:6] = y23_derivs
	B[6:, :] = np.eye(6)
	return A, B

def simulate(dt, num_steps, x0, u):
	x = np.zeros((x0.shape[0], num_steps))
	x[:, 0] = x0
	for i in range(1, num_steps):
		x[:, i] = x[:, i-1] + dt * quadrotor_dynamics(x[:, i-1], u)
		# A, B = linear_quadrotor_dynamics(x[:, i-1], u)
		# x[:, i] = x[:, i-1] + dt * (A @ x[:, i-1] + B @ u)
	return x
	
if __name__ == '__main__':
	x0 = np.array([1, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0])
	u = np.array([0, 0, -1, 1, 0, 0])
	dt = 0.1
	num_steps = 100
	x = simulate(dt, num_steps, x0, u)

	x0 = np.array([0, 0, 0, 0])
	x = np.zeros((x0.shape[0], num_steps))
	u = np.ones((2, num_steps))
	u[:, :int(num_steps/2)]*=0.2
	u[:, int(num_steps/2):]*=-0.2
	for i in range(num_steps-1):
		x[:, i+1] = x[:, i] + dt * single_quadrotor_dynamics(x[:, i], u[:, i])
	print(x[:, -1])

	fig = plt.figure(1)
	fig.set_dpi(100)
	fig.set_size_inches(7, 6.5)

	ax = plt.axes(xlim=(-10, 10), ylim=(-10, 10))
	quad1 = plt.Circle((x[0, 0], x[1, 0]), 1, fc='b')
	# spring1 = plt.Line2D((x[0, 0], x[2, 0]), (x[1, 0], x[3, 0]), lw=2.5)
	# quad2 = plt.Circle((x[2, 0], x[3, 0]), 1, fc='r')
	# spring2 = plt.Line2D((x[2, 0], x[4, 0]), (x[3, 0], x[5, 0]), lw=2.5)
	# quad3 = plt.Circle((x[4, 0], x[5, 0]), 1, fc='g')

	def init():
		ax.add_patch(quad1)
		# ax.add_patch(quad2)
		# ax.add_patch(quad3)
		# ax.add_line(spring1)
		# ax.add_line(spring2)
		return quad1,# spring1, quad2, spring2, quad3

	def animate(i):
		quad1.center = (x[0, i], x[1, i])
		# quad2.center = (x[2, i], x[3, i])
		# quad3.center = (x[4, i], x[5, i])
		# spring1.set_xdata((x[0, i], x[2, i]))
		# spring1.set_ydata((x[1, i], x[3, i]))
		# spring2.set_xdata((x[2, i], x[4, i]))
		# spring2.set_ydata((x[3, i], x[5, i]))
		return quad1,# spring1, quad2, spring2, quad3

	anim = animation.FuncAnimation(fig, animate, 
                           init_func=init, 
                           frames=num_steps, 
                           interval=60,
                           blit=True)

	# anim.save('quad_spring.gif')

	plt.show()
