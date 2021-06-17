import numpy as np
from matplotlib import pyplot as plt

def linear_pendulum_dynamics(theta, thetadot, b):
	A = np.array([[0, 1], [-np.cos(theta), -b]])
	B = np.array([[0], [1]])
	return A, B

def pendulum_dynamics(x, u):
	xdot = np.array([[x[1]], [u - (x[0]**2 - 1) * x[1] - np.sin(x[0])]])
	return xdot

num_iter = 100
num_steps = 100
dt = 0.1
b = 0.3

Q = np.eye(2)
R = np.eye(1)
Qf = 3*np.eye(2)

x0 = np.array([[0], [0]])
xd = np.array([[np.pi], [0]])
u0 = 0
ud = 0

theta = np.zeros((num_iter+1, num_steps))
thetadot = np.zeros((num_iter+1, num_steps))
u = np.zeros((num_iter+1, num_steps))
theta[0, 0] = x0[0]
thetadot[0, 0] = x0[1]
u[0, 0] = u0

error = np.inf
it = 1

# for it in range(1, num_iter):
while error > 1e-6 and it < num_iter-1:

# compute improved nominal control sequence by linearizing dynamics about nominal trajectory and performing Ricatti recursion

	dx = xd-x0
	du = ud-u0

	S2 = [Qf]
	S1 = [-2 * Qf @ dx]
	S0 = [dx.T @ Qf @ dx]
	
	for i in range(num_steps):
	    t = num_steps-i-1
	    At, Bt = linear_pendulum_dynamics(theta[it, t], thetadot[it, t], b)
	    # At, Bt = linear_pendulum_dynamics(theta[it-1, t], thetadot[it-1, t], b)
	    S2t = S2[i]
	    S1t = S1[i]
	    S0t = S0[i]
	    S2dot = Q - S2t @ Bt @ np.linalg.inv(R) @ Bt.T @ S2t + S2t @ At + At.T @ S2t
	    S1dot = (At.T - S2t @ Bt @ np.linalg.inv(R) @ Bt.T) @ S1t
	    S0dot = 0.25 * S1t.T @ Bt @ np.linalg.inv(R) @ Bt.T @ S1t
	    S2.append(S2t + dt * S2dot)
	    S1.append(S1t + dt * S1dot)
	    S0.append(S0t + dt * S0dot)

	S2 = S2[::-1]
	S1 = S1[::-1]
	S0 = S0[::-1]

	# compute improved nominal trajectory by forward simulating improved nominal control sequence

	for t in range(num_steps-1):
	    [At, Bt] = linear_pendulum_dynamics(theta[it+1, t], thetadot[it+1, t], b)
	    u[it+1, t] = -np.linalg.inv(R) @ Bt.T @ (S2[t] @ np.array([[theta[it+1, t] - theta[it, t]], [thetadot[it+1, t] - thetadot[it, t]]]) + 0.5 * S1[t])
	    xdot = pendulum_dynamics(np.array([theta[it+1, t], thetadot[it+1, t]]), u[it+1, t])
	    theta[it+1, t+1] = theta[it+1, t] + dt * xdot[0]
	    thetadot[it+1, t+1] = thetadot[it+1, t] + dt * xdot[1]

	    # [At, Bt] = linear_pendulum_dynamics(theta[it, t], thetadot[it, t], b)
	    # u[it, t] = -np.linalg.inv(R) @ Bt.T @ (S2[t] @ np.array([[theta[it, t] - theta[it-1, t]], [thetadot[it, t] - thetadot[it-1, t]]]) + 0.5 * S1[t])
	    # xdot = pendulum_dynamics(np.array([theta[it, t], thetadot[it, t]]), u[it, t])
	    # theta[it, t+1] = theta[it, t] + dt * xdot[0]
	    # thetadot[it, t+1] = thetadot[it, t] + dt * xdot[1]

	error = np.linalg.norm(np.array([theta[it+1, :] - theta[it, :], thetadot[it+1, :] - thetadot[it, :]]))
	# error = np.linalg.norm(np.array([theta[it, :] - theta[it-1, :], thetadot[it, :] - thetadot[it-1, :]]))
	
	# error += np.linalg.norm(xd - np.array([[theta[it, -1]], [thetadot[it, -1]]]))
	# error += np.linalg.norm(xd - np.array([[theta[it-1, -1]], [thetadot[it-1, -1]]]))
	# print(np.linalg.norm(xd - np.array([[theta[it, -1]], [thetadot[it, -1]]])))
	# print(np.linalg.norm(xd - np.array([[theta[it-1, -1]], [thetadot[it-1, -1]]])))
	it += 1

print(it)
print(np.linalg.norm(np.array([theta[it+1, :], thetadot[it+1, :]])))
print(np.linalg.norm(np.array([theta[it, :], thetadot[it, :]])))
print(np.linalg.norm(np.array([theta[it-1, :], thetadot[it-1, :]])))

T = num_steps * dt
time = np.linspace(0, T-dt, num_steps)

plt.figure(1)
legend1 = []
# for i in range(num_iter+1):
for i in range(it):
	plt.plot(theta[i, :], thetadot[i, :])
	legend1.append('iter ' + str(i))
plt.xlabel('theta (rad)')
plt.ylabel('thetadot (rad/s)')
# plt.legend(legend1)

plt.figure(2)
legend2 = []
# for i in range(num_iter+1):
for i in range(it):
	plt.plot(time, u[i, :])
	legend2.append('iter ' + str(i))
plt.xlabel('time (s)')
plt.ylabel('control')
# plt.legend(legend2)

plt.show()