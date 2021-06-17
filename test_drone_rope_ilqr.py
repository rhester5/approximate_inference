import sys
sys.path.append('../ilqr/ilqr')

'''
not outside of my capabilities to get ilqr/ddp working by final report
make it educational, tell the story, exceptional presentation skills, very thoughtful guy
explain difference between ilqr and ddp
ilqr and dpp complicated but Ricatti recursion is most complicated theoretically, rest is just implementation details
just requires you to put in the time to debug
understand you lost your partner, know you've been working hard
everyone else has put in the work to get something working, so not unreasonable to expect that you do too
not about getting a piece of software working, would rather you do your own implementation
'''

from ilqr import iLQR

import numpy as np
import matplotlib.pyplot as plt

from Simulate import simulate_system
from EKF import EKF
from iLQRNonlinearDynamics import iLQR_nonlinear_dynamics
from iLQRCost import iLQR_cost
from NonlinearModel import create_nonlinear_model_parameters
from MotionModel import MotionModel
from MeasurementModel import MeasurementModel
from iLQRController import iLQR_Controller
from Utils import actuator_limits, state_limits, hover_traj, step_response, circle_traj, line_traj

if __name__ == '__main__':

	# TODO
	# - copy needed files over
	# - create initial and desired trajectories (really just initial = all zeros and desired = just the goal vector)
	# - make iLQR cost for drone rope
	# - make simulate for drone rope
	# - change all matrix sizes for drone rope
	# LATER
	# - generalize to n drones and n-1 springs
	# - generalize to arbitrary dynamics

	np.random.seed(21)

	T = 0.05 # seconds
	time = 10.0 # seconds
	K = int(time/T) # num time steps
	init_P = 0.0

	noise = True
	filt = True
	if filt and not noise:
		raise TypeError('noise must be True to use Kalman Filter')

	controller = 'iLQR'

	reoptimize_traj = True
	reoptimize_freq = 2.0 # Hz
	reoptimize_k = int(1.0/T/reoptimize_freq)

	trajectory = hover_traj(0.0, 0.0, 5.0, K, linear)
	trajectory = step_response(trajectory, 0.0, 5.0, int(2.5/T), int(time/T))
	# trajectory = circle_traj(1.0, 1.0, K, linear)
	# trajectory = line_traj('x', 5.0, 0.0, 60.0, K, linear)

	x0 = trajectory[0, :]
		
	if noise:
		model_cov = [1e-3, 1e-3, 1e-5, 1e-5]
		meas_cov = [1e-1, 1e-1, 1e-3, 1e-3]
	else:
		model_cov = [0, 0, 0, 0]
		meas_cov = [0, 0, 0, 0]
	f, F, h, H, Q, R = create_nonlinear_model_parameters(T, mass, g, J, d, model_cov, meas_cov)
	motion_model = MotionModel([f, Q])
	meas_model = MeasurementModel(h, R, linear)
	if filt:
		P0 = init_P*np.eye(18)
		kalman_filter = EKF(f, F, h, H, Q, R, x0, P0)
		est_state = np.zeros((K, 18))
		est_cov = np.zeros((K, 18, 18))
		error = np.zeros((K, 18))

	if controller == 'iLQR':
		dynamics = iLQR_nonlinear_dynamics(T, mass, g, J)
		us_init = np.zeros((K, 4))
		us_init[:, 0] = us_init[:, 0] + g
		xs = np.zeros((K, 18))
		us = np.zeros((K, 4))

		# chunks = 10
		# num_steps = int(K/chunks)
		# for i in range(chunks):
		# 	print('chunk ', i)
		# 	x0_temp = trajectory[i*num_steps, :]
		# 	if i == chunks-1:
		# 		goal = trajectory[-1, 0:3]
		# 	else:
		# 		goal = trajectory[(i+1)*num_steps, 0:3]
		# 	cost = iLQR_cost(goal, linear)
		# 	print('run iLQR')
		# 	ilqr = iLQR(dynamics, cost, num_steps-1)
		# 	xs_chunk, us_chunk = ilqr.fit(x0_temp, us_init[i*num_steps:(i+1)*num_steps-1, :])
		# 	xs[i*num_steps:(i+1)*num_steps, :] = xs_chunk
		# 	us[i*num_steps:(i+1)*num_steps-1, :] = us_chunk

		goal = trajectory[-1, 0:3]
		cost = iLQR_cost(goal, linear)
		print('run iLQR')
		ilqr = iLQR(dynamics, cost, K)
		xs, us = ilqr.fit(x0, us_init)

		trajectory = xs
		if controller == 'iLQR':
			# controller = iLQR_Controller(us)
			controller = iLQR_Controller(dynamics, us, cost, us_init, K)
			# controller = iLQR_Controller(dynamics, trajectory, us, linear, 20)

	state, meas, setpoint, control = simulate_system(K, x0, trajectory, controller, motion_model, meas_model)
	tracking_error = np.sum(np.abs(setpoint-state))/K
	measurement_error = np.sum(np.abs(setpoint-meas))/K
	print(tracking_error)
	print(measurement_error)
	if filt:
		x = x0
		P = P0
		for k in range(K):
			kp1 = k+1 if k < K-1 else k
			u = controller(x, trajectory[k], trajectory[kp1], k)

			kalman_filter.predict(u)
			kalman_filter.update(meas[k, :])
			(x, P) = kalman_filter.get_state()
			est_state[k, :] = x
			est_cov[k, ...] = P
			error[k, :] = x - state[k, :]
		filtered_tracking_error = np.sum(np.abs(setpoint-est_state))/K
		print(filtered_tracking_error)


	t = np.linspace(0, K-1, K)

	plt.figure(1)
	if filt: # and linear:
		plt.plot(meas[:, 0], meas[:, 1], 'rx')
	plt.plot(setpoint[:, 0], setpoint[:, 1], '-go')
	plt.plot(state[:, 0], state[:, 1], '-bo')
	if filt:
		plt.plot(est_state[:, 0], est_state[:, 1], '-ko')
	plt.xlabel('x [m]')
	plt.ylabel('y [m]')
	if filt: # and linear:
		plt.legend(['observed measurement', 'setpoint', 'ground truth', 'state estimate'])
	# elif filt and not linear:
	# 	plt.legend(['setpoint', 'ground truth', 'state estimate'])
	else:
		plt.legend(['setpoint', 'ground truth'])

	plt.figure(2)
	if filt: # and linear:
		plt.plot(t, meas[:, 2], 'rx')
	# elif filt and not linear:
	# 	plt.plot(t, meas[:, 0], 'rx')
	plt.plot(t, setpoint[:, 2], '-go')
	plt.plot(t, state[:, 2], '-bo')
	if filt:
		plt.plot(t, est_state[:, 2], '-ko')
	plt.xlabel('t [s]')
	plt.ylabel('z [m]')
	if filt:
		plt.legend(['observed measurement', 'setpoint', 'ground truth', 'state estimate'])
	else:
		plt.legend(['setpoint', 'ground truth'])

	plt.figure(3)
	if filt: # and linear:
		plt.plot(t, meas[:, 0], 'rx')
	plt.plot(t, setpoint[:, 0], '-go')
	plt.plot(t, state[:, 0], '-bo')
	if filt:
		plt.plot(t, est_state[:, 0], '-ko')
	plt.xlabel('t [s]')
	plt.ylabel('x [m]')
	if filt: # and linear:
		plt.legend(['observed measurement', 'setpoint', 'ground truth', 'state estimate'])
	# elif filt and not linear:
	# 	plt.legend(['setpoint', 'ground truth', 'state estimate'])
	else:
		plt.legend(['setpoint', 'ground truth'])
	

	# plt.figure(4)
	# if filt and linear:
	# 	plt.plot(t, meas[:, 1], 'rx')
	# plt.plot(t, setpoint[:, 1], '-go')
	# plt.plot(t, state[:, 1], '-bo')
	# if filt:
	# 	plt.plot(t, est_state[:, 1], '-ko')
	# plt.xlabel('t [s]')
	# plt.ylabel('y [m]')
	# if filt and linear:
	# 	plt.legend(['observed measurement', 'setpoint', 'ground truth', 'state estimate'])
	# elif filt and not linear:
	# 	plt.legend(['setpoint', 'ground truth', 'state estimate'])
	# else:
	# 	plt.legend(['setpoint', 'ground truth'])

	# plt.figure(5)
	# plt.plot(t, np.log(state[:, 6]))
	# plt.plot(t, np.log(state[:, 7]))
	# plt.plot(t, np.log(state[:, 8]))
	# plt.plot(t, np.log(state[:, 9]))
	# plt.plot(t, np.log(state[:, 10]))
	# plt.plot(t, np.log(state[:, 11]))
	# plt.plot(t, np.log(state[:, 12]))
	# plt.plot(t, np.log(state[:, 13]))
	# plt.plot(t, np.log(state[:, 14]))
	# plt.legend(['Rx1', 'Rx2', 'Rx3', 'Ry1', 'Ry2', 'Ry3', 'Rz1', 'Rz2', 'Rz3'])

	# plt.figure(6)
	# plt.plot(t, state[:, 3])
	# plt.plot(t, state[:, 4])
	# plt.plot(t, state[:, 5])
	# plt.legend(['vx', 'vy', 'vz'])

	# plt.figure(7)
	# plt.plot(t, state[:, 15])
	# plt.plot(t, state[:, 16])
	# plt.plot(t, state[:, 17])
	# plt.legend(['wx', 'wy', 'wz'])

	# plt.figure(8)
	# plt.plot(t, control[:, 0])
	# plt.plot(t, control[:, 1])
	# plt.plot(t, control[:, 2])
	# plt.plot(t, control[:, 3])
	# plt.legend(['xdd', 'ydd', 'zdd'])
	# plt.legend(['f', 'M1', 'M2', 'M3'])

	plt.show()