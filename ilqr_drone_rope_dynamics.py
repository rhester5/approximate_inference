import sys
sys.path.append('../ilqr/ilqr')

from ilqr import iLQR 
import theano.tensor as T
from ilqr.dynamics import AutoDiffDynamics
import numpy as np

def get_distance(x1, y1, x2, y2):
	return np.sqrt((x2-x1)**2 + (y2-1)**2)

def get_angle(x1, y1, x2, y2):
	return np.arctan2(y2-y1, x2-x1)

def iLQR_nonlinear_dynamics(dt):
	# spring constant
	k = 1

	# rest length
	rl = 1

	# state
	x1 = T.dscalar("x1")  # Position
	y1 = T.dscalar("y1")  # Position
	x2 = T.dscalar("x2")  # Position
	y2 = T.dscalar("y2")  # Position
	x3 = T.dscalar("x3")  # Position
	y3 = T.dscalar("y3")  # Position
	x1dot = T.dscalar("x1_dot")  # Velocity
	y1dot = T.dscalar("y1_dot")  # Velocity
	x2dot = T.dscalar("x2_dot")  # Velocity
	y2dot = T.dscalar("y2_dot")  # Velocity
	x3dot = T.dscalar("x3_dot")  # Velocity
	y3dot = T.dscalar("y3_dot")  # Velocity

	# control
	x1ddot = T.dscalar("x1_ddot")  # Acceleration
	y1ddot = T.dscalar("y1_ddot")  # Acceleration
	x2ddot = T.dscalar("x2_ddot")  # Acceleration
	y2ddot = T.dscalar("y2_ddot")  # Acceleration
	x3ddot = T.dscalar("x3_ddot")  # Acceleration
	y3ddot = T.dscalar("y3_ddot")  # Acceleration

	# compute distances and angles
	dist21 = get_distance(x2, y2, x1, y1)
	dist12 = get_distance(x1, y1, x2, y2)
	dist32 = get_distance(x3, y3, x2, y2)
	dist23 = get_distance(x2, y2, x3, y3)
	angle21 = get_angle(x2, y2, x1, y1)
	angle12 = get_angle(x1, y1, x2, y2)
	angle32 = get_angle(x3, y3, x2, y2)
	angle23 = get_angle(x2, y2, x3, y3)

	# Discrete dynamics model definition
	f = T.stack([
	    x1 + dt * x1dot,
	    y1 + dt * y1dot,
	    x2 + dt * x2dot,
	    y2 + dt * y2dot,
	    x3 + dt * x3dot,
	    y3 + dt * y3dot,
	    x1dot + dt * (x1dd - k * (dist21 - rl) * np.cos(angle21))
	    y1dot + dt * (y1dd - k * (dist21 - rl) * np.sin(angle21))
	    x2dot + dt * (x2dd - k * (dist12 - rl) * np.cos(angle12) - k * (dist32 - rl) * np.cos(angle32))
	    y2dot + dt * (y2dd - k * (dist12 - rl) * np.sin(angle12) - k * (dist32 - rl) * np.sin(angle32))
	    x3dot + dt * (x3dd - k * (dist23 - rl) * np.cos(angle23))
	    y3dot + dt * (y3dd - k * (dist23 - rl) * np.sin(angle23))
	])

	x_inputs = [x1, y1, x2, y2, x3, y3, x1dot, y1dot, x2dot, y2dot, x3dot, y3dot]  # State vector
	u_inputs = [x1ddot, y1ddot, x2ddot, y2ddot, x3ddot, y3ddot]  # Control vector

	# Compile the dynamics
	# NOTE: This can be slow as it's computing and compiling the derivatives
	# But that's okay since it's only a one-time cost on startup
	print('compiling dynamics')
	dynamics = AutoDiffDynamics(f, x_inputs, u_inputs)
	return dynamics