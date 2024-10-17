import numpy as np
import matplotlib.pyplot as plt
import json
import GeometeryClass
import Flow
from scipy.optimize import newton


class Main:
	"""
	A class to represent the main application logic for the airfoil project.

	Attributes
	----------
	config_file : str
		The path to the input json  file.
	radius : float
		The radius of the cylinder.
	x_low_val : float
		The lower limit for the x-axis in the plot.
	x_up_val : float
		The upper limit for the x-axis in the plot.
	geometry : GeometeryClass.Geometery
		An instance of the Geometery class.

	Methods
	-------
	load_config():
		Loads the configuration from the JSON file.
	setup_geometry():
		Initializes the Geometery object and calculates the camber, upper surface, and lower surface.
	plot():
		Plots the camber line, upper surface, and lower surface.
	surface_tangent(x):
		Calculate the surface tangent vectors at a given x-coordinate.
	surface_normal(x):
		Calculate the surface normal vectors at a given x-coordinate.
	run():
		Executes the main logic of the application.
	"""

	def __init__(self, config_file):
		"""
		Constructs all the necessary attributes for the Main object.

		Parameters
		----------
		config_file : str
			The path to the configuration file.
		"""
		self.config_file = config_file
		

	def load_config(self):
		"""
		Loads the configuration from the JSON file.
		"""
		with open(self.config_file, 'r') as file:
			json_vals = json.load(file)
		self.radius = json_vals['geometry']['cylinder_radius']
		self.z0 = json_vals['geometry']['z0']

		self.free_stream_velocity = json_vals['operating']['freestream_velocity']
		self.angle_of_attack = json_vals['operating']['angle_of_attack[deg]']
		self.gamma = json_vals['operating']['vortex_strength']
		
		self.x_start = json_vals['plot']['x_start']
		self.x_low_val = json_vals['plot']['x_lower_limit']
		self.x_up_val = json_vals['plot']['x_upper_limit']
		self.delta_s = json_vals['plot']['delta_s']
		self.n_lines = json_vals['plot']['n_lines']
		self.delta_y = json_vals['plot']['delta_y']

		self.LE = -1 * self.radius+ self.z0[0]
		self.TE = self.radius + self.z0[0]
		


	def setup_Geometry(self):
		"""
		Initializes the Geometery object and calculates the camber, upper surface, and lower surface.
		"""
		self.geometry = GeometeryClass.Geometery(self.radius, self.z0)
		
	def load_flow_field(self, V_inf, alpha):
		"""
		Loads the flow field parameters.

		Parameters:
		V_inf (float): The free stream velocity.
		alpha (float): The angle of attack.
		"""
		self.flow = Flow.Flow(self.radius, V_inf, alpha, self.x_low_val, self.x_up_val, self.gamma, self.z0)

	def surface_tangent(self, x):
		"""
        Calculate the surface normal vectors at a given x-coordinate.

        Parameters:
        surface (np.ndarray): A 2xN array where the first row contains x-coordinates and the second row contains y-coordinates.
        x (float): The x-coordinate at which to calculate the surface normals.

        Returns:
        tuple: A tuple containing the unit normal vectors for the upper and lower surfaces.
        """
		delta = 1e-5  # Small perturbation for numerical differentiation

        # Interpolate to find the y-coordinate at x
		#Trailing Edge Case
		if np.abs(x-self.TE) < delta:
			_, y_upper_minus, y_lower_minus = self.geometry.circle(x-delta)
			tangent_upper = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), -1*(y_upper_minus[1]-y_lower_minus[1])])
			tangent_lower = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), -1*(y_upper_minus[1]-y_lower_minus[1])])
		#Leading Edge Case	
		elif np.abs(x-self.LE) < delta:
			_,y_upper_plus, y_lower_plus = self.geometry.circle(x+delta)

			tangent_upper = np.array([-1*np.abs(y_upper_plus[0]-y_upper_plus[0]), y_upper_plus[1]-y_lower_plus[1]])
			tangent_lower = np.array([-1*np.abs(y_upper_plus[0]-y_upper_plus[0]), y_upper_plus[1]-y_lower_plus[1]])
		#General Case
		else:
			_,y_upper_plus, y_lower_plus = self.geometry.circle(x+delta)

			_,y_upper_minus, y_lower_minus = self.geometry.circle(x-delta)

			tangent_upper = np.array([np.abs(y_upper_plus[0]-y_upper_minus[0]), y_upper_plus[1] - y_upper_minus[1]])
			tangent_lower = np.array([-1 * np.abs(y_lower_plus[0]-y_lower_minus[0]), -1*(y_lower_plus[1] - y_lower_minus[1])])

			


		unit_tangent_upper = tangent_upper / np.linalg.norm(tangent_upper)	
		unit_tangent_lower = tangent_lower / np.linalg.norm(tangent_lower)	


		return unit_tangent_upper, unit_tangent_lower

	def surface_normal(self, x):
		"""""
        Calculate the surface tangent vectors at a given x-coordinate.

        Parameters:
        surface (np.ndarray): A 2xN array where the first row contains x-coordinates and the second row contains y-coordinates.
        x (float): The x-coordinate at which to calculate the surface tangents.

        Returns:
        np.ndarray: The unit tangent vector at the given x-coordinate.
        """
		delta = 1e-5  # Small perturbation for numerical differentiation

        # Interpolate to find the y-coordinate at x
		if np.abs(x-self.radius) < delta:
			_,y_upper_minus, y_lower_minus = self.geometry.circle(x-delta)
			tangent_upper = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), y_upper_minus[1]-y_lower_minus[1]])
			tangent_lower = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), y_upper_minus[1]-y_lower_minus[1]])

			normal_upper = np.array([tangent_upper[1], tangent_upper[0]])
			normal_lower = np.array([tangent_lower[1], -tangent_lower[0]])		
		elif np.abs(x+self.radius) < delta:
			_,y_upper_plus, y_lower_plus = self.geometry.circle(x+delta)

			tangent_upper = np.array([-1 * np.abs(y_upper_minus[0]-y_upper_minus[0]), y_upper_plus[1]-y_lower_plus[1]])
			tangent_lower = np.array([-1 * np.abs(y_upper_minus[0]-y_upper_minus[0]), y_upper_plus[1]-y_lower_plus[1]])

			normal_upper = np.array([-tangent_upper[1], tangent_upper[0]])
			normal_lower = np.array([-tangent_lower[1], -tangent_lower[0]]) 
		else:
			_,y_upper_plus, y_lower_plus = self.geometry.circle(x+delta)

			_,y_upper_minus, y_lower_minus = self.geometry.circle(x-delta)

			tangent_upper = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), y_upper_plus[1] - y_upper_minus[1]])
			tangent_lower = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), (y_lower_plus[1] - y_lower_minus[1])])


			normal_upper = np.array([-tangent_upper[1], tangent_upper[0]])
			normal_lower = np.array([tangent_lower[1], -tangent_lower[0]])

		unit_normal_upper = normal_upper / np.linalg.norm(normal_upper)	
		unit_normal_lower = normal_lower / np.linalg.norm(normal_lower)

		return unit_normal_upper, unit_normal_lower
	
	def surface_tangential_velocity(self, x, upper=True):
		"""
		Calculate the tangential velocity at a given x-coordinate.

		Parameters:
		x (float): The x-coordinate at which to calculate the tangential velocity.

		Returns:
		tuple: A tuple containing the tangential velocity for the upper and lower surfaces.
		"""
		unit_tangent_upper, unit_tangent_lower = self.surface_tangent(x)

		if upper:
			_, point, _ = self.geometry.circle(x)
			velocity_upper = np.dot(self.flow.flow_over_cylinder_complex(point[0], point[1]), unit_tangent_upper)
			return velocity_upper
		else:
			_, _, point = self.geometry.circle(x)
			velocity = self.flow.flow_over_cylinder_complex(point[0], point[1])
			velocity_lower = np.dot(velocity, unit_tangent_lower)
			return velocity_lower

	
	def stagnation_point(self):
		"""
        Calculate the stagnation point of the flow field.

        Returns:
        tuple: A tuple containing the x and y coordinates of the forward and aft stagnation points.
        """
		epsilon = 1e-5  # Small perturbation for numerical differentiation

        # Check leading edge (LE) and trailing edge (TE)
		velocity_LE = self.surface_tangential_velocity(self.LE)
		velocity_TE = self.surface_tangential_velocity(self.TE, upper=False)

		if np.abs(velocity_LE) < epsilon:
			forward_stagnation_point = [self.LE, self.z0[1]]
		else:
			if velocity_LE < 0:
                # Use upper surface
				upper_flag = True
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				x_stag = self.bisection_method(velocity_function, self.LE, self.TE-self.radius)
				_, forward_stagnation_point, _ = self.geometry.circle(x_stag)
			else:
                # Use lower surface
				upper_flag = False
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				x_stag = self.bisection_method(velocity_function, self.LE, self.TE-self.radius)
				_, _, forward_stagnation_point = self.geometry.circle(x_stag)

		if np.abs(velocity_TE) < epsilon:
			aft_stagnation_point = [self.TE, self.z0[1]]
		else:
			if velocity_TE > 0:
                # Use lower surface
				upper_flag = False
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				x_stag = self.bisection_method(velocity_function, self.LE+self.radius, self.TE)
				_, _, aft_stagnation_point = self.geometry.circle(x_stag)
			else:
                # Use upper surface
				upper_flag = True
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				x_stag = self.bisection_method(velocity_function, self.LE+self.radius, self.TE)
				_, aft_stagnation_point, _ = self.geometry.circle(x_stag)

		return forward_stagnation_point, aft_stagnation_point
	
	def bisection_method(self,f, a, b, tol=1e-10, max_iter=100):
		"""
		Bisection method to find the root of a function f(x) = 0.
		
		Parameters:
		f: function
			The function for which the root is to be found.
		a: float
			Left endpoint of the interval.
		b: float
			Right endpoint of the interval.
		tol: float, optional
			Tolerance for the root. The method stops when |b - a| < tol.
		max_iter: int, optional
			Maximum number of iterations to perform.
		
		Returns:
		float
			The approximation for the root, or None if no root was found within the interval.
		"""
		
		# Check if a root is guaranteed in the initial interval
		if f(a) * f(b) >= 0:
			print(f(a), f(b))
			print("Bisection method fails. The function must have different signs at the endpoints.")
			return None
		
		for i in range(max_iter):
			# Calculate midpoint
			c = (a + b) / 2
			
			# Check if the function value at the midpoint is close enough to zero
			if abs(f(c)) < tol:
				print(f"Root found after {i+1} iterations: {f(c)}")
				return c
			
			# Update the interval based on the sign of f(c)
			if f(c) * f(a) < 0:
				b = c  # Root is in the left half
			else:
				a = c  # Root is in the right half
		
		print("Warning: Maximum number of iterations reached.")
		return None


	def plot(self, x_start, x_lower_limit, x_upper_limit, delta_s, n_lines, delta_y):
		"""
		Plot the streamlines for the flow field.

		Parameters:
		x_start (float): The x-coordinate at which to start the streamlines.
		x_lower_limit (float): The lower limit for the x-axis.
		x_upper_limit (float): The upper limit for the x-axis.
		delta_s (float): The step size for the streamlines.
		n_lines (int): The number of streamlines to plot.
		delta_y (float): The spacing between streamlines.
		"""
		x = np.linspace(self.LE, self.TE, 1000)
		camber = np.empty((2, 0))
		upper_surface = np.empty((2, 0))
		lower_surface = np.empty((2, 0))

		for i in range(len(x)):
			camber_temp, upper_surface_temp, lower_surface_temp = self.geometry.circle(x[i])

			# Append the new values to the arrays
			camber = np.hstack((camber, camber_temp.reshape(2, 1)))
			upper_surface = np.hstack((upper_surface, upper_surface_temp.reshape(2, 1)))
			lower_surface = np.hstack((lower_surface, lower_surface_temp.reshape(2, 1)))

		#calculate stagnation points
		forward_stagnation_point, aft_stagnation_point = self.stagnation_point()
		forward_normal_stag_up, forward_normal_stag_down = self.surface_normal(forward_stagnation_point[0])
		aft_normal_stag_up, aft_normal_stag_down = self.surface_normal(aft_stagnation_point[0])
		if forward_stagnation_point[1] < 0:
			forward_normal_stag = forward_normal_stag_down
		else:
			forward_normal_stag = forward_normal_stag_up
		if aft_stagnation_point[1] < 0:
			aft_normal_stag = aft_normal_stag_down
		else:
			aft_normal_stag = aft_normal_stag_up

		# Set up the plot
		plt.figure()
		plt.plot(camber[0, :], camber[1, :], label='Camber', color='red')
		plt.plot(upper_surface[0, :], upper_surface[1, :], label='Upper Surface', color='blue')
		plt.plot(lower_surface[0, :], lower_surface[1, :], label='Lower Surface',color='blue')
		plt.xlabel('X')
		plt.ylabel('Y')
		plt.title('Streamlines')
		# Calculate the streamlines
		forward_stag_streamline = self.flow.streamlines(forward_stagnation_point[0]+forward_normal_stag[0] * 1e-3, forward_stagnation_point[1]+forward_normal_stag[1] * 1e-3, -1*delta_s)
		aft_stag_streamline = self.flow.streamlines(aft_stagnation_point[0]+aft_normal_stag[0] * 1e-3, aft_stagnation_point[1] + aft_normal_stag[1]*1e-3, delta_s)
		plt.plot(forward_stag_streamline[:, 0], forward_stag_streamline[:, 1],color='black')
		plt.plot(aft_stag_streamline[:, 0], aft_stag_streamline[:, 1],color='black')

		for i in range(n_lines):
			x = x_start
			y = -delta_y * (i+1) + forward_stag_streamline[-1,1]
			streamline = self.flow.streamlines(x, y, delta_s)
			plt.plot(streamline[:, 0], streamline[:, 1],color='black')
			y = delta_y * (i+1) + forward_stag_streamline[-1,1]
			streamline = self.flow.streamlines(x, y, delta_s)
			plt.plot(streamline[:, 0], streamline[:, 1],color='black')
		# plt.plot(forward_stagnation_point[0], forward_stagnation_point[1], 'ro', label='Forward Stagnation Point')
		# plt.plot(aft_stagnation_point[0], aft_stagnation_point[1], 'ro', label='Aft Stagnation Point')


		plt.xlim(x_lower_limit, x_upper_limit)
		plt.ylim(x_lower_limit, x_upper_limit)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()


	def run(self):
		"""
		Executes the main logic of the application.
		"""
		self.load_config()
		self.setup_Geometry()
		self.load_flow_field(self.free_stream_velocity, self.angle_of_attack)

		self.plot(self.x_start, self.x_low_val, self.x_up_val, self.delta_s, self.n_lines, self.delta_y)

		# streamlines_upper = self.flow.streamlines(-4.5,-.25,0.01)

		# plt.plot(streamlines_upper[:,0], streamlines_upper[:,1],label='Streamlines Upper')
		# plt.show()
		# forward_stagnation_point, aft_stagnation_point = self.stagnation_point()
		# print("Forward Stagnation Point: ", forward_stagnation_point)
		# print("Aft Stagnation Point: ", aft_stagnation_point)

		

		# x = np.linspace(-1 * self.radius, self.radius, 1000)
		# camber = np.empty((2, 0))
		# upper_surface = np.empty((2, 0))
		# lower_surface = np.empty((2, 0))

		# for i in range(len(x)):
		# 	camber_temp, upper_surface_temp, lower_surface_temp = self.geometry.circle(x[i])

		# 	# Append the new values to the arrays
		# 	camber = np.hstack((camber, camber_temp.reshape(2, 1)))
		# 	upper_surface = np.hstack((upper_surface, upper_surface_temp.reshape(2, 1)))
		# 	lower_surface = np.hstack((lower_surface, lower_surface_temp.reshape(2, 1)))


		# x_coord = 2
		# _,point,_ = self.geometry.circle(x_coord)
		# y_coord = point[1]
		# tangent_upper, tangent_lower = self.surface_tangent(x_coord)
		# normal_upper, normal_lower = self.surface_normal(x_coord)

		

		# print("tangent_upper: ", tangent_upper)
		# print("tangent_lower: ", tangent_lower)
		# print("normal_upper: ", normal_upper)
		# print("normal_lower: ", normal_lower)

		# # Plotting the results
		# plt.figure()
		# plt.plot(camber[0, :], camber[1, :], label='Camber')
		# plt.plot(upper_surface[0, :], upper_surface[1, :], label='Upper Surface')
		# plt.plot(lower_surface[0, :], lower_surface[1, :], label='Lower Surface')
		# # Plot tangent and normal vectors

		# plt.quiver(x_coord, y_coord, tangent_upper[0], tangent_upper[1], color='green', scale=10, label='Tangent Upper')
		# plt.quiver(x_coord, y_coord, normal_upper[0], normal_upper[1], color='purple', scale=10, label='Normal Upper')
		# plt.quiver(x_coord, -y_coord, tangent_lower[0], tangent_lower[1], color='red', scale=10, label='Tangent Lower')
		# plt.quiver(x_coord, -y_coord, normal_lower[0], normal_lower[1], color='yellow', scale=10, label='Normal Lower')
		# plt.legend(loc='upper right')
		# plt.xlabel('X')
		# plt.ylabel('Y')
		# plt.title('Airfoil Geometry')
		# plt.show()

if __name__ == "__main__":
	main = Main('input.json')
	main.run()
		