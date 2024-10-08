import numpy as np
import matplotlib.pyplot as plt
import json
import GeometeryClass
from Flow import Flow
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
		self.elements = []
		self.load_config()
		

	def load_config(self):
		"""
		Loads the configuration from the JSON file.
		"""
		with open(self.config_file, 'r') as file:
			config = json.load(file)
			self.elements = list(config['elements'].values())
			self.plot_config = config['plot']
			
	def load_flow_field(self):
		"""
		Loads the flow field parameters.
		"""
		self.flow = Flow(self.elements, self.plot_config["x_lower_limit"], self.plot_config["x_upper_limit"])
		
		# Initialize an empty dictionary to store elements
		self.element_dict = {}
		
		# Loop through the elements and save each element in the dictionary
		for index, element in enumerate(self.elements):
			self.element_dict[f'element{index}'] = element


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
		if np.abs(x-self.radius) < delta:
			_, y_upper_minus, y_lower_minus = self.geometry.circle(x-delta)
			tangent_upper = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), -1*(y_upper_minus[1]-y_lower_minus[1])])
			tangent_lower = np.array([np.abs(y_upper_minus[0]-y_upper_minus[0]), -1*(y_upper_minus[1]-y_lower_minus[1])])
		elif np.abs(x+self.radius) < delta:
			_,y_upper_plus, y_lower_plus = self.geometry.circle(x+delta)

			tangent_upper = np.array([-1*np.abs(y_upper_plus[0]-y_upper_plus[0]), y_upper_plus[1]-y_lower_plus[1]])
			tangent_lower = np.array([-1*np.abs(y_upper_plus[0]-y_upper_plus[0]), y_upper_plus[1]-y_lower_plus[1]])
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
			velocity_upper = np.dot(self.flow.flow_over_cylinder_circulation(point[0], point[1]), unit_tangent_upper)
			return velocity_upper
		else:
			_, _, point = self.geometry.circle(x)
			velocity_lower = np.dot(self.flow.flow_over_cylinder_circulation(point[0], point[1]), unit_tangent_lower)
			return velocity_lower

	def velocity_derivative(self, x, upper=True):
		"""
		Calculate the derivative of the velocity at a given x-coordinate.

		Parameters:
		x (float): The x-coordinate at which to calculate the velocity derivative.

		Returns:
		tuple: A tuple containing the velocity derivative for the upper and lower surfaces.
		"""
		delta = 1e-5
		if upper:
			velocity_upper_plus = self.surface_tangential_velocity(x + delta, upper=True)
			velocity_upper_minus = self.surface_tangential_velocity(x - delta, upper=True)
			velocity_derivative_upper = (velocity_upper_plus - velocity_upper_minus) / (2 * delta)
			return velocity_derivative_upper
		else:
			velocity_lower_plus = self.surface_tangential_velocity(x + delta, upper=False)
			velocity_lower_minus = self.surface_tangential_velocity(x - delta, upper=False)
			velocity_derivative_lower = (velocity_lower_plus - velocity_lower_minus) / (2 * delta)
			return velocity_derivative_lower

	
	def stagnation_point(self):
		"""
        Calculate the stagnation point of the flow field.

        Returns:
        tuple: A tuple containing the x and y coordinates of the forward and aft stagnation points.
        """
		epsilon = 1e-5  # Small perturbation for numerical differentiation

        # Check leading edge (LE) and trailing edge (TE)
		velocity_LE = self.surface_tangential_velocity(self.LE)
		velocity_TE = self.surface_tangential_velocity(self.TE)

		if np.abs(velocity_LE) < epsilon:
			forward_stagnation_point = [self.LE, 0]
		else:
			if velocity_LE < 0:
                # Use upper surface
				upper_flag = True
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				def velocity_prime(x):
					return self.velocity_derivative(x, upper=upper_flag)
				x_stag = newton(velocity_function, self.LE + 0.01, fprime=velocity_prime, tol=epsilon)
				_, forward_stagnation_point, _ = self.geometry.circle(x_stag)
			else:
                # Use lower surface
				upper_flag = False
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				def velocity_prime(x):
					return self.velocity_derivative(x, upper=upper_flag)
				x_stag = newton(velocity_function, self.LE + 0.01, fprime=velocity_prime, tol=epsilon)
				_, _, forward_stagnation_point = self.geometry.circle(x_stag)

		if np.abs(velocity_TE) < epsilon:
			aft_stagnation_point = [self.TE, 0]
		else:
			if velocity_LE < 0:
                # Use lower surface
				upper_flag = True
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				def velocity_prime(x):
					return self.velocity_derivative(x, upper=upper_flag)
				x_stag = newton(velocity_function, self.TE - 0.01, fprime=velocity_prime, tol=epsilon)
				_, _, aft_stagnation_point = self.geometry.circle(x_stag)
			else:
                # Use upper surface
				upper_flag = False
				def velocity_function(x):
					return self.surface_tangential_velocity(x, upper=upper_flag)
				def velocity_prime(x):
					return self.velocity_derivative(x, upper=upper_flag)
				x_stag = newton(velocity_function, self.TE - 0.01, fprime=velocity_prime, tol=epsilon)
				_, aft_stagnation_point, _ = self.geometry.circle(x_stag)

		return forward_stagnation_point, aft_stagnation_point

	
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
		# Set up the plot
		plt.figure()
		plt.xlabel('X')
		plt.ylabel('Y')
		plt.title('Streamlines')

		x = x_start
		y = 0
		streamline = self.flow.streamlines(x, y, delta_s)
		plt.plot(streamline[:, 0], streamline[:, 1],color='black')

		for i in range(n_lines):
			x = x_start
			y = -delta_y * (i+1)
			streamline = self.flow.streamlines(x, y, delta_s)
			plt.plot(streamline[:, 0], streamline[:, 1],color='black')
			y = delta_y * (i+1)
			streamline = self.flow.streamlines(x, y, delta_s)
			plt.plot(streamline[:, 0], streamline[:, 1],color='black')

		for element in self.elements:
			if element['type'] == 'freestream':
				continue  # Ignore freestream elements

			x = element.get('x', 0)
			y = element.get('y', 0)
			if element['type'] == 'source':
				magnitude = element['lambda']
				marker = '*'
			elif element['type'] == 'doublet':
				magnitude = element['kappa']
				marker = 'D'
			elif element['type'] == 'vortex':
				magnitude = element['gamma']
				marker = 'x'

			color = 'blue' if magnitude > 0 else 'red'
			plt.scatter(x, y, color=color, s=10, marker=marker)
		
		plt.xlim(x_lower_limit, x_upper_limit)
		plt.ylim(x_lower_limit, x_upper_limit)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()



	def run(self):
		"""
		Executes the main logic of the application.
		"""
		self.load_config()
		
		self.load_flow_field()

		self.plot(self.plot_config['x_start'], self.plot_config['x_lower_limit'], self.plot_config['x_upper_limit'], self.plot_config['delta_s'], self.plot_config['n_lines'], self.plot_config['delta_y'])
		
		
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
	main = Main('A11/a11_input.json')
	main.run()
		