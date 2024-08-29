import numpy as np
import matplotlib.pyplot as plt
import json
import GeometeryClass
import Flow


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
		self.x_low_val = json_vals['plot']['x_lower_limit']
		self.x_up_val = json_vals['plot']['x_upper_limit']

	def setup_Geometry(self):
		"""
		Initializes the Geometery object and calculates the camber, upper surface, and lower surface.
		"""
		self.geometry = GeometeryClass.Geometery(self.radius)
		
	def load_flow_field(self, V_inf, alpha):
		"""
		Loads the flow field parameters.

		Parameters:
		V_inf (float): The free stream velocity.
		alpha (float): The angle of attack.
		"""
		self.flow = Flow.Flow(self.radius, V_inf, alpha)

	def surface_normal(self, x):
		"""
		Calculate the surface normal vectors at a given x-coordinate.

		Parameters:
		x (float): The x-coordinate at which to calculate the surface normals.

		Returns:
		tuple: A tuple containing the unit normal vectors for the upper and lower surfaces.
		"""
		normal_upper = np.array([x, np.sqrt(self.radius**2 - x**2)])
		normal_lower = np.array([x, np.sqrt(self.radius**2 - x**2)])

		unit_normal_upper = normal_upper / np.linalg.norm(normal_upper)
		unit_normal_lower = normal_lower / np.linalg.norm(normal_lower)

		return unit_normal_upper, unit_normal_lower
	
	def surface_tangent(self, x):
		"""
		Calculate the surface tangent vectors at a given x-coordinate.

		Parameters:
		x (float): The x-coordinate at which to calculate the surface tangents.

		Returns:
		tuple: A tuple containing the unit tangent vectors for the upper and lower surfaces.
		"""
		unit_tangent_upper, unit_tangent_lower = self.surface_normal(x)

		z = np.array([0, 0, 1])

		tangent_upper = np.cross(unit_tangent_upper, z)
		tangent_lower = np.cross(-unit_tangent_lower, z)

		unit_tangent_upper = tangent_upper / np.linalg.norm(tangent_upper)
		unit_tangent_lower = tangent_lower / np.linalg.norm(tangent_lower)

		unit_tangent_lower = unit_tangent_lower[:2]
		unit_tangent_upper = unit_tangent_upper[:2]

		return unit_tangent_upper, unit_tangent_lower
	
	def surface_tangential_velocity(self, x):
		"""
		Calculate the tangential velocity at a given x-coordinate.

		Parameters:
		x (float): The x-coordinate at which to calculate the tangential velocity.

		Returns:
		tuple: A tuple containing the tangential velocity for the upper and lower surfaces.
		"""
		unit_tangent_upper, unit_tangent_lower = self.surface_tangent(x)

		velocity_upper = np.dot(self.flow.flow_over_cylinder_cartesin(x, np.sqrt(self.radius**2 - x**2)), unit_tangent_upper)
		velocity_lower = np.dot(self.flow.flow_over_cylinder_cartesin(x, -np.sqrt(self.radius**2 - x**2)), unit_tangent_lower)

		return velocity_upper, velocity_lower

	
	def run(self):
		"""
		Executes the main logic of the application.
		"""
		self.load_config()
		self.setup_Geometry()
		self.load_flow_field(1, 0)

		Velocity = self.flow.flow_over_cylinder_cartesin(4, 0)
		print(Velocity)

		velocity_upper, velocity_lower = self.surface_tangential_velocity(1)
		print("Velocity Upper: ", velocity_upper)
		print("Velocity Lower: ", velocity_lower)




		# # set up X array
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

		# unit_normal_upper, unit_normal_lower = self.surface_normal(-2)
		# unit_tangent_upper, unit_tangent_lower = self.surface_tangent(-2)

		# # Print the results
		# print("normal_upper: ", unit_normal_upper)
		# print("normal_lower: ", unit_normal_lower)
		# print("tangent_upper: ", unit_tangent_upper)
		# print("tangent_lower: ", unit_tangent_lower)

		# Plotting the results
		# plt.figure()
		# plt.plot(camber[0, :], camber[1, :], label='Camber')
		# plt.plot(upper_surface[0, :], upper_surface[1, :], label='Upper Surface')
		# plt.plot(lower_surface[0, :], lower_surface[1, :], label='Lower Surface')
		# plt.legend(loc='upper right')
		# plt.xlabel('X')
		# plt.ylabel('Y')
		# plt.title('Airfoil Geometry')
		# plt.show()

if __name__ == "__main__":
	main = Main('input.json')
	main.run()
		