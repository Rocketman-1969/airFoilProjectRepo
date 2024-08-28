import numpy as np
import matplotlib.pyplot as plt
import json
import GeometeryClass


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
		

	def surface_tangent(self, x):
		
		tangent_upper = np.array([x, np.sqrt(self.radius**2 - x**2)])
		tangent_lower = np.array([x, -np.sqrt(self.radius**2 - x**2)])

		unit_tangent_upper = tangent_upper / np.linalg.norm(tangent_upper)
		unit_tangent_lower = tangent_lower / np.linalg.norm(tangent_lower)

		return unit_tangent_upper, unit_tangent_lower
	
	def surface_normal(self, x):
		
		unit_tangent_upper, unit_tangent_lower = self.surface_tangent(x)

		z = np.array([0, 0, 1])

		tangent_upper = np.cross(unit_tangent_upper, z)
		tangent_lower = np.cross(unit_tangent_lower, z)

		unit_tangent_upper = tangent_upper / np.linalg.norm(tangent_upper)
		unit_tangent_lower = tangent_lower / np.linalg.norm(tangent_lower)

		return unit_tangent_upper, unit_tangent_lower

		
	def plot(self, camber, upper_surface, lower_surface):
		plt.figure()
		plt.plot(camber[0], camber[1], label='Camber Line')
		plt.plot(upper_surface[0], upper_surface[1], color='blue', label='Upper Surface')
		plt.plot(lower_surface[0], lower_surface[1], color='red', label='Lower Surface')
		plt.xlim(self.x_low_val, self.x_up_val)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.legend(loc='upper right')
		plt.xlabel('X')
		plt.ylabel('Y')
		plt.title('Airfoil Geometry')
		plt.show()

	def run(self):
		"""
		Executes the main logic of the application.
		"""
		self.load_config()
		self.setup_Geometry()

		# set up X array
		x = np.linspace(-1 * self.radius, self.radius, 1000)
		camber = np.empty((2, 0))
		upper_surface = np.empty((2, 0))
		lower_surface = np.empty((2, 0))

		for i in range(len(x)):
			camber_temp, upper_surface_temp, lower_surface_temp = self.geometry.circle(x[i])

			# Append the new values to the arrays
			camber = np.hstack((camber, camber_temp.reshape(2, 1)))
			upper_surface = np.hstack((upper_surface, upper_surface_temp.reshape(2, 1)))
			lower_surface = np.hstack((lower_surface, lower_surface_temp.reshape(2, 1)))

		# Plotting the results
		plt.figure()
		plt.plot(camber[0, :], camber[1, :], label='Camber')
		plt.plot(upper_surface[0, :], upper_surface[1, :], label='Upper Surface')
		plt.plot(lower_surface[0, :], lower_surface[1, :], label='Lower Surface')
		plt.legend(loc='upper right')
		plt.xlabel('X')
		plt.ylabel('Y')
		plt.title('Airfoil Geometry')
		plt.show()

		
		

if __name__ == "__main__":
	main = Main('input.json')
	main.run()
		



		
		
		

if __name__ == "__main__":
	main = Main('input.json')
	main.run()

