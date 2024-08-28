import numpy as np
import pandas as pd
import scipy as sp
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

	def setup_geometry(self):
		"""
		Initializes the Geometery object and calculates the camber, upper surface, and lower surface.
		"""
		self.geometry = GeometeryClass.Geometery(self.radius)
		self.camber, self.upper_surface, self.lower_surface = self.geometry.circle()

	def surface_tangent(self, x):
		
		y_upper = np.interp(x, self.upper_surface[0], self.upper_surface[1])
		y_lower = np.interp(x, self.lower_surface[0], self.lower_surface[1])

		idx_upper = np.searchsorted(self.upper_surface[0], x)
		idx_lower = np.searchsorted(self.lower_surface[0], x)
		
		self.upper_surface = np.insert(self.upper_surface, idx_upper, [x, y_upper], axis=1)
		self.lower_surface = np.insert(self.lower_surface, idx_lower, [x, y_lower], axis=1)
		
		dx_upper = self.upper_surface[0][idx_upper + 1] - self.upper_surface[0][idx_upper-1]
		dy_upper = self.upper_surface[1][idx_upper + 1] - self.upper_surface[1][idx_upper-1]

		dx_lower = self.lower_surface[0][idx_lower + 1] - self.lower_surface[0][idx_lower-1]
		dy_lower = self.lower_surface[1][idx_lower + 1] - self.lower_surface[1][idx_lower-1]

		tangent_upper = np.array([dx_upper, dy_upper])
		tangent_lower = np.array([dx_lower, dy_lower])
		
		unit_tangent_upper = tangent_upper / np.linalg.norm(tangent_upper)
		unit_tangent_lower = tangent_lower / np.linalg.norm(tangent_lower)

		return unit_tangent_upper, unit_tangent_lower
	
	def surface_normal(self, x):
		
		y_upper = np.interp(x, self.upper_surface[0], self.upper_surface[1])
		y_lower = np.interp(x, self.lower_surface[0], self.lower_surface[1])

		idx_upper = np.searchsorted(self.upper_surface[0], x)
		idx_lower = np.searchsorted(self.lower_surface[0], x)
		
		self.upper_surface = np.insert(self.upper_surface, idx_upper, [x, y_upper], axis=1)
		self.lower_surface = np.insert(self.lower_surface, idx_lower, [x, y_lower], axis=1)
		
		dx_upper = self.upper_surface[0][idx_upper + 1] - self.upper_surface[0][idx_upper-1]
		dy_upper = self.upper_surface[1][idx_upper + 1] - self.upper_surface[1][idx_upper-1]

		dx_lower = self.lower_surface[0][idx_lower + 1] - self.lower_surface[0][idx_lower-1]
		dy_lower = self.lower_surface[1][idx_lower + 1] - self.lower_surface[1][idx_lower-1]

		normal_upper = np.array([dy_upper, dx_upper])
		normal_lower = np.array([-dy_lower, -dx_lower])
		
		unit_tangent_upper = normal_upper / np.linalg.norm(normal_upper)
		unit_tangent_lower = normal_lower / np.linalg.norm(normal_lower)

		return unit_tangent_upper, unit_tangent_lower

		
	def plot(self):
		"""
		Plots the camber line, upper surface, and lower surface.
		"""
		plt.figure()
		plt.plot(self.camber[0], self.camber[1], label='Camber Line')
		plt.plot(self.upper_surface[0], self.upper_surface[1], color='blue', label='Upper Surface')
		plt.plot(self.lower_surface[0], self.lower_surface[1], color='red', label='Lower Surface')
		plt.xlim(self.x_low_val, self.x_up_val)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.legend(loc='upper right')
		plt.show()

	def run(self):
		"""
		Executes the main logic of the application.
		"""
		self.load_config()
		self.setup_geometry()
		#self.plot()
		unit_tangent_upper, unit_tangent_lower = self.surface_tangent(-1.99999)
		unit_normal_upper, unit_normal_lower = self.surface_normal(0)
		print(f'Unit tangent vector on the upper surface at x=0: {unit_tangent_upper}')
		print(f'Unit tangent vector on the lower surface at x=0: {unit_tangent_lower}')
		print(f'Unit normal vector on the upper surface at x=0: {unit_normal_upper}')
		print(f'Unit normal vector on the lower surface at x=0: {unit_normal_lower}')
		
		

if __name__ == "__main__":
	main = Main('input.json')
	main.run()

