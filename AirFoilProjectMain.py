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
		self.plot()

if __name__ == "__main__":
	main = Main('input.json')
	main.run()

