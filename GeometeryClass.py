import numpy as np


class Geometery:
    def __init__(self, radius):
        self.radius = radius

    def circle(self):
        # Define the geometery of a circle
        angles = np.linspace(0, np.pi, 100)
        upper_surface = np.array([self.radius*np.cos(angles), self.radius*np.sin(angles)])
        lower_surface = np.array([self.radius*np.cos(angles), -self.radius*np.sin(angles)])
        camber = np.array([[self.radius,-self.radius], [0,0]])

        return camber, upper_surface, lower_surface
