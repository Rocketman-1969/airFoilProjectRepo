import numpy as np


class Geometery:
    """
    A class to represent the geometry of a circle.

    Methods
    -------
    circle():
        Calculates the camber line, upper surface, and lower surface of the circle.
    """

    def __init__(self, radius):
        """
        Constructs all the necessary attributes for the Geometery object.

        Parameters
        ----------
        radius : float
            The radius of the circle.
        """
        self.radius = radius

    def circle(self):
        """
        Calculates the camber line, upper surface, and lower surface of the circle.

        Returns
        -------
        tuple
            A tuple containing the camber line, upper surface, and lower surface as numpy arrays.
        """
        x_values = np.linspace(-self.radius, self.radius, 10000)
        upper_surface = np.array([x_values, np.sqrt(self.radius**2 - x_values**2)])
        lower_surface = np.array([x_values, -np.sqrt(self.radius**2 - x_values**2)])
        camber = np.array([x_values, np.zeros_like(x_values)])

        return camber, upper_surface, lower_surface
