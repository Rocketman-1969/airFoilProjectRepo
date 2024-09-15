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

    def circle(self,x):
        """
        Calculates the camber line, upper surface, and lower surface of the circle.

        Returns
        -------
        tuple
            A tuple containing the camber line, upper surface, and lower surface as numpy arrays.
        """

        upper_surface = np.array([x, np.sqrt(self.radius**2 - x**2)])
        lower_surface = np.array([x, -np.sqrt(self.radius**2 - x**2)])
        camber = np.array([x, 0])

        return camber, upper_surface, lower_surface
