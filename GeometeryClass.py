import numpy as np


class Geometery:
    """
    A class to represent the geometry of a circle.

    Attributes
    ----------
    radius : float
        The radius of the circle.

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
        angles = np.linspace(0, np.pi, 100)
        upper_surface = np.array([self.radius * np.cos(angles), self.radius * np.sin(angles)])
        lower_surface = np.array([self.radius * np.cos(angles), -self.radius * np.sin(angles)])
        camber = np.array([[self.radius, -self.radius], [0, 0]])

        return camber, upper_surface, lower_surface
