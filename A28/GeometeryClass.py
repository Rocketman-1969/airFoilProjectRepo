import numpy as np

class Geometery:
    """
    A class to represent the geometry of a circle.

    Methods
    -------
    circle():
        Calculates the camber line, upper surface, and lower surface of the circle.
    """

    def __init__(self, radius, z0):
        """
        Constructs all the necessary attributes for the Geometery object.

        Parameters
        ----------
        radius : float
            The radius of the circle.
        """
        self.radius = radius
        self.z0 = z0

    def circle(self,x):
        """
        Calculates the camber line, upper surface, and lower surface of the circle.

        Returns
        -------
        tuple
            A tuple containing the camber line, upper surface, and lower surface as numpy arrays.
        """

        # if abs(x - self.z0[0]) > self.radius:
        #     raise ValueError("x must be within the range of the circle's radius.")

        # Calculating upper and lower surface points relative to z0 (the circle's origin)
        upper_surface = np.array([x, self.z0[1] + np.sqrt(self.radius**2 - (x - self.z0[0])**2)])
        lower_surface = np.array([x, self.z0[1] - np.sqrt(self.radius**2 - (x - self.z0[0])**2)])
        camber = np.array([x, self.z0[1]])  # Camber line is constant at z0[1]

        return camber, upper_surface, lower_surface
