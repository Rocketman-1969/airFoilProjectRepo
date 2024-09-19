import numpy as np


class Geometery:
    """
    A class to represent the geometry of a circle.

    Methods
    -------
    circle():
        Calculates the camber line, upper surface, and lower surface of the circle.
    """

    def __init__(self, radius, NACA, n_points):
        """
        Constructs all the necessary attributes for the Geometery object.

        Parameters
        ----------
        radius : float
            The radius of the circle.
        """
        self.radius = radius
        self.NACA = str(int(NACA))
        self.n_points = n_points

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
    
    def Cose_cluster(self, n_points):
            # Define step size for odd or even number of points
        if n_points % 2 == 1:  # Odd number of points
            # Equation (4.2.16): Calculate delta_theta
            delta_theta = np.pi / (n_points // 2)
            indices = np.arange(1, n_points // 2 + 1)
            x_cos = 0.5 * (1 - np.cos(indices * delta_theta))
            x_cos = np.insert(x_cos, 0, 0.0)
        else:  # Even case
             # Equation (4.2.19): Calculate delta_theta
            delta_theta = np.pi / (n_points / 2 - 0.5)
            indices = np.arange(1, n_points // 2 + 1)
            x_cos = 0.5 * (1 - np.cos(indices * delta_theta - 0.5 * delta_theta))
            
        return x_cos


    def generate_naca4_airfoil(self, naca, x):
        # Convert naca to string if it's an integer
        naca = str(naca)
        x = np.atleast_1d(x)  # Ensure x is an array

        # Extract NACA parameters
        m = int(naca[0]) / 100.0  # Maximum camber
        p = int(naca[1]) / 10.0   # Position of maximum camber
        t = int(naca[2:]) / 100.0 # Thickness

        # Thickness distribution
        yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)

        # Camber line
        if m == 0:
            yc = np.zeros_like(x)
            dyc_dx = np.zeros_like(x)
        else:
            yc = np.where(x < p, m / (p**2) * (2 * p * x - x**2), m / ((1 - p)**2) * ((1 - 2 * p) + 2 * p * x - x**2))
            dyc_dx = np.where(x < p, 2 * m / (p**2) * (p - x), 2 * m / ((1 - p)**2) * (p - x))

        # Angle of the camber line
        theta = np.arctan(dyc_dx)

        # Upper and lower surface coordinates for each x value
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)

        # If a single x value was provided, return just the corresponding upper and lower surface points
        if x.size == 1:
            return np.array([xu[0], yu[0]]), np.array([xl[0], yl[0]])

        # Combine coordinates for the full airfoil, including both leading and trailing edges
        elif self.n_points % 2 == 1:
            x_coords = np.concatenate([xu[::-1], xl[1:]])
            y_coords = np.concatenate([yu[::-1], yl[1:]])
            
        else:
            x_coords = np.concatenate([xu[::-1], xl])
            y_coords = np.concatenate([yu[::-1], yl])

        return x_coords, y_coords
