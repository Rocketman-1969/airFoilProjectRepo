import numpy as np

class Flow:

    def __init__(self, radius, V_inf, alpha):
        self.radius = radius
        self.V_inf = V_inf
        self.alpha = alpha

    def flow_over_cylinder_cartesin(self, x, y):
        """
        Calculates the camber line, upper surface, and lower surface of the circle.

        Returns
        -------
        tuple
            A tuple containing the camber line, upper surface, and lower surface as numpy arrays.
        """

        f_1 = self.V_inf * (1 - (self.radius**2 / (x**2 + y**2)))*((x / np.sqrt(x**2 + y**2))*np.cos(self.alpha) + (y / np.sqrt(x**2 + y**2))*np.sin(self.alpha))
        f_2 = self.V_inf * (1 + (self.radius**2 / (x**2 + y**2)))*((y / np.sqrt(x**2 + y**2))*np.cos(self.alpha) - (y / np.sqrt(x**2 + y**2))*np.cos(self.alpha))

        velocity = np.array([f_1*(x/(np.sqrt(x**2 + y**2))) - f_2*(y/np.sqrt(x**2 + y**2)), f_1*(y/(np.sqrt(x**2 + y**2))) + f_2*(x/np.sqrt(x**2 + y**2))])

        return velocity