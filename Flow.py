import numpy as np

class Flow:

    def __init__(self, radius, V_inf, alpha, x_low_val, x_up_val):
        self.radius = radius
        self.V_inf = V_inf
        self.alpha = alpha
        self.x_low_val = x_low_val
        self.x_up_val = x_up_val


    def flow_over_cylinder_cartesin(self, x, y):
        """
        Calculates the camber line, upper surface, and lower surface of the circle.

        Returns
        -------
        tuple
            A tuple containing the camber line, upper surface, and lower surface as numpy arrays.
        """
        r = np.sqrt(x**2 + y**2)
        costheta = x / r
        sintheta = y / r

        f_1 = self.V_inf * (1 - (self.radius**2 / r**2))*((costheta)*np.cos(self.alpha) + (sintheta)*np.sin(self.alpha))
        f_2 = -self.V_inf * (1 + (self.radius**2 / r**2))*((sintheta)*np.cos(self.alpha) - (costheta)*np.cos(self.alpha))

        velocity = np.array([f_1*(costheta) - f_2*(sintheta), f_1*(sintheta) + f_2*(costheta)])

        # f_1 = self.V_inf * (1 - (self.radius**2 / (x**2 + y**2)))*((x / np.sqrt(x**2 + y**2))*np.cos(self.alpha) + (y / np.sqrt(x**2 + y**2))*np.sin(self.alpha))
        # f_2 = self.V_inf * (1 + (self.radius**2 / (x**2 + y**2)))*((y / np.sqrt(x**2 + y**2))*np.cos(self.alpha) - (y / np.sqrt(x**2 + y**2))*np.cos(self.alpha))

        #velocity = np.array([f_1*(x/(np.sqrt(x**2 + y**2))) - f_2*(y/np.sqrt(x**2 + y**2)), f_1*(y/(np.sqrt(x**2 + y**2))) + f_2*(x/np.sqrt(x**2 + y**2))])

        return velocity
    
    def test_flow(self, x, y):

        velocity = np.array([x, 0])

        return velocity
    
    def streamlines(self, x,y, delta_s):
        """
        Calculate the streamlines at a given x-coordinate.
    
        Parameters:
        x (float): The x-coordinate at which to calculate the streamlines.
        delta_s (float): The step size for the streamlines.
    
        Returns:
        tuple: A tuple containing the streamlines for the upper and lower surfaces.
        """
        streamline = []

        while True:
            k1 = delta_s * self.flow_over_cylinder_cartesin(x, y)
            k2 = delta_s * self.flow_over_cylinder_cartesin(x + 0.5*k1[0], y + 0.5*k1[1])
            k3 = delta_s * self.flow_over_cylinder_cartesin(x + 0.5*k2[0], y + 0.5*k2[1])
            k4 = delta_s * self.flow_over_cylinder_cartesin(x + k3[0], y + k3[1])

            x_new = x + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
            y_new = y + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
            
            streamline.append([x_new, y_new])

            x = x_new
            y = y_new

            print(x_new, y_new)
            if x_new < self.x_low_val or x_new > self.x_up_val:
                break

            
        return streamline