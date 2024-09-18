import numpy as np

class Flow:

    def __init__(self, elements, x_low_val, x_up_val):
        #print(self.elements)
        self.elements = elements
        self.x_low_val = x_low_val
        self.x_up_val = x_up_val
    
    def freestream(self, velocity, alpha, point):
        """
        Set the freestream conditions.
    
        Parameters:
        velocity (float): The freestream velocity.
        alpha (float): The angle of attack in degrees.
        """
        r = np.sqrt((point[0])**2 + (point[1])**2)
        theta = np.arctan2(point[1], point[0])
        alpha = np.deg2rad(alpha)

        r_dot = velocity * np.cos(theta - alpha)
        theta_dot = -velocity * np.sin(theta - alpha)

        x_dot = r_dot * np.cos(theta) - theta_dot * np.sin(theta)
        y_dot = r_dot * np.sin(theta) + theta_dot * np.cos(theta)
        

        return x_dot, y_dot
    
    def source(self, strength, xs, ys, point):
        """
        Set the source strength and location.
    
        Parameters:
        strength (float): The strength of the source.
        xs (float): The x-coordinate of the source.
        ys (float): The y-coordinate of the source.
        """
        r = np.sqrt((point[0] - xs)**2 + (point[1] - ys)**2)
        theta = np.arctan2(point[1] - ys, point[0] - xs)

        r_dot = strength / (2 * np.pi * r)
        theta_dot = 0

        x_dot = r_dot * np.cos(theta) - theta_dot * np.sin(theta)
        y_dot = r_dot * np.sin(theta) + theta_dot * np.cos(theta)     

        return x_dot, y_dot
    
    def vortex(self, strength, xs, ys, point):
        """
        Set the vortex strength and location.
    
        Parameters:
        strength (float): The strength of the vortex.
        xs (float): The x-coordinate of the vortex.
        ys (float): The y-coordinate of the vortex.
        """
        r = np.sqrt((point[0] - xs)**2 + (point[1] - ys)**2)
        theta = np.arctan2(point[1] - ys, point[0] - xs)

        r_dot = 0
        theta_dot = -strength / (2 * np.pi * r)

        x_dot = r_dot * np.cos(theta) - theta_dot * np.sin(theta)
        y_dot = r_dot * np.sin(theta) + theta_dot * np.cos(theta)
        

        return x_dot, y_dot
    
    def doublet(self, strength, xs, ys, alpha, point):
        """
        Set the doublet strength and location.
    
        Parameters:
        strength (float): The strength of the doublet.
        xs (float): The x-coordinate of the doublet.
        ys (float): The y-coordinate of the doublet.
        alpha (float): The angle of the doublet.
        """
        alpha = np.deg2rad(alpha)

        r = np.sqrt((point[0] - xs)**2 + (point[1] - ys)**2)
        theta = np.arctan2(point[1] - ys, point[0] - xs)

        r_dot = -(strength*np.cos(theta-alpha)) / (2 * np.pi * r**2)
        theta_dot = -(strength*np.sin(theta-alpha)) / (2 * np.pi * r**2)

        x_dot = r_dot * np.cos(theta) - theta_dot * np.sin(theta)
        y_dot = r_dot * np.sin(theta) + theta_dot * np.cos(theta)

        

        return x_dot, y_dot
    
    def velocity_field(self, X, Y):
        """
        Calculate the velocity field.
    
        Parameters:
        X (2D array): The x-coordinates of the mesh grid.
        Y (2D array): The y-coordinates of the mesh grid.
    
        Returns:
        tuple: A tuple containing the x-component and y-component of the velocity field.
        """
        Vx = 0
        Vy = 0
        
        for element in self.elements:
            if element['type'] == 'freestream':
                u, v = self.freestream(element['velocity'], element['angle_of_attack'], [X, Y])
            elif element['type'] == 'source':
                u, v = self.source(element['lambda'], element['x'], element['y'], [X, Y])
            elif element['type'] == 'vortex':
                u, v = self.vortex(element['gamma'], element['x'], element['y'], [X, Y])
            elif element['type'] == 'doublet':
                u, v = self.doublet(element['kappa'], element['x'], element['y'], element['alpha'], [X, Y])
            Vx += u
            Vy += v
            
        return np.array([Vx, Vy])        
    
    def unit_velocity(self, x, y):
        velocity = self.velocity_field(x, y)
        
        return velocity

    
    def streamlines(self, x, y, delta_s):
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
            k1 = delta_s * self.unit_velocity(x, y)
            k2 = delta_s * self.unit_velocity(x + 0.5*k1[0], y + 0.5*k1[1])
            k3 = delta_s * self.unit_velocity(x + 0.5*k2[0], y + 0.5*k2[1])
            k4 = delta_s * self.unit_velocity(x + k3[0], y + k3[1])

            x_new = x + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
            y_new = y + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
            
            streamline.append([x_new, y_new])

            x = x_new
            y = y_new

            if x_new < self.x_low_val or x_new > self.x_up_val:
                break

            
        return np.array(streamline)