import numpy as np

class Geometery:
    """
    A class to represent the geometry of a circle.

    Methods
    -------
    circle():
        Calculates the camber line, upper surface, and lower surface of the circle.
    """

    def __init__(self, radius, z0, epsillon):
        """
        Constructs all the necessary attributes for the Geometery object.

        Parameters
        ----------
        radius : float
            The radius of the circle.
        """
        self.radius = radius
        self.z0 = z0
        self.epsillon = epsillon
    
    def zeta_2_z(self, zeta):

        return zeta + ((self.radius-self.epsillon)**2)/zeta

    def geometery_zeta(self,x):
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

    def geometery_zplane(self, x):

        camber_zeta, upper_surface_zeta, lower_surface_zeta = self.geometery_zeta(x)
        camber_zeta = camber_zeta[0]+1j*camber_zeta[1]
        upper_surface_zeta = upper_surface_zeta[0]+1j*upper_surface_zeta[1]
        lower_surface_zeta = lower_surface_zeta[0]+1j*lower_surface_zeta[1]

        camber_z = self.zeta_2_z(camber_zeta)
        upper_surface_z = self.zeta_2_z(upper_surface_zeta)
        lower_surface_z = self.zeta_2_z(lower_surface_zeta)

        
        upper_surface_z = np.array([upper_surface_z.real, upper_surface_z.imag])
        lower_surface_z = np.array([lower_surface_z.real, lower_surface_z.imag])
        camber_z = np.array([(upper_surface_z[0]+lower_surface_z[0])/2, (upper_surface_z[1]+lower_surface_z[1])/2])
        
        return camber_z, upper_surface_z, lower_surface_z

    