import numpy as np
import matplotlib.pyplot as plt


class VortexPannelMethod:

    def __init__(self, chord, velocity, alpha):
        self.velocity = velocity
        self.alpha = alpha
        self.chord = chord
    
    def get_control_points(self, x, y):
        x_cp = (x[:-1] + x[1:]) / 2
        y_cp = (y[:-1] + y[1:]) / 2
        return x_cp, y_cp

    def get_length_of_jth_pannel(self, x, y, j):
        return np.sqrt((x[j+1] - x[j])**2 + (y[j+1] - y[j])**2)
    
    def get_xi_eta(self, x, y, x_cp, y_cp, l_j, j):
        matrix1 = np.array([[x[j+1] - x[j], y[j+1] - y[j]], 
                            [-(y[j+1] - y[j]), x[j+1] - x[j]]])
        matrix2 = np.array([x_cp - x[j], y_cp - y[j]])
        xi_eta = (1 / l_j) * np.matmul(matrix1, matrix2)
        return xi_eta[0], xi_eta[1]
    
    def get_phi(self, eta, xi, l_j):
        return np.arctan2(eta * l_j, eta**2 + xi**2 - xi * l_j)
    
    def get_psi(self, eta, xi, l_j):
        return 0.5 * np.log((xi**2 + eta**2) / (((xi - l_j)**2) + eta**2))
    
    def get_P_matrix(self, x, y, x_cp, y_cp, i, j):
        l_j = self.get_length_of_jth_pannel(x, y, j)
        xi, eta = self.get_xi_eta(x, y, x_cp, y_cp, l_j, j)
        phi = self.get_phi(eta, xi, l_j)
        psi = self.get_psi(eta, xi, l_j)
        
        matrix1 = np.array([[x[j+1] - x[j], -1 * (y[j+1] - y[j])],
                            [y[j+1] - y[j], x[j+1] - x[j]]])
        matrix2 = np.array([[(l_j - xi) * phi + eta * psi, xi * phi - eta * psi],
                            [eta * phi - (l_j - xi) * psi - l_j, 
                             -eta * phi - xi * psi + l_j]])
        P = (1 / (2 * np.pi * l_j**2)) * np.matmul(matrix1, matrix2)
        return P
    
    def get_A_matrix(self, x_all, y_all, x_cp_all, y_cp_all):
        n_total = sum(len(x) - 1 for x in x_all)  # Total number of panels
        A = np.zeros((n_total, n_total))
        offset_i = 0  # Offset for the rows
        for k_i, (x_i, y_i, x_cp_i, y_cp_i) in enumerate(zip(x_all, y_all, x_cp_all, y_cp_all)):
            offset_j = 0  # Offset for the columns
            for k_j, (x_j, y_j) in enumerate(zip(x_all, y_all)):
                for i in range(len(x_i) - 1):  # Loop through panels on airfoil i
                    for j in range(len(x_j) - 1):  # Loop through panels on airfoil j
                        P = self.get_P_matrix(x_j, y_j, x_cp_i[i], y_cp_i[i], i, j)
                        l_i = self.get_length_of_jth_pannel(x_i, y_i, i)

                        if offset_i + i < n_total and offset_j + j < n_total:
                            A[offset_i + i, offset_j + j] += ((x_i[i+1] - x_i[i]) / l_i) * P[1, 0] \
                                                            - ((y_i[i+1] - y_i[i]) / l_i) * P[0, 0]
                        if offset_i + i < n_total and offset_j + j + 1 < n_total:
                            A[offset_i + i, offset_j + j + 1] += ((x_i[i+1] - x_i[i]) / l_i) * P[1, 1] \
                                                           - ((y_i[i+1] - y_i[i]) / l_i) * P[0, 1]
                offset_j += len(x_j) - 1
        offset_i += len(x_i) - 1

        # Apply Kutta condition for each airfoil
        offset = 0
        for x in x_all:
            A[offset + len(x) - 2, offset] = 1.0
            A[offset + len(x) - 2, offset + len(x) - 2] = 1.0
            offset += len(x) - 1

        A += 1e-10 * np.eye(A.shape[0])  # Regularization
        return A

    def get_B_matrix(self, x_all, y_all):
        n_total = sum(len(x) - 1 for x in x_all)
        B = np.zeros((n_total, 1))
        alpha = np.radians(self.alpha)
        offset = 0
        for x, y in zip(x_all, y_all):
            for i in range(len(x) - 1):
                l_j = self.get_length_of_jth_pannel(x, y, i)
                B[offset + i, 0] = self.velocity * ((y[i+1] - y[i]) * np.cos(alpha) - (x[i+1] - x[i]) * np.sin(alpha)) / l_j
            offset += len(x) - 1
        return B

    def get_gamma(self, A, B):
        return np.linalg.solve(A, B)

    def run(self, x_all, y_all):
        x_cp_all, y_cp_all = zip(*[self.get_control_points(x, y) for x, y in zip(x_all, y_all)])
        A = self.get_A_matrix(x_all, y_all, x_cp_all, y_cp_all)
        B = self.get_B_matrix(x_all, y_all)
        gamma_global = self.get_gamma(A, B)

        # Split gamma into separate arrays for each airfoil
        gamma_split = []
        offset = 0
        for x in x_all:
            n_panels = len(x) - 1
            gamma_split.append(gamma_global[offset:offset + n_panels].flatten())
            offset += n_panels

        return gamma_split


    def induced_velocity(self, x, y, x_panel, y_panel, gamma):
        u, v = 0, 0
        for i in range(len(x_panel) - 1):
            x1, y1 = x_panel[i], y_panel[i]
            x2, y2 = x_panel[i + 1], y_panel[i + 1]
            dx, dy = x2 - x1, y2 - y1
            r1 = np.sqrt((x - x1)**2 + (y - y1)**2)
            r2 = np.sqrt((x - x2)**2 + (y - y2)**2)
            theta1 = np.arctan2(y - y1, x - x1)
            theta2 = np.arctan2(y - y2, x - x2)
            cross = dx * (y - y1) - dy * (x - x1)
            u += gamma[i] * cross * (1 / r1 - 1 / r2) / (2 * np.pi)
            v += gamma[i] * cross * (theta2 - theta1) / (2 * np.pi)
        return u, v

    def plot_streamlines(self, x_all, y_all, gamma_split, grid_size=100):
        """
        Plot streamlines around the airfoil(s).
        
        Parameters:
            x_all, y_all: List of x and y coordinates of all airfoils.
            gamma_split: List of circulation strengths for each airfoil.
            grid_size: Number of points in each direction for the velocity grid.
        """
        x_min, x_max = min(min(x) for x in x_all), max(max(x) for x in x_all)
        y_min, y_max = min(min(y) for y in y_all), max(max(y) for y in y_all)
        x_range = x_max - x_min
        y_range = y_max - y_min
        x_min -= 0.5 * x_range
        x_max += 0.5 * x_range
        y_min -= 0.5 * y_range
        y_max += 0.5 * y_range
        x_grid = np.linspace(x_min, x_max, grid_size)
        y_grid = np.linspace(y_min, y_max, grid_size)
        X, Y = np.meshgrid(x_grid, y_grid)
        U = np.full_like(X, self.velocity * np.cos(np.radians(self.alpha)))
        V = np.full_like(Y, self.velocity * np.sin(np.radians(self.alpha)))
        for x, y, g in zip(x_all, y_all, gamma_split):
            for i in range(grid_size):
                for j in range(grid_size):
                    u_induced, v_induced = self.induced_velocity(X[i, j], Y[i, j], x, y, g)
                    U[i, j] += u_induced
                    V[i, j] += v_induced
        plt.figure(figsize=(10, 6))
        for x, y in zip(x_all, y_all):
            plt.plot(x, y, color='black', linewidth=2)
        plt.streamplot(X, Y, U, V, color=np.sqrt(U**2 + V**2), linewidth=1.5, cmap='viridis')
        plt.colorbar(label='Velocity Magnitude')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Streamlines Around the Airfoils')
        plt.axis('equal')
        plt.grid()
        plt.show()



# Example Usage
if __name__ == "__main__":
    # Example airfoil geometries
    x1 = np.linspace(0, 1, 50)
    y1 = 0.1 * (0.2969 * np.sqrt(x1) - 0.126 * x1 - 0.3516 * x1**2 + 0.2843 * x1**3 - 0.1015 * x1**4)
    x2 = np.linspace(0, 1, 50)
    y2 = 0.08 * (0.2969 * np.sqrt(x2) - 0.126 * x2 - 0.3516 * x2**2 + 0.2843 * x2**3 - 0.1015 * x2**4) - 0.2

    # Combine airfoils
    x_all = [x1, x2]
    y_all = [y1, y2]

    # Initialize and run vortex panel method
    vpm = VortexPannelMethod(chord=1.0, velocity=1.0, alpha=5.0)
    gamma_split = vpm.run(x_all, y_all)

    # Plot streamlines
    vpm.plot_streamlines(x_all, y_all, gamma_split)

