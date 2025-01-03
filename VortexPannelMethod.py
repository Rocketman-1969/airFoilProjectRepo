import numpy as np

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

        # Regularize the matrix to avoid numerical issues
        A += 1e-10 * np.eye(A.shape[0])
        print(f"Matrix A:\n{A}")
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
        gamma = np.linalg.solve(A,B)
        return gamma
    
    def run(self, x_all, y_all):
        x_cp_all, y_cp_all = zip(*[self.get_control_points(x, y) for x, y in zip(x_all, y_all)])
        A = self.get_A_matrix(x_all, y_all, x_cp_all, y_cp_all)
        B = self.get_B_matrix(x_all, y_all)
        gamma = self.get_gamma(A, B)
        return gamma
