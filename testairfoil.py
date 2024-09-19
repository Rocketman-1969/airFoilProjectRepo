import numpy as np

def generate_naca4_airfoil(naca, total_points):
    # Extract NACA parameters
    m = int(naca[0]) / 100.0  # Maximum camber
    p = int(naca[1]) / 10.0   # Position of maximum camber
    t = int(naca[2:]) / 100.0 # Thickness

    # Define step size for odd or even number of points
    if total_points % 2 == 1:  # Odd case
        delta_theta = np.pi / (total_points // 2)
        indices = np.arange(1, total_points // 2 + 1)
        x_cos = 0.5 * (1 - np.cos(indices * delta_theta))
        x_cos = np.insert(x_cos, total_points // 2, 0.0)
    else:  # Even case
        delta_theta = np.pi / ((total_points // 2) - 0.5)
        indices = np.arange(1, total_points // 2 + 1)
        x_cos = 0.5 * (1 - np.cos((indices - 0.5) * delta_theta))

    # Thickness distribution
    yt = 5 * t * (0.2969 * np.sqrt(x_cos) - 0.1260 * x_cos - 0.3516 * x_cos**2 + 0.2843 * x_cos**3 - 0.1015 * x_cos**4)

    # Camber line
    yc = np.where(x_cos < p, m / (p**2) * (2 * p * x_cos - x_cos**2), m / ((1 - p)**2) * ((1 - 2 * p) + 2 * p * x_cos - x_cos**2))
    dyc_dx = np.where(x_cos < p, 2 * m / (p**2) * (p - x_cos), 2 * m / ((1 - p)**2) * (p - x_cos))

    # Angle of the camber line
    theta = np.arctan(dyc_dx)

    # Upper and lower surface coordinates
    xu = x_cos - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x_cos + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    # Combine coordinates, avoiding duplication at leading edge
    x_coords = np.concatenate([xu[::-1], xl[1:]])
    y_coords = np.concatenate([yu[::-1], yl[1:]])

    return x_coords, y_coords

def save_airfoil_to_file(naca, total_points, filename):
    x_coords, y_coords = generate_naca4_airfoil(naca, total_points)
    
    # Write to file
    with open(filename, 'w') as file:
        file.write(f"NACA {naca} Airfoil\n")
        for x, y in zip(x_coords, y_coords):
            file.write(f"{x:.6f} {y:.6f}\n")

# Input from user
naca_number = input("Enter the NACA 4-digit series (e.g., 2412): ")
total_points = int(input("Enter the total number of points for the airfoil: "))
filename = input("Enter the output filename (e.g., airfoil.dat): ")

# Generate and save airfoil
save_airfoil_to_file(naca_number, total_points, filename)
print(f"NACA {naca_number} airfoil coordinates have been saved to {filename}.")
