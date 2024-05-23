import matplotlib
from tkinter import Tk
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numba import jit

# Constants
MU_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
I = 1  # Current (A)
a = 0.01  # Radius of the conductors (m)
d = 0.04  # Distance between the centers of the conductors (m)
Lx = 0.06  # Width of the grid (m)
Ly = 0.06  # Height of the grid (m)

@jit(nopython=True)
def r_tau(x, y, d):
    """
    Calculate the distances from a point to the centers of two conductors.
    """
    r_tau1 = np.sqrt((x - d / 2) ** 2 + y ** 2)
    r_tau2 = np.sqrt((x + d / 2) ** 2 + y ** 2)
    return r_tau1, r_tau2

@jit(nopython=True)
def vector_potential(x, y, a, d, I):
    """
    Calculate the magnetic vector potential at a point due to two cylindrical conductors.
    """
    r_tau1, r_tau2 = r_tau(x, y, d)
    if r_tau1 <= a and r_tau2 >= a:
        A = -MU_0 * I / (np.pi * a**2) * r_tau1**2 / 4 + MU_0 * I / (2 * np.pi) * np.log(r_tau2 / a) + MU_0 * I / (4 * np.pi)
    elif r_tau2 <= a and r_tau1 >= a:
        A = MU_0 * I / (np.pi * a**2) * r_tau2**2 / 4 - MU_0 * I / (2 * np.pi) * np.log(r_tau1 / a) - MU_0 * I / (4 * np.pi)
    elif r_tau1 >= a and r_tau2 >= a:
        A = MU_0 * I / (2 * np.pi) * np.log(r_tau2 / r_tau1)
    else:
        A = 0
    return A

def calculate_b_field(A, dx, dy):
    """
    Calculate the magnetic B field components from the vector potential.
    """
    Bx, By = np.gradient(A, dx, dy)
    Bx = -Bx
    return Bx, By

def center_window(fig):
    """
    Center the plot window on the screen.
    """
    root = Tk()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.withdraw()

    window_width, window_height = fig.get_size_inches() * fig.dpi
    position = (
        (screen_width - window_width) / 2,
        (screen_height - window_height) / 2
    )
    fig.canvas.manager.window.wm_geometry("+%d+%d" % position)

def plot_vector_potential_surface(A_normalized, X, Y):
    """
    Plot the surface of the normalized vector potential.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, A_normalized, cmap='viridis', edgecolor='none')
    ax.set_zlim(np.min(A_normalized), np.max(A_normalized))  
    ax.set_title('Normalized Vector Potential Surface')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('Az (Normalized)')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    center_window(fig)
    plt.show()

def plot_equipotential_lines(A_normalized, X, Y):
    """
    Plot the equipotential lines of the normalized vector potential.
    """
    fig, ax = plt.subplots(figsize=(10, 10))    
    levels = np.linspace(np.min(A_normalized), np.max(A_normalized), 40)
    contourf = ax.contourf(X, Y, A_normalized, levels=levels, cmap='viridis', extend='both')
    plt.colorbar(contourf, ax=ax, shrink=0.5, aspect=5)
    contour = ax.contour(X, Y, A_normalized, levels=levels, colors='black', linestyles='solid')
    ax.clabel(contour, inline=True, fontsize=8, fmt="%.2f")
    ax.set_title('Equipotential Lines')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(True)
    ax.set_aspect('equal', adjustable='box')
    center_window(fig)
    plt.show()

def plot_b_field(X, Y, Bx, By):
    """
    Plot the magnetic field lines using streamlines.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    stream = ax.streamplot(X, Y, Bx, By, cmap='viridis', linewidth=1, density=1)
    ax.set_title('Magnetic Field Lines')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(True)
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    ax.set_aspect('equal', adjustable='box')
    center_window(fig)
    plt.tight_layout() 
    plt.show()

def generate_meshgrid(x_range, y_range):
    """
    Generate a meshgrid for given x and y ranges.
    """
    x = np.linspace(x_range[0], x_range[1], 200)
    y = np.linspace(y_range[0], y_range[1], 200)
    return np.meshgrid(x, y)

def main():
    x_range = (-Lx / 2, Lx / 2)
    y_range = (-Ly / 2, Ly / 2)
    X, Y = generate_meshgrid(x_range, y_range)
    A = np.array([[vector_potential(x, y, a, d, I) for x in X[0]] for y in Y[:, 0]])
    A_normalized = A / (MU_0 * I / (2 * np.pi))
    dx = x_range[1] - x_range[0]
    dy = y_range[1] - y_range[0]
    Bx, By = calculate_b_field(A_normalized, dx, dy)
    plot_b_field(X, Y, Bx, By)
    plot_vector_potential_surface(A_normalized, X, Y)
    plot_equipotential_lines(A_normalized, X, Y)

if __name__ == "__main__":
    main()

