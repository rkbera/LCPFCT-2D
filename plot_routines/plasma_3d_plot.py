import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import os

# Load data from .dat file
def load_data(file_path):
    data = np.loadtxt(file_path)
    return data

# Create a 3D color plot
def create_3D_color_plot(data):
     # Apply a style
    plt.style.use('classic')

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Assuming data has three columns: x, y, z
    x = data[:, 0]
    y = data[:, 1]
    density = data[:, 2]
    density_beam = data[:, 3]
    Ex = data[:, 4]


    # Create a grid for the surface plot
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)
    density = griddata((x, y), density, (X, Y), method='cubic')
    density_beam = griddata((x, y), density_beam, (X, Y), method='cubic')
    Ex = griddata((x, y), Ex, (X, Y), method='cubic')


    # Create a surface plot with color based on z-values
    surf = ax.plot_surface(X, Y, density, cmap='plasma', edgecolor='none')
    
    # Add color bar which maps values to colors
    cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
    cbar.set_label('density', rotation=270, labelpad=15)

    # Enhancing plot aesthetics
    #ax.set_title('Beautiful 3D Color Plot', fontsize=15)
    ax.set_xlabel('X ', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('density', fontsize=12)

    # Adjusting the view angle for better visualization
    ax.view_init(elev=20, azim=-60)

    # Adding grid and customizing ticks
    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)

    plt.show()


if __name__ == "__main__":
    # Absolute or relative path to the .dat file
    file_path = '../PLASMA_DATA/plasma_data_0100.dat'
    
    # Make sure the path exists
    if not os.path.exists(file_path):
        print(f"The file {file_path} does not exist.")
    else:
        data = load_data(file_path)
        create_3D_color_plot(data)
