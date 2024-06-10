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
    energy_particle = data[:, 4]*0.5 #MeV
    
    ax.scatter(x, y, energy_particle, c='r', marker='o')

    ax.set_xlabel('X ')
    ax.set_ylabel('Y')
    ax.set_zlabel('Energy, MeV')

    plt.show()


if __name__ == "__main__":
    # Absolute or relative path to the .dat file
    file_path = '../TEST_PARTICLE_DATA/tp_data_0000.dat'
    
    # Make sure the path exists
    if not os.path.exists(file_path):
        print(f"The file {file_path} does not exist.")
    else:
        data = load_data(file_path)
        create_3D_color_plot(data)
