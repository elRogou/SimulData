import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
import plotly.graph_objs as go

def plot_3d_density_dependent(x, y, z):
    """
    Plots the joint Gaussian density of three arrays in 3D using a contour plot,
    assuming a dependency between the variables.

    Args:
        x, y, z (array-like): The arrays to plot.

    Returns:
        None
    """
    # Combine the arrays into a single matrix
    data = np.vstack((x, y, z)).T

    # Calculate the multivariate Gaussian density
    kde = gaussian_kde(data.T)

    # Create a grid of points to plot the density on
    x_grid, y_grid, z_grid = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
    pos = np.empty(x_grid.shape + (3,))
    pos[:, :, :, 0] = x_grid
    pos[:, :, :, 1] = y_grid
    pos[:, :, :, 2] = z_grid

    # Calculate the density for each point in the grid
    data_eval = np.vstack((x_grid.flatten(), y_grid.flatten(), z_grid.flatten()))
    density = kde(data_eval)

    # Reshape the density to match the grid
    density = density.reshape(x_grid.shape)

    # Set up the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the density as a contour plot
    ax.contour(x_grid, y_grid, density, cmap='viridis')

    # Add labels and a title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Density')
    ax.set_title('Gaussian Density in 3D')

    # Show the plot
    plt.show()


def plot_3d_density_independent(x, y, z):
    """
    Plots the joint Gaussian density of three arrays in 3D using a contour plot.

    Args:
        x, y, z (array-like): The arrays to plot.

    Returns:
        None
    """
    # Calculate the Gaussian density for each dimension
    kde_x = gaussian_kde(x)
    kde_y = gaussian_kde(y)
    kde_z = gaussian_kde(z)

    # Create a grid of points to plot the density on
    x_grid, y_grid, z_grid = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
    pos = np.empty(x_grid.shape + (3,))
    pos[:, :, :, 0] = x_grid
    pos[:, :, :, 1] = y_grid
    pos[:, :, :, 2] = z_grid

    # Calculate the density for each point in the grid
    density = kde_x(x_grid.flatten()) * kde_y(y_grid.flatten()) * kde_z(z_grid.flatten())
    density = density.reshape(x_grid.shape)

    # Create a 3D contour plot
    fig = go.Figure(data=go.Contour(x=x_grid.flatten(), y=y_grid.flatten(), z=density, contours=dict(z=dict(show=True, start=density.min(), end=density.max(), size=0.1)), colorscale='Viridis', showscale=True))

    # Add labels and a title
    fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Density',), title='Gaussian Density in 3D')

    # Show the plot
    fig.show()

def plot_3d_surface_independent(x, y, z):
    """
    Plots the joint Gaussian density of three arrays in 3D using a surface plot.

    Args:
        x, y, z (array-like): The arrays to plot.

    Returns:
        None
    """
    # Calculate the Gaussian density for each dimension
    kde_x = gaussian_kde(x)
    kde_y = gaussian_kde(y)
    kde_z = gaussian_kde(z)
    print("KDEs:")
    print("kde_x:", kde_x)
    print("kde_y:", kde_y)
    print("kde_z:", kde_z)

    # Create a grid of points to plot the density on
    x_grid, y_grid, z_grid = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
    print(f"x_grid shape: {x_grid.shape}")
    print(f"y_grid shape: {y_grid.shape}")
    print(f"z_grid shape: {z_grid.shape}")

    pos = np.empty(x_grid.shape + (3,))
    pos[:, :, :, 0] = x_grid
    pos[:, :, :, 1] = y_grid
    pos[:, :, :, 2] = z_grid

    # Calculate the density for each point in the grid
    density = kde_x(x_grid.flatten()) * kde_y(y_grid.flatten()) * kde_z(z_grid.flatten())
    density = density.reshape(x_grid.shape)
    print("Density:")
    print("density:", density)

    # Create a 3D surface plot
    fig = go.Figure(data=go.Surface(x=x_grid, y=y_grid, z=density, colorscale='Viridis'))

    # Add labels and a title
    fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Density',), title='Gaussian Density in 3D')

    # Show the plot
    fig.show()

    # return fig