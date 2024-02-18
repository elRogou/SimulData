from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import plotly.graph_objs as go

# Here is a function that gets an x,y and z variables as input,
# calculates gaussian density and plots it in 3D.
# the calculation of density is done using the gaussian_kde function from scipy.stats
# with an assumption of dependency between the variables.

def plot_3d_density_dependent(x, y, z):
    # combine the arrays into a single matrix
    data = np.vstack((x, y, z)).T

    # calculate the multivariate gaussian density
    kde = gaussian_kde(data.T)

    # create a grid of points to plot the density on
    x_grid, y_grid, z_grid = np.mgrid[x.min():x.max():100j, y.min():y.max():100j, z.min():z.max():100j]
    pos = np.empty(x_grid.shape + (3,))
    pos[:, :, :, 0] = x_grid
    pos[:, :, :, 1] = y_grid
    pos[:, :, :, 2] = z_grid

    # calculate the density for each point in the grid
    data_eval = np.vstack((x_grid.flatten(), y_grid.flatten(), z_grid.flatten()))
    density = kde(data_eval)

    # reshape the density to match the grid
    density = density.reshape(x_grid.shape)

    # set up the 3D plot
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # plot the density as a contour plot
    # ax.contour(x_grid, y_grid, density, cmap='viridis')

    # add labels and a title
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Density')
    # ax.set_title('Gaussian Density in 3D')

    # show the plot
    # plt.show()

    # plot using plotly
    fig = go.Figure(data=[go.Surface(z=density)])
    fig.update_layout(title='Gaussian Density in 3D', autosize=False,
                      width=500, height=500,
                      margin=dict(l=65, r=50, b=65, t=90))
    fig.show()
    return density


def plot_3d_contour(x,y,z):
    # calculate the gaussian kde density
    kde_x = gaussian_kde(x)
    kde_y = gaussian_kde(y)
    kde_z = gaussian_kde(z)

    # create a grid of points
    x_grid = np.linspace(x.min(), x.max(), 100)
    y_grid = np.linspace(y.min(), y.max(), 100)
    z_grid = np.linspace(z.min(), z.max(), 100)

    # evaluate the gaussian kde on the grid
    x_density = kde_x(x_grid)
    y_density = kde_y(y_grid)
    z_density = kde_z(z_grid)

    # plot the results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_grid, y_density, z_density, color='blue')
    ax.plot(x_density, y_grid, z_density, color='red')
    ax.plot(x_density, y_density, z_grid, color='green')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Density')
    ax.set_title('Gaussian Density in 3D')
    plt.show()

    # plot the results with plotly
    fig = go.Figure(data=[go.Scatter3d(x=x_grid, y=y_density, z=z_density, mode='lines', line=dict(color='blue')),
                            go.Scatter3d(x=x_density, y=y_grid, z=z_density, mode='lines', line=dict(color='red')),
                            go.Scatter3d(x=x_density, y=y_density, z=z_grid, mode='lines', line=dict(color='green'))])
    fig.show()
