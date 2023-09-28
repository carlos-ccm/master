import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# load the coordinates from file
data_coords = np.loadtxt('coords.txt')
x = data_coords[:,0]
y = data_coords[:,1]
z = data_coords[:,2]

# load the coefficients from file
coeff = np.loadtxt('coeffs.txt')

# iterate over the columns of the coefficient array and create a plot for each one
for i in range(coeff.shape[1]):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x, y, z, color='black')

    # create the molecular orbitals as spheres with varying size and orientation
    for j in range(len(x)):
        r = abs(coeff[j,i])
        if coeff[j,i] >= 0:
            ax.plot([x[j]], [y[j]], [z[j]+0.1*r], marker='o', markersize=50*r, color='lightblue')
        else:
            ax.plot([x[j]], [y[j]], [z[j]-0.1*r], marker='o', markersize=50*r, color='grey')

    # connect the atoms with lines based on their distance
    for j in range(len(x)):
        for k in range(j+1, len(x)):
            dist = np.sqrt((x[j]-x[k])**2 + (y[j]-y[k])**2 + (z[j]-z[k])**2)
            if dist <= 2:
                ax.plot([x[j], x[k]], [y[j], y[k]], [z[j], z[k]], color='black', linewidth=1)

    # set the axis limits
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    min_z, max_z = np.min(z), np.max(z)
    if min_z == max_z:
        min_z -= 0.5
        max_z += 0.5
    ax.set_zlim(min_z, max_z)

    # hide the axis
    ax.set_axis_off()

    # set the title
    ax.set_title('Ψ {}'.format(i+1))

    # show the plot
    plt.show()
