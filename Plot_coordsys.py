import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

def plot_coordinatframe(X, Y, Z, Q, titel=None, crystalframe = None, crystalframe2 = None, I=None, rotation = None):

    fig = plt.figure(titel)


    ax = fig.gca(projection='3d')

    # plot Specimen frame
    ax.plot([0,1], [0,0], 0, label='Specimen frame', color='blue')
    ax.text(0.8, 0.1, 0, "x_P", color='blue')

    ax.plot([0,0], [0,1], 0, color='blue')
    ax.text(0.1, 0.8, 0, "y_P", color='blue')

    ax.plot([0,0], [0,0], [0,1],  color='blue')
    ax.text(0, 0.1, 0.8, "z_P", color='blue')

    # plot Rotated frame
    ax.plot([0,X[0]], [0,X[1]], [0, X[2]], label='rotated frame', color='red')
    ax.text(X[0]*0.8, X[1]*0.8, X[2]*0.8, "x_L", color='red')

    ax.plot([0,Y[0]], [0,Y[1]], [0, Y[2]], color='red')
    ax.text(Y[0]*0.8, Y[1]*0.8, Y[2]*0.8, "y_L", color='red')

    ax.plot([0,Z[0]], [0,Z[1]], [0, Z[2]], color='red')
    ax.text(Z[0]*0.8, Z[1]*0.8, Z[2]*0.8, "z_L", color='red')

    # q
    ax.plot([0,Q[0]], [0,Q[1]], [0, Q[2]],label='scattering vector', color='green')
    ax.text(Q[0], Q[1], Q[2], "Q", color='green')

    # plot crystallframe frame
    if (crystalframe is not None):
        X, Y, Z, rot = crystalframe
        ax.plot([0, X[0]], [0, X[1]], [0, X[2]], label='crystal frame', color='magenta')
        ax.text(X[0]*0.6, X[1]*0.6, X[2]*0.6, "x_C, %.0f"%(rot), color='magenta')

        ax.plot([0, Y[0]], [0, Y[1]], [0, Y[2]], color='magenta')
        ax.text(Y[0]*0.6, Y[1]*0.6, Y[2]*0.6, "y_C, %.0f"%(rot), color='magenta')

        ax.plot([0, Z[0]], [0, Z[1]], [0, Z[2]], color='magenta')
        ax.text(Z[0]*0.6, Z[1]*0.6, Z[2]*0.6, "z_C, %.0f"%(rot), color='magenta')

    if (crystalframe2 is not None):
        X, Y, Z, rot = crystalframe2

        ax.plot([0, X[0]], [0, X[1]], [0, X[2]], label='Euler frame', color='black')
        ax.text(X[0]*0.5, X[1]*0.5, X[2]*0.5, "x_E, %.0f"%(rot), color='black')

        ax.plot([0, Y[0]], [0, Y[1]], [0, Y[2]], color='black')
        ax.text(Y[0]*0.5, Y[1]*0.5, Y[2]*0.5, "y_E, %.0f"%(rot), color='black')

        ax.plot([0, Z[0]], [0, Z[1]], [0, Z[2]], color='black')
        ax.text(Z[0]*0.5, Z[1]*0.5, Z[2]*0.5, "z_E, %.0f"%(rot), color='black')

    if (I is not None):
        X = I
        ax.plot([0, X[0]], [0, X[1]], [0, X[2]], label='schnittgerade frame', color='darkblue')
        ax.text(X[0]*1, X[1]*1, X[2]*1, "I, %.0f"%(rotation), color='darkblue')

    # ax.arrow(0, 0, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)

    ax.legend()

# X = np.array([[1], [-1], [-1]])
# Y = np.array([[-1], [1], [-1]])
# Z = np.array([[-1], [-1], [1]])
# plot_coordinatframe(X, Y, Z)
# plt.show()
