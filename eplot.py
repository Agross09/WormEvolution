# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sim

# Graphs Vmem vs time. 'Which_cells' says which cells to graph.
# Really, it isn't limited to Vmem; it can do any piece of data that's
# stored in a 1D array[N_CELLS])
# Inputs:
#   - t_shots is a list of times
#   - data_shots is a list (of the same size as t_shots) of items. Each item
#     corresponds to a time in t_shots, and is the value of something at that
#     time. It is a 1D Numpy array. It is often an array[N_CELLS] of cell Vmem,
#     or a single row from cc_cells (which is also an array[N_CELLS]).
#   - which_cells is an array of indices, saying which indices from each item
#     in data_shots to actually plot. In the common case when the items are
#     [N_CELLS], it tells which cells to plot.
#   - ylabel is the graph label for the y axis (the x axis label is "Time(s)").
#   - filename is an optional argument. If empty, then the graph is drawn to
#     the screen. If provided, then it should be a string that is a filename,
#     and the graph is written to that file. Be sure that the file has a
#     reasonable file extension (e.g., ".jpg"), so that plot_Vmem_graph() knows
#     what file format to use.
def plot_Vmem_graph (t_shots, data_shots, which_cells, ylabel, filename=None):
    plt.figure()

    for cell in which_cells:
        data = np.asarray([shot[cell] for shot in data_shots])
        print ('Data for cell #',cell, ' is in [',data.min(),':',data.max(),']')
        plt.plot(t_shots, data, linewidth=2.0, label=str(cell))

    plt.xlabel('Time (s)', fontsize = 20)
    plt.ylabel(ylabel)

    # The legend is the list of which cells we're graphing & which color
    # line each uses.
    leg = plt.legend(loc='best',ncol=2,mode="expand",shadow=True,fancybox=True)
    leg.get_frame().set_alpha(0.5)

    if (filename is None):
        plt.show()
    else:
        plt.savefig(filename)

######################################################
# Pretty plotting.
# A pretty plot is a pretty representation of the network of cells at one point
# in time.
# - Each cell is drawn as a circle on the screen, with the cell's Vmem
#   written inside. Furthermore, the cells are colored according to their Vmem
#   (like a heat map).
# - GJs are drawn as lines between the appropriate cells.
######################################################

# Set the shape of a future output plot of the cell network.
# Inputs: shape
#	'shape' is a two-element array that sets the shape of an eventual
#	pretty plot. Each cell will be drawn as a circle; the circles will be
#	in a rectangular grid, with the GJs draw as lines connecting the cells.
#	Shape[0] tells how many rows in the plot; [1] tells how many columns.
#	With two rows of three columns, the cells will be drawn as
#	3 4 5
#	0 1 2
# Outputs: none, but the side effect of saving shape for a later plot call.
g_shape = None
def set_network_shape (shape):
    global g_shape
    g_shape = shape

# Draw a network with each cell as a circle, and each GJ as a line connecting
# its two cells.
# Label each cell with its index
# Color each cell based on whatever data we're given.
# Inputs: data[num_cells]
#	A 1D array of numbers, one per cell, saying what to plot for each cell.
def pretty_plot (data):
    global g_shape
    num_layers = g_shape[0]
    cells_per_layer = g_shape[1]

    # Assign a plot-worthy x,y coordinate pair to each cell.
    # Specifically, build xypts[N_CELLS,2]: each row #r is the (x,y) coordinates
    # of where to plot cell #r. Each layer of cells is a row in the plot, and
    # row #0 is at the bottom (with cell #0 at the left).
    # So if there are 2 cells per layer and 3 layers, then xypts is
    #	[[0. 0.] [1. 0.] [0. 1.] [1. 1.] [0. 2.] [1. 2.] ]
    # I.e., cell #0 is drawn at the lower left; then go left to right across
    # the bottom row, then left to right one row up, etc.
    # And num_layers(cells_per_layer) is the number of rows(column)
    y = np.linspace(0, num_layers-1, num_layers)
    x = np.linspace(0, cells_per_layer-1, cells_per_layer)
    X, Y = np.meshgrid(x, y)
    xypts = np.asarray([X.ravel(), Y.ravel()]).T

    # Line segments defining gap junctions (for plotting)
    # After this fancy indexing, gj_segs is a 3D array [N_GJ,2,2]
    # Any piece gj_segs[g] is a 2x2 array that describes how to draw GJ #g
    # as a line segment. I.e., the 1st row of the 2x2 array is the (x,y)
    # location of the cell for the GJ's input, and the 2nd is the (x,y) location
    # of the cell for the GJ's output.
    gj_from_to = np.stack((sim.gj_connects['from'],sim.gj_connects['to']),1)
    gj_segs = xypts[gj_from_to]

    plt.figure()
    ax = plt.subplot(111)

    # Each cell is a circle:
    # - 'c' gives the colors; plt.colorbar() picks a color for each element of
    #   Vm (i.e., maps from our Vmem range to nice colors)
    # - 's' is the area of the circles.
    plt.scatter(xypts[:,0], xypts[:,1], c=data, s=500)
    # Draw a bar showing the mapping from Vmem to color
    plt.colorbar(orientation='horizontal')

    # Label each cell with its index number.
    for i, (xi, yi) in enumerate(xypts):
        label = "Cell #{:d}\n({:.2f})".format (i, (data[i]))
        print (label)
        plt.text(xi, yi, label, fontsize=14, fontweight='bold')

    # Draw the gap junctions
    GJs = LineCollection(gj_segs)
    ax.add_collection(GJs)

    plt.axis('equal')	# ensure that circles are circular and not ellipses.
    plt.axis('off')
    plt.show()
