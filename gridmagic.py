# This script serves the reason to create an evenly spaced grid, map logarithmically spaced data onto the grid, interpolate the empty grid cells and create a streamline plot of the data.
# This is done by setting an evenly spaced grid onto the data and computing the nearest grid intersection point for each log-spaced data point (weighted mean values for more than one point.) This 'new' data can then be fed into matplotlib routines for streamline plots.

# imports
import numpy as __np__
import matplotlib.pyplot as __plt__
import matplotlib.mlab as __mlab__

def BaseGrid2D(xdata, ydata, ncells=50):
    """
    This function uses the given xdata and ydata to create an evenly spaced x and y list with borders similar to the given data. The number of cells controls the smoothness of the grid. Iterating over the lists resembles a grid.

    returns: two 1dim arrays for x and y
    """

    x_bounds = [__np__.min(xdata), __np__.max(xdata)]  # lower and upper limit for the x axis of the grid
    y_bounds = [__np__.max(ydata), __np__.min(ydata)]  # same for the y axis of the grid

    return __np__.linspace(*x_bounds, num=ncells), __np__.linspace(*y_bounds, num=ncells)

def FindMinVals(vals, quantity, offset=0, absolute=False):
    """
    This function finds the smallest values in a given array 'vals'. The amount of values to be found is specified by quantity. The offset is subtracted from the given values. If abs is set True, the absolute value of vals-offset is used to find the minimum. This can be useful, if you're searching for the nearest value(s) around a specific fixed value.

    returns array of the minimum values
    """
    # if an offset is given, subtract it:
    if(offset):
        if(absolute):
            vals_abs = __np__.abs(vals-offset)
        else:
            vals = vals-offset

    # check, if absolute value is wanted:
    if(not absolute):
        # check the lenght of the given array
        if(len(vals) < quantity):
            return __np__.sort(vals, kind='mergesort')

        else:
            min_vals = []  # initialize empty list

            for i in range(0, quantity):  # find the smallest <quantity> items and append them to the vals list
                min_index = vals.argmin()
                min_vals.append(vals[min_index])
                vals = __np__.delete(vals, min_index)

            # typecast the list to an array
            return __np__.array(min_vals)

    if(absolute):
        # check the lenght of the given array
        if(len(vals) < quantity):
            return __np__.sort(vals, kind='mergesort')

        else:
            min_vals = []  # initialize empty list

            for i in range(0, quantity):  # find the smallest <quantity> items and append them to the vals list
                min_index = vals_abs.argmin()
                min_vals.append(vals[min_index])
                vals_abs = __np__.delete(vals_abs, min_index)

            # typecast the list to an array
            return __np__.array(min_vals)

# def NearestDataPoints2D(data, xgrid, ygrid, Ndata):
#     """
#     This function takes a 2-dim data array, array([xdata, ydata]), and two 1-dim grid-coordinate-lists to map the data on. The nearest data point for a grid intersection point is found via the FindMinVals function and referencing the corresponding indices to the x,y data value pairs. Ndata is the amount of values that should be in the new gridcell.
#     """
#
#     #xdata, ydata = data  # extracting the x and y arrays
#     #sorted_data = __np__.sort(data, axis=0, kind='mergesort')  # sorting the data along the x-axis
#
#     # iterating over the grid lists:
#     for i in xgrid:
#         for j in ygrid:
#             nearest_points = FindMinVals(data[0], )
