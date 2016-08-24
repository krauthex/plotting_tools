# This script serves the reason to create an evenly spaced grid, map logarithmically spaced data onto the grid, interpolate the empty grid cells and create a streamline plot of the data.
# This is done by setting an evenly spaced grid onto the data and computing the nearest grid intersection point for each log-spaced data point (weighted mean values for more than one point.) This 'new' data can then be fed into matplotlib routines for streamline plots.

# imports
import numpy as __np__
import matplotlib.pyplot as __plt__
# import matplotlib.mlab as __mlab__
from scipy.interpolate import interp1d as __interp1d__  # used in FindCellDensity function

def BaseGrid2D(xdata, ydata, ncells=50):
    """
    This function uses the given xdata and ydata to create an evenly spaced x and y list with borders similar to the given data. The number of cells controls the smoothness of the grid. Iterating over the lists resembles a grid.

    returns: four 1dim arrays: two for the x and y grid points and two for the x and y value points
    """

    x_bounds = [__np__.min(xdata), __np__.max(xdata)]  # lower and upper limit for the x axis of the grid
    y_bounds = [__np__.max(ydata), __np__.min(ydata)]  # same for the y axis of the grid

    # create the x and y grid indices
    x_grid_index = __np__.linspace(*x_bounds, num=ncells)
    y_grid_index = __np__.linspace(*y_bounds, num=ncells)

    # calculate the intermediate values between two grid indices as data value indices. also for x and y. This index lists contains (N-1)x(N-1) entries when using a NxN grid.
    x_vals_index = (x_grid_index[1:]+x_grid_index[:-1])/2
    y_vals_index = (y_grid_index[1:]+y_grid_index[:-1])/2

    return x_grid_index, y_grid_index, x_vals_index, y_vals_index

def FindMinVals(vals, quantity, offset=0, absolute=False):
    """
    This function finds the smallest values in a given array 'vals'. The amount of values to be found is specified by quantity. The offset is subtracted from the given values. If absolute is set True, the absolute value of (vals-offset) is used to find the minimum. This can be useful, if you're searching for the nearest value(s) around a specific fixed value.

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

# TODO:
#   - find a method, to calculate the mean $value for a given value grid point with respect to the surrounding grid cells:
#       * Method A: row/column 'density' and calculate mean value of both and interpolate/extrapolate missing grid cells by using the previously found values.
#       * Method B: sort values along x (or y) axis. use FindMinVals to find e.g. the 100 smallest values. filter out all values who's y (or x) value is too far away from the grid point. If the amount of filtered values is too low, rerun the FindMinVals function with a higher quantity.
#       * Method C: calculate r_ij = sqrt((x_data-i)^2 + (y_data-j)^2) and (merge)sort the r_ij values via FindMinVals. use e.g. the smallest 100 values (or up to a certain maximal radius).

# a class for raising errors
class ShapeMismatchError(Exception):
    def __init__(self, msg):
        self.message = msg

# Method A
def FindCellDensity(data, ncells, sort_axis=0, interpol_threshold=0):
    """
    This function takes a 2xN-dim data-array ([[coordinates], [values]]) and splits the data in (ncells-1) cells, after sorting them along the sort_axis. Then the mean value for each cell is calculated. If a cell contains less values than the interpolation threshold (interpol_threshold), the value for that cell is interpolated (or extrapolated, if there are no following cells that contain a value).

    Returns: (sorted) array of the calculated mean values and the interpolation function (or parameters?)
    """
    data = __np__.sort(data, axis=sort_axis, kind='mergesort')  # sorting the data

    # grid
    bounds = [__np__.min(data), __np__.max(data)]  # lower and upper limit for the x axis of the grid
    grid_index = __np__.linspace(*bounds, num=ncells)  # create the x grid indices
    vals_index = (grid_index[1:]+grid_index[:-1])/2  # intermediate values between two grid indices as data value indices.

    # initalize empty lists
    coordinate_block = []
    value_block = []
    means = []
    uncerts = []
    emptyCells = []

    # partition the data in blocks
    # the following part finds the index, where to place j in data[0] and from that position, the next grid cell edge. with these boundaries, the value blocks to calculate mean values from can be sliced from the original data.
    for i,j in enumerate(grid_index[:-1]):
        # data block boundaries
        lBound = __np__.searchsorted(data[0], j)
        uBound = __np__.searchsorted(data[0], grid_index[i+1])

        coordinate_block.append(data[0, lBound:uBound])
        value_block.append(data[1, lBound:uBound])

    for i,j in enumerate(value_block):
        if(len(j) > interpol_threshold):
            means.append(__np__.mean(j))
            uncerts.append(__np__.std(j)/__np__.sqrt(len(j)))  # calculating the uncertainty of the mean value

        else:
            emptyCells.append(i)  # appending the indices of empty cells to the list emptyCells

    # finding the indices, where two or more empty cells are next to each other by shifting the array by 1 and subtracting it from the original one. positions where empty cells are next to each other just differ by 1.
    emptyCells = __np__.array(emptyCells)
    emptyCellsShifted = __np__.roll(emptyCells, -1)
    diff = np.abs(emptyCells-emptyCellsShifted)

    # EdgeStretch !
    # expand the first value found next to the empty cells at the edges of the grid to all the empty edge cells
    # [][][][2][...] ---> [2][2][2][2][...]
    # 'left' side:
    if(emptyCells[0] == 0):
        emptyCellCounter = 1
        while(diff[emptyCellCounter-1] == 1):
            emptyCellCounter += 1

        emptyCells = emptyCells[emptyCellCounter:]  # shorten the list of empty cells on left edge
        MeanToExpand = means[0]
        UncertToExpand = uncerts[0]
        mean_to_add_left = []  # maybe exchange this part by : __np__.empty(emptyCellCounter); array.fill(MeanToExpand); and typecast the outcome to list
        uncert_to_add_left = []
        for i in range(emptyCellCounter):
            mean_to_add_left.append(MeanToExpand)
            uncert_to_add_left.append(UncertToExpand)
        means = mean_to_add_left + means
        uncerts = uncert_to_add_left + uncerts

    # 'right' side
    if(emptyCells[-1] == (len(value_block)-1)):
        rightDiff = __np__.abs(emptyCells-__np__.roll(emptyCells, 1))
        emptyCellCounter = 1
        while(rightDiff[-emptyCellCounter] == 1):
            emptyCellCounter += 1

        emptyCells = emptyCells[:-emptyCellCounter]  # shorten the list of empty cells on right edge
        MeanToExpand = means[-1]
        UncertToExpand = uncerts[-1]
        mean_to_add_right = []
        uncert_to_add_right = []
        for i in range(emptyCellCounter):
            mean_to_add_right.append(MeanToExpand)
            uncert_to_add_right.append(UncertToExpand)
        means = means+mean_to_add_right
        uncerts = uncerts+uncert_to_add_right

    # recalculate the 'diff' array:
    emptyCellsShifted = __np__.roll(emptyCells, -1)
    diff = np.abs(emptyCells-emptyCellsShifted)

    # Linear interpolation of empty cells between filled cells.
    # emptyCellPairCounter = 1
    for i,j in enumerate(diff):
        if(j == 1):
            emptyCellPairCounter = 1

            while(diff[i+emptyCellPairCounter] == 1):  # if there are more than one empty cells in a row, increase the counter
                emptyCellPairCounter = emptyCellPairCounter + 1

            # interp1d is piecewise linear
            LinInt = __interp1d__([vals_index[k] for k in [i, i+emptyCellPairCounter]], [means[l] for l in [i-1, i]])  # linear interpolation between the known values
            for k in emptyCells[i:i+emptyCellPairCounter+1]:
                means.insert(k, LinInt(vals_index[k]))
                uncerts.insert(k, __np__.mean([uncerts[l] for l in [k-1, k]]))

        else:
            LinInt = __interp1d__([vals_index[k] for k in [i, i+1]], [means[l] for l in [i-1, i]])
            means.insert(i, LintInt(vals_index[i]))
            uncerts.insert(i, __np__.mean([uncerts[l] for l in [k-1, k]]))
    # sanity checks:
    if(len(means) != len(value_index)):
        raise ShapeMismatchError("Number of mean values and value indices is not the same!")

    gridded_data = __np__.array([vals_index, means, uncerts])

    return gridded_data
