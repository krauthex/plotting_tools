####################################################################################################
###     This script is a collection of functions used for a specific kind of data analysis,      ###
###     mainly making a coordinate system conversion and doing plots.                            ###
###                                                                                              ###
###     Author: Fabian Krautgasser, krautgasser@mpia-hd.mpg.de      2016                         ###
###     Python Version: >2.*                                                                     ###
####################################################################################################
# TODO: cleanup the code in the plotting functions
# TODO: add the option to plot the mass flux

# imports
import numpy as __np__
import matplotlib.pyplot as __plt__
import matplotlib.mlab as __mlab__
import datetime as __dt__

"""
parameter dictionary:
r_hill is the radius to find boxes on, Here it is the Hill Sphere for the given planet radius
delta is the factor around the r_hill to find boxes in
box is the radial size of a simulation box
r0_planet is the physical radius of the planet
planet_location_x is the location of the planet in (cartesian) x-direction in the Star-centered-Coordinate-System
"""
params = dict(r_hill=16685629566.462387, delta=0.1, box=1.4e9, r0_planet=2*6.371E8, planet_location_x=-748000011879.0)


def timestamp():
        """
        returns: a string of the current time in the format YYYY-MM-DD_HH:MM:SS
        """
        return str(__dt__.datetime.now()).split(sep='.')[0].replace(' ', '_')  # timestamp


def coordinate_trafo(r_star, theta_star, phi_star, plan_loc_x=params['planet_location_x']):
    """
    This function transforms the given arrays of the spherical coordinates in the Star-centered-Coordinate-System
    to spherical coordinates in the Planet-centered-Coordinate-System.
    Spherical_star -> cartesian_star -> cartesian_planet -> spherical_planet

    returns: r_pl, theta_pl, phi_pl
    """
    # computing the cartesian coordinates in the star coordinate system
    x_star = []
    y_star = []
    z_star = []

    for i in r_star:
        for j in theta_star:
            for k in phi_star:
                x_star.append(i*__np__.sin(j)*__np__.cos(k))
                y_star.append(i*__np__.sin(j)*__np__.sin(k))
                z_star.append(i*__np__.cos(j))

    x_star = __np__.array(x_star)
    y_star = __np__.array(y_star)
    z_star = __np__.array(z_star)

    # cartesian coordinates in the planet coordinate system
    x_plan = x_star - plan_loc_x
    y_plan = y_star
    z_plan = z_star

    # computing spherical coordinates in the planet coordinate system
    r_plan = __np__.sqrt(x_plan**2 + y_plan**2 + z_plan**2)
    theta_plan = __np__.arccos(z_plan / r_plan)
    phi_plan = __np__.arctan2(y_plan, x_plan)

    return r_plan, theta_plan, phi_plan


def velocity_trafo(vr_star, vtheta_star, vphi_star, r_star, theta_star, phi_star, r_plan, theta_plan, phi_plan):
    """
    This function transforms the given arrays of the velocity-vector-components in the star-centered-coordinate-system
    to the velocity-vector-components in the Planet-centered-Coordinate-System.

    returns: vr_plan, vtheta_plan, vphi_plan
    """
    coords = []
    for i in r_star:
        for j in theta_star:
            for k in phi_star:
                coords.append([i,j,k])

    coords = __np__.array(coords).T  # first entry are radii, second theta, third phi
    r_star, theta_star, phi_star = coords

    # transforming the velocities in the star-centered coordinate system into cartesian coordinates x,y,z
    x_dot_star = (vr_star*__np__.sin(theta_star)*__np__.cos(phi_star) +
                  vtheta_star*__np__.cos(theta_star)*__np__.cos(phi_star) -
                  vphi_star*__np__.sin(phi_star))

    y_dot_star = (vr_star*__np__.sin(theta_star)*__np__.sin(phi_star) +
                  vtheta_star*__np__.cos(theta_star)*__np__.sin(phi_star) +
                  vphi_star*__np__.cos(phi_star))

    z_dot_star = vr_star*__np__.cos(theta_star) - vtheta_star*__np__.sin(theta_star)

    # transforming the velocities in the (cartesian) planets coordinate system into spherical coordinates r, teta, phi
    vr_plan = (x_dot_star*__np__.sin(theta_plan)*__np__.cos(phi_plan) +
               y_dot_star*__np__.sin(theta_plan)*__np__.sin(phi_plan) +
               z_dot_star*__np__.cos(theta_plan))

    vtheta_plan = (x_dot_star*__np__.cos(theta_plan)*__np__.cos(phi_plan) +
                   y_dot_star*__np__.cos(theta_plan)*__np__.sin(phi_plan) -
                   z_dot_star*__np__.sin(theta_plan))

    vphi_plan = (y_dot_star*__np__.cos(phi_plan) -
                 x_dot_star*__np__.sin(phi_plan))

    return vr_plan, vtheta_plan, vphi_plan


def extract_shell(vr, vtheta, vphi, r_pl, shell=params['r_hill'], delta=params['delta'], box=params['box']):
    """
    This function takes the values of vr, vtheta, vphi and extracts the velocity vectors who reside
    in the given shell +/- delta*box

    returns: vr, vtheta, vphi (extracted)
    """
    # NOTE: find a good value for delta, so that the plots are looking nice
    indices = [i for i,j in enumerate(r_pl) if ((shell + delta*box) > j > (shell- delta*box))]
    print('\n:: Shell Extractor:\n:: Number of Values in shell:', len(indices), '\n')

    return __np__.array([vr[indices], vtheta[indices], vphi[indices]]), indices


def hammer_projection(theta, phi):
    """
    theta, phi should be given in radians
    converts phi, theta to latitude and longitude (in rad)
    does a hammer projection from lat and lon

    returns: x, y coordinates
    """
    # lat, lon calculation
    lat = phi  # -pi <= phi <= pi
    lon = theta - __np__.pi/2  # -pi/2 <= theta <= pi/2

    # hammer projection
    x = 2*__np__.sqrt(2)*__np__.cos(lon)*__np__.sin(lat/2) / __np__.sqrt(1+__np__.cos(lon)*__np__.cos(lat/2))
    y = __np__.sqrt(2)*__np__.sin(lon) / __np__.sqrt(1+__np__.cos(lon)*__np__.cos(lat/2))

    return x,y


def contour_plot_velocity(x, y, vr, num_levels=15, levels=[], colorbar=True, figsize=(12,9)):
    """
    This function takes the values of x (phi) and y (theta) (for example from the hammer projection) and
    plots a contour plot of the given values for vr. This is done by a weighted 2D Histogram. The
    Points in between are interpolated linearly.
    if levels are given, num_levels is omitted.

    returns: the figure and the axis objects
    """
    # defining a x,y-grid to interpolate z-data on ([-pi, pi] and [-pi/2, pi/2])
    xg = __np__.linspace(-4, 4, 1000)
    yg = __np__.linspace(-2, 2, 1000)
    zg = __mlab__.griddata(x, y, vr, xg, yg, 'linear')  # interpolate the z-data (vr)

    # Plotting Section
    color_map = __plt__.cm.seismic
    fig = __plt__.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Contour Plot
    if(len(levels)):
        # if levels are given, use them
        con = ax.contourf(xg, yg, zg, levels, cmap=color_map, vmax=abs(zg).max(), vmin=-abs(zg).max())
    else:
        # if no levels are given, use num_levels. levels are chosen automatically.
        con = ax.contourf(xg, yg, zg, num_levels, cmap=color_map, vmax=abs(zg).max(), vmin=-abs(zg).max())

    # Color Bar
    cax = fig.add_axes([0.8, 0.25, 0.025, 0.5])  # creating a new axis for the colorbar
    cb = __plt__.colorbar(mappable=con, cax=cax, orientation='vertical')
    cb.formatter.set_powerlimits((0,0))  # setting the scientific notation of the color bar ticks
    cb.update_ticks()
    cb.set_label('$v_R$', rotation='horizontal')

    return fig, ax, cax


def scatter_contour(x, quantity, x_norm=params['r0_planet'], nbins=80, num_levels=15, levels=[],
                    colorbar=True, figsize=(12,9), xlog=True):
    """
    This function takes a quantity (density, temperature) and does a scatter plot and a dot-density
    contour plot on top. This is done via a 2D-Histogram. The x-Axis is normalized by the factor x_norm.
    if levels are given, num_levels is omitted.

    returns: fig, ax, cax
    """
    # normalizing x
    x = x/x_norm

    # depending on xlog, computing the 2d histogram
    if(xlog):
        H, xedges, yedges = __np__.histogram2d(__np__.log(x), quantity, bins=nbins)
        xcenter = __np__.exp((xedges[:-1] + xedges[1:]) / 2)
        ycenter = (yedges[:-1] + yedges[1:]) / 2
    else:
        H, xedges, yedges = __np__.histogram2d(x, quantity, bins=nbins)
        xcenter = (xedges[:-1] + xedges[1:]) / 2
        ycenter = (yedges[:-1] + yedges[1:]) / 2

    # the histogram has to be flipped and rotated
    H = __np__.rot90(H)
    H = __np__.flipud(H)

    # Plotting Section
    color_map = __plt__.cm.seismic
    fig = __plt__.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Dot Plot
    ax.plot(x, quantity, 'g.', alpha=0.3, zorder=0, markeredgewidth=0.0)
    if(xlog):
        ax.set_xscale('log')
    else: ax.set_xscale('linear')

    # Contour Plot
    if(len(levels)):
        # if levels are given, use them
        con = ax.contour(xcenter, ycenter, H, levels, cmap=color_map, zorder=1, linewidths=2)
    else:
        # if no levels are given, use num_levels. levels are chosen automatically.
        con = ax.contour(xcenter, ycenter, H, num_levels, cmap=color_map, zorder=1, linewidths=2)

    # Color bar
    cax = fig.add_axes([0.8, 0.25, 0.025, 0.5])  # creating a new axis for the colorbar
    cb = __plt__.colorbar(mappable=con, cax=cax, orientation='vertical')
    cb.formatter.set_powerlimits((0,0))  # setting the scientific notation of the color bar ticks
    cb.update_ticks()
    cb.set_label('$Dot\ Density$')

    return fig, ax, cax
