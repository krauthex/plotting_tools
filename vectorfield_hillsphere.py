# this testscript is used to slice the data around the planet
# to achive this, the cartesian coordinate data is used

import main
import numpy as __np__
import matplotlib.pyplot as __plt__

# data
print('\n:: Loading Data...\n')
vr, vtheta, vphi = __np__.loadtxt('data/prepared/VEL_PL', skiprows=1, unpack=True)  # velocities
r, theta, phi = __np__.loadtxt('data/prepared/SPH_COOR_PL.txt', skiprows=1, unpack=True)  # spherical coordinates
rho, temp = __np__.loadtxt('data/raw/RHO_TEMP.txt', skiprows=1, unpack=True)  # density and temperature

def cartesian_coord_planet(r_star, theta_star, phi_star, plan_loc_x=main.params['planet_location_x']):
    """
    This function calculates the cartesian coordinates of the planet's data from the given spherical coordinates
    of the star.

    returns x_plan, y_plan, z_plan
    """
    # computing the cartesian coordinates in the star coordinate system
    x_star = []
    y_star = []
    z_star = []

    # creating the matching value pairs
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

    return x_plan, y_plan, z_plan

# reading the raw data and calculating the cartesian coordinates for the planet
print("\n:: Reading in the spherical coordinates in the star's system.\n")
r_star = __np__.loadtxt('data/raw/GRID_I_COORD.txt', unpack=True, usecols=[1])
theta_star = __np__.loadtxt('data/raw/GRID_J_COORD.txt', unpack=True, usecols=[1])
phi_star = __np__.loadtxt('data/raw/GRID_K_COORD.txt', unpack=True, usecols=[1])
# print(__np__.shape(r_star), __np__.shape(theta_star), __np__.shape(phi_star))

print("\n:: Calculating the cartesian coordinates for the planet's coordinate system.\n")
PLANET_RADIUS = main.params['r0_planet']  # default radius
HILL_SPHERE = main.params['r_hill']  # radius of the hill sphere
x_plan, y_plan, z_plan = cartesian_coord_planet(r_star, theta_star, phi_star)

# slicing the data to a box shape around the planet with z = PLANET_RADIUS, x=y=2*HILL_SPHERE
print('\n:: Slicing the dataset.\n')
dataset = __np__.array([x_plan, y_plan, z_plan])
indices = [i for i,j in enumerate(dataset.T) if __np__.abs(j[0]) <= 2*HILL_SPHERE  # x-condition
                                             if __np__.abs(j[1]) <= 2*HILL_SPHERE  # y-condition
                                             if __np__.abs(j[2]) <= PLANET_RADIUS]  # z-condition
print('\n:: Number of values found: ', len(indices), '\n')

# using these indices to filter the needed values for the spherical coordinates and velocity components
sph_coord_sliced = __np__.array([r, theta, phi]).T[indices]
v_sliced = __np__.array([vr, vtheta, vphi]).T[indices]
dataset_sliced = dataset.T[indices]

# plotting these values should give and overview where the velocity values are found.
print('\n:: Ready to plot!\n')
fig =  __plt__.figure()
ax = fig.add_subplot(111)
ax.plot(dataset_sliced.T[0], dataset_sliced.T[1], 'k.', zorder=0)  # coordinates of values
hill = __plt__.Circle((0, 0), radius=HILL_SPHERE, color='r', fill=False, lw=2.0)  # Hillsphere
planet = __plt__.Circle((0, 0), radius=PLANET_RADIUS, color='b', fill=True)  # Drawing a dot for the planet
#outer_hill = __plt__.Circle((0, 0), radius=HILL_SPHERE+main.params['delta']*main.params['box'], color='r', fill=False, ls='dashed')
#inner_hill = __plt__.Circle((0, 0), radius=HILL_SPHERE-main.params['delta']*main.params['box'], color='r', fill=False, ls='dashed')
ax.add_artist(hill)
ax.add_artist(planet)
# the distance between the circles is so small, that they're indistinguishable from each other
#ax.add_artist(outer_hill)
#ax.add_artist(inner_hill)
fig.savefig('plots/sliced_values_test.pdf', format='pdf', dpi=200)
print('\n:: Saved the first plot, now for the second one...\n')

# Streamline/vector plot
def velocity_trafo_sph2car(r_star, theta_star, phi_star, vr_star, vtheta_star, vphi_star):
    """
    This function transforms the velocity-vector-components in the spherical coordinate system of the
    star to the planets cartesion coordinate system.

    returns: vx_plan, vy_plan, vz_plan
    """
    # computation of the velocity-vector-components
    coords = []
    for i in r_star:
        for j in theta_star:
            for k in phi_star:
                coords.append([i,j,k])

    coords = __np__.array(coords).T  # first entry are radii, second theta, third phi
    r_star, theta_star, phi_star = coords

    x_dot_star = (vr_star*__np__.sin(theta_star)*__np__.cos(phi_star) +
                  vtheta_star*__np__.cos(theta_star)*__np__.cos(phi_star) -
                  vphi_star*__np__.sin(phi_star))

    y_dot_star = (vr_star*__np__.sin(theta_star)*__np__.sin(phi_star) +
                  vtheta_star*__np__.cos(theta_star)*__np__.sin(phi_star) +
                  vphi_star*__np__.cos(phi_star))

    z_dot_star = vr_star*__np__.cos(theta_star) - vtheta_star*__np__.sin(theta_star)

    return x_dot_star, y_dot_star, z_dot_star

# load the values for vr_star, vtheta_star, vphi_star
vr_star, vtheta_star, vphi_star = __np__.loadtxt('data/raw/VEL_STAR.txt', unpack=True, skiprows=1)

# compute the velocity componets in cartesian coordinates and slice the array
print('\n:: Computing the velocity-vector-components in cartesian coordinates..\n')
vx_plan, vy_plan, vz_plan = velocity_trafo_sph2car(r_star, theta_star, phi_star, vr_star, vtheta_star, vphi_star)
vel_cart = __np__.array([vx_plan, vy_plan, vz_plan]).T[indices]

# do the vector plot
fig2 = __plt__.figure()
ax = fig2.add_subplot(111)
ax.quiver(dataset_sliced.T[0], dataset_sliced.T[1], vel_cart.T[0], vel_cart.T[1], cmap='seismic')
fig2.savefig('plots/vector_plot_test.pdf', format='pdf', dpi=200)
print('\n:: Everything went smooth, yay!\n')
