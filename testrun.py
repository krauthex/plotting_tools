import main
import numpy as np

# data
print('\nLoading Data...\n')
vr, vtheta, vphi = np.loadtxt('data/prepared/VEL_PL', skiprows=1, unpack=True)  # velocities
r, theta, phi = np.loadtxt('data/prepared/SPH_COOR_PL.txt', skiprows=1, unpack=True)  # spherical coordinates
rho, temp = np.loadtxt('data/raw/RHO_TEMP.txt', skiprows=1, unpack=True)  # density and temperature

# computation
print('\nComputation/Data Extraction...\n')
velocities, indices = main.extract_shell(vr, vtheta, vphi, r)  # shell extractor
x, y = main.hammer_projection(theta[indices], phi[indices])  # hammer projection
print('Done.\nPlotting...')
"""
fig, ax, cax = main.contour_plot_velocity(x, y, velocities[0], figsize=(12,9))  # contour plot

# plotting Section
ax.set_xlabel('$\\phi\ \mathrm{[rad]}$')
ax.set_ylabel('$\\theta\ \mathrm{[rad]}$')
ax.set_xlim(-3.5, 4)
ax.set_ylim(-1.5, 1.5)

fig.savefig('test.pdf', dpi=200, format='pdf')
print('Saved the first plot.\nSecond plot coming...\n')

fig2, ax2, cax2 = main.scatter_contour(r, np.log(rho), num_levels=3)
# plottin section 2
ax2.set_xlabel('$\mathrm{log}\ \\frac{R}{R_{Planet}}$')
ax2.set_ylabel('$Density\ \mathrm{[g/cm^3]}$')

# cosmetics
ax2.set_xlim(0.1, 1e3)

fig2.savefig('test3.png', format='png')
"""
