###############################################################
#
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#             
###############################################################

#import :
import sys
sys.path.append('/home/sq1/Build/Cantera-2.1.1/lib/python3.6/site-packages/cantera')
from cantera import *
from matplotlib.pylab import *
import numpy as np
import os
#Functions :

#################################################################
# Prepare your run
#################################################################
#Parameter values :
	
	#General
p          =  20.0*one_atm                 # pressure
tin        =   800.0               # unburned gas temperature
phi        =   1.0
if os.path.exists("baseflame000.xml"):
  os.remove("baseflame000.xml")

	#Initial grids, chosen to be 0.02cm long : 
		# - Refined grid at inlet and outlet, 6 points in x-direction :
initial_grid = [0.0, 1e-6, 1.1e-6, 1.e-5, 2.e-4, 0.98e-3, 0.99e-3, 1e-3] # m
		# - Uniform grid, 6 points in x-direction (import numpy):
#initial_grid = 0.02*array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],'d') # m
		# - Uniform grid of 300 points using numpy :
#initial_grid = np.linspace(0,0.02 , 10)

#Set tolerance properties
tol_ss    = [1.0e-5, 1.0e-8]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-5, 1.0e-8]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)
				    
refine_grid = True                  # True to enable refinement, False to
                                    # disable 				   

#Import gas phases with mixture transport model
gas = Solution('dme39.xml')
#################
#Stoechiometry :

fuel_species = 'CH3OCH3'
m=gas.n_species
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

x = zeros(m,'d')
x[ifuel] = phi
x[io2] = 3.0
x[in2] = 3.0*3.76

#################
#Assembling objects :

	#Set gas state to that of the unburned gas
gas.TPX = tin, p, x

	#Create the free laminar premixed flame
f = FreeFlame(gas, initial_grid)
#f.set_fixed_temperature(650)

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

f.inlet.X = x
f.inlet.T = tin

#################################################################
# Program starts here
#################################################################
#First flame:

	#No energy for starters
f.energy_enabled = False

	#Refinement criteria
f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)

	#Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(50, 50)

	#Set time steps whenever Newton convergence fails
f.set_time_step(5.e-06, [10, 20, 80]) #s

	#Calculation
f.solve(loglevel, refine_grid)

#################
#Second flame:

	#Energy equation enabled
f.energy_enabled = True

	#Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)

	#Calculation and save of the results
f.solve(loglevel, refine_grid)

#################
#Third flame and so on ...:

	#Refinement criteria should be changed ...
f.set_refine_criteria(ratio = 5.0, slope = 0.3, curve = 0.3)

f.solve(loglevel, refine_grid)

#################
#Third flame and so on ...:

	#Refinement criteria should be changed ...
f.set_refine_criteria(ratio = 3.0, slope = 0.1, curve = 0.1)

f.solve(loglevel, refine_grid)

#################
f.set_refine_criteria(ratio = 2.0, slope = 0.05, curve = 0.05, prune = 0.01)

f.solve(loglevel, refine_grid)

#Fourth flame and so on ...

f.set_refine_criteria(ratio = 2.0, slope = 0.02, curve = 0.02, prune = 0.01)

f.solve(loglevel, refine_grid)

print ('mixture averaged flamespeed = ',f.u[0])
#################################################################
# Save your results if needed
#################################################################
#Write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('dme_adiabatic.csv', quiet=False)
f.save('baseflame000.xml','energy_mix')
#f.write_avbp('Sol-CAN2AV_P-'+str(p)+'-T-'+str(tin)+'-Phi-'+str(phi)+'.csv', quiet=False)

#################################################################
# Plot your results
#################################################################
#Plot the velocity, temperature, density
z = f.flame.grid
T = f.T
u = f.u

fig=figure(1)

	# create first subplot - adiabatic flame temperature
a=fig.add_subplot(221)
a.plot(z,T)
title(r'$T_{adiabatic}$ vs. Position')
xlabel(r'Position [m]', fontsize=15)
ylabel("Adiabatic Flame Temperature [K]")
a.xaxis.set_major_locator(MaxNLocator(10)) # this controls the number of tick marks on the axis

	# create second subplot - velocity
b=fig.add_subplot(222)
b.plot(z,u)
title(r'Velocity vs. Position')
xlabel(r'Position [m]', fontsize=15)
ylabel("velocity [m/s]")
b.xaxis.set_major_locator(MaxNLocator(10)) 

        # create third subplot - rho
c=fig.add_subplot(223)
p = zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    p[n]= gas.density_mass
c.plot(z,p)
title(r'Rho vs. Position')
xlabel(r'Position [m]', fontsize=15)
ylabel("Rho [kg/m^3]")
c.xaxis.set_major_locator(MaxNLocator(10)) 


        # create fourth subplot - specie CH4
d=fig.add_subplot(224)
ch4 = zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    ch4[n]= gas.Y[ifuel]
d.plot(z,ch4)
title(r'CH4 vs. Position')
xlabel(r'Position [m]', fontsize=15)
ylabel("CH4 Mole Fraction")
d.xaxis.set_major_locator(MaxNLocator(10))

        # Set title
fig.text(0.5,0.95,r'Adiabatic $CH_{4}$ + Air Free Flame at Phi = 1 Ti = 300K and P = 1atm',fontsize=22,horizontalalignment='center')

subplots_adjust(left=0.08, right=0.96, wspace=0.25)

show()

f.show_stats

