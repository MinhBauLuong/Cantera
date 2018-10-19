#Script for transfering the cantera solution based on adaptive mesh refinement to a solution with regularly spaced grid 

import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import math
import sys
import cantera as ct
import fileinput
import csv

dns = ct.Solution('dme30.xml')
dnsnames = dns.species_names
runid = '115'
zstr = '0.00116'
f = ct.FreeFlame(dns) #defining the gas object
f.restore('baseflame000.xml','energy_mix') #restoring the previously calculated cantera solution
#f.restore('flame'+str(runid)+'-'+str(zstr)+'.xml','energy_mix') #restoring the previously calculated cantera solution
grid_size = f.grid.size #declaring the variable for storing the grid size of the domain containing the flame
MF = f.Y #mass fractions of all 9 species
NPR = f.net_production_rates #net production rates of all 9 species
S = list() #chemical symbols of all 30 species
for i in range (len(dnsnames)):
  S.insert(i, dnsnames[i])
TS = len(dnsnames) #total number of species
Z = f.grid #position of each grid point in the flame in m
V = f.u #flame speed at each grid point in m/s
T = f.T #flame temperature at each grid point in K
D = f.density_mass #flame density at each grid point in kg/m^3
P = f.P #flame pressure at each grid point in Pa
HRR = (-1.0)*np.dot(f.net_production_rates.transpose(), f.partial_molar_enthalpies)
TAS = 2*TS+6 #number of columns in the old/new array
Names = ['z (m)','u (m/s)','T (K)','rho (kg/m3)', 'P (Pa)'] #list containing the titles excluding species symbols to be printed in the first row of the output csv file
FR = Names + S + S + ['HRR (W/m^2)'] #list containing all the titles including species symbols to be printed in the first row of the output csv file
oldarray = np.zeros((grid_size,TAS)) #array containing the flame solution based on cantera's adaptive grid
max_production_rate = np.zeros((TS,1)) #array containing the maximum net production rate for each species
min_production_rate = np.zeros((TS,1)) #array containing the minimum net production rate for each species
z1 = np.zeros((TS,1)) #array containing the position of the maximum net production rate for each species    
z2 = np.zeros((TS,1)) #array containing the position of the minimum net production rate for each species 
P2P = np.zeros((TS,1)) #array containing the distance between the positions of maximum and the minimum net production rate for each species in m

for n in range(grid_size): #control loop for storing the flame solution in the 'oldarray' defined above
  oldarray[n,0] = Z[n] #z(m)
  oldarray[n,1] = V[n] #u (m/s)
  oldarray[n,2] = T[n] #T (K)
  oldarray[n,3] = D[n] #rho (kg/m3)
  oldarray[n,4] = P #P (Pa)
  for m in range (TS):
    oldarray[n,m+5] = MF[m,n] #Mass fractions of all 30 species
    oldarray[n,m+35] = NPR[m,n] #Net production rates of all 30 species
    oldarray[n,65] = HRR[n][0]
            
for i in range (TS): #control loop for finding the maximum and minimum net production rates, their difference, their positions on
                     #the grid as well as the distance between their respective positions
    s = dns.species_name(i)
    print('\n')
    print('Species =',s)
    max_production_rate[i] = np.amax(NPR[i,:])
    min_production_rate[i] = np.amin(NPR[i,:])
    for n in range (grid_size):
      if (NPR[i,n] == max_production_rate[i]):
        z1[i] = f.grid[n]
      if (NPR[i,n] == min_production_rate[i]):
        z2[i] = f.grid[n] 
    print('Max. Net Prod. Rate =',float(np.round(max_production_rate[i],12)))
    print('Min. Net Prod. Rate =',float(np.round(min_production_rate[i],12)))
    print('Z at Max-NPR =',float(np.round(z1[i],8)))
    print('Z at Min-NPR =',float(np.round(z2[i],8)))
    P2P[i] = abs(z2[i] - z1[i])
    print ('Distance from Peak to Peak =',float(np.round(P2P[i],8)))

r = np.zeros(TS-1) #temperory array for storing the distance between the maximum and the minimum net production rate for each species excluding N2 in m

for j in range(TS):
    if (P2P[j]>0):
        r[j] = P2P[j]

min_P2P = np.amin(r) #least distance between the maximum and the minimum net production rate amongst all the species excluding N2 in m 

for k in range(TS): #control loop for finding the species (excluding N2) with the least distance between the maximum and the minimum net production rate  
    if (P2P[k] == min_P2P):
        l = dns.species_name(k)
        z_max = z2[k]
        z_min = z1[k]

print('\n')
print('Species with Min. P2P Distance =',l)
print('Min. P2P Distance = %.5e'%float(np.round(min_P2P,8)))
print('Z at Max-NPR =',float(np.round(z_min,8)))
print('Z at Min-NPR =',float(np.round(z_max,8)))
print('Cantera Spacing =',float(np.round(z_max-z_min,8)))
print('Min. Spacing Required =',float(np.round(min_P2P/12,8)))

user_input = input('\nEnter the minimum spacing for the new grid (should be less than the minimum spacing required):') #user input for specifying the spacing for the new equi-spaced grid
new_grid_spacing=float(user_input) #reading the user input for the new regular grid spacing
allowance = new_grid_spacing*1e-2 #allowance that will be added to the end of the control structure generating the new grid so that all the points in the old cantera grid are 
                                  #accomodated in the new grid

newgrid = np.arange(Z[0], Z[grid_size-1]+allowance, new_grid_spacing) #defining the new grid with the user-defined spacing
m = np.size(newgrid) #size of the new grid
newarray = np.zeros((m,TAS)) #array containing the flame solution based on newly defined equi-spaced grid
newarray[:,0] = newgrid #position of each grid point in the new grid
c = 5 #index for counting the number of species in the new grid
d = 35 #index for counting the number of species in the new grid
for i in range(1,TAS,1): #control loop for interpolating the data contained in the old array based on cantera's adaptive grid and storing the interpolated data in the new array
                         #based on equi-spaced grid 
    g = interpolate.interp1d(oldarray[:,0],oldarray[:,i], kind='linear', fill_value='extrapolate', bounds_error=False)
    if (i == 1):
        newarray[:,1] = g(newgrid) #Interpolated Velocity
    if (i==2):
        newarray[:,2] = g(newgrid) #Interpolated Temperature
    if (i==3):
        newarray[:,3] = g(newgrid) #Interpolated Density
    if (i==4):
        newarray[:,4] = g(newgrid) #Interpolated Pressure
    if (i>4 and i<35):
        newarray[:,c] = g(newgrid) #Interpolated Mass Fractions
        c+=1
    if (i>34 and i<65):
        newarray[:,d] = g(newgrid) #Interpolated Net Production Rates
        d+=1
    if (i>64):
        newarray[:,i] = g(newgrid) #Interpolated heat release rates
        

csv_file = ('Mapped_Solution.csv') #output excel file containing the flame solution based on new equi-spaced grid
with open(csv_file, 'w') as outfile: #control loop for writing information in the above output file
    writer = csv.writer(outfile) #converting the data into delimited strings
    writer.writerow(FR) #writing the first row in the output file
    for i in range(m): #control loop for writing the data for each species at each grid point position and temperature
        writer.writerow([newarray[i,0], newarray[i,1], newarray[i,2], newarray[i,3], newarray[i,4]] + list(newarray[i,5:65]) + [newarray[i,65]]) #writing the data row-wise
    w = newarray[-1,0]
    while(m!= int(math.ceil(m/100.0)*100)):
        writer.writerow([w+new_grid_spacing, newarray[-1,1], newarray[-1,2], newarray[-1,3], newarray[-1,4]] + list(newarray[-1,5:65]) + [newarray[-1,65]]) #writing the data row-wise
        w = w+new_grid_spacing
        m = m+1        
print('\nOutput written to {0}'.format(csv_file))

print('\nNote:')
print('Old grid had',grid_size, 'points') #printing the number of grid points contained in cantera's adaptive grid
print('New grid has',m, 'points') #printing the number of grid points contained in the newly defined equi-spaced grid

plt.clf() #plotting the temperature profile as well as the net production rate profile of species (excluding N2) with the least distance between the
          #maximum and the minimum net production rate based on cantera's adaptive grid as well as the new equi-spaced grid for comparison
plt.figure(1)
plt.subplot(211)
plt.plot(oldarray[:,0], oldarray[:,2],'blue', marker='p') 
plt.title('Temperature profile - Cantera grid')
plt.xlabel('Z (m)')
plt.ylabel('T (K)')
#plt.xlim(0, 0.04)
plt.subplot(212)
plt.plot(newarray[:,0], newarray[:,2],'crimson', marker='H')
plt.title('Temperature profile - New regular grid')
plt.xlabel('Z (m)')
plt.ylabel('T (K)')
#plt.xlim(0, 0.04)
plt.figure(2)
plt.subplot(211)
plt.plot(oldarray[:,0], oldarray[:,dns.species_index(l)+len(Names)],'blue', marker='p')
plt.title('Mass fraction/Net production rate of the species with minimum peak-to-peak distance - Cantera grid')
plt.xlabel('Z (m)')
plt.ylabel('Net production rate (kmol/m^3-s)')
#plt.xlim(0, 0.04)
plt.subplot(212)
plt.plot(newarray[:,0],newarray[:,dns.species_index(l)+len(Names)],'crimson', marker='H')
plt.title('Mass fraction/Net production rate of the species with minimum peak-to-peak distance - New regular grid')
plt.xlabel('Z (m)')
plt.ylabel('Mass fraction/Net production rate (kmol/m^3-s)')
#plt.xlim(0, 0.04)
plt.show()
