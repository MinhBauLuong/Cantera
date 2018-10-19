import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import math
import sys
import cantera as ct
import fileinput
import csv
import string
import matplotlib.animation as animation
from scipy.signal import savgol_filter
from sympy import *

def init_plotting():
    plt.rcParams['figure.figsize'] = (10, 5)
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.minor.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

init_plotting()

dns = ct.Solution('dme30.xml')
species = 'O2'
species_index = dns.species_index(species)
MW_species = dns.molecular_weights[species_index]
rows1= [4200] #3300, 2300, 5600, 4200 
rows = rows1[0]
Yk_diff_coeff_mass = np.zeros((rows))
Yk_diff_coeff_therm = np.zeros((rows))
Vk_mass = np.zeros((rows))
Vk_therm = np.zeros((rows))
Vk = np.zeros((rows))
D = np.zeros((rows))
dD_dz = np.zeros((rows))
Z = np.zeros((rows))
P = np.zeros((rows))
T = np.zeros((rows))
U = np.zeros((rows))
Y = np.zeros((rows))
rho = np.zeros((rows))
WDOT = np.zeros((rows))
Y_row = np.zeros((rows, dns.n_species))
Convection = np.zeros(rows)
Diffusion = np.zeros(rows)
Reaction = np.zeros(rows)
time = [1,1] # O2
for l in range(len(time)-1):
        str_time = ("%.5E"%time[l])
        filename = "1060_Extended_4200.csv"
        print('\nReading from file:%s'%filename)
        with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter=',')
                f.readline()
                j = -1
                for row in reader:
                        if (not row or j>rows):
                                break;
                        else:
                                j += 1
                                Z[j] = float(row[0])*1e3
                                U[j] = float(row[1])
                                T[j] = float(row[2])
                                rho[j] = float(row[3])
                                P[j] = float(row[4])
                                Y[j] = float(row[5+species_index])
                                WDOT[j] = float(row[35+species_index])*MW_species
                                for w in range (dns.n_species):
                                        Y_row[j,w] = float(row[5+w])


        delx = (Z[1] - Z[0])/1e3
        dYdx = np.gradient(Y,delx, edge_order=2)
        dTdx = np.gradient(T,delx, edge_order=2)
        for i in range (rows):
            dns.TPY = T[i],P[i],Y_row[i,:]
            Yk_diff_coeff_mass[i] = dns.mix_diff_coeffs_mass[species_index]
            Yk_diff_coeff_therm[i] = dns.thermal_conductivity/(dns.density_mass*dns.cp_mass)
            Vk_mass[i] = -1.0/Y[i]*Yk_diff_coeff_mass[i]*dYdx[i]
            Vk_therm[i] = (-1.0/(rho[i]*Y[i]*T[i]))*Yk_diff_coeff_therm[i]*dTdx[i]
            Vk[i] = Vk_mass[i] + Vk_therm[i]

        for i in range (rows):
            Convection[i] = rho[i]*Y[i]*U[i]
            Reaction[i] = WDOT[i]
            Diffusion[i] = -rho[i]*Y[i]*Vk[i]
            

        Convection = np.gradient(Convection, delx, edge_order=2)
        Convection = savgol_filter(Convection, 51, 11)
        Diffusion = np.gradient(Diffusion, delx, edge_order=2)
        #Diffusion = Convection - Reaction
        #Diffusion = savgol_filter(Diffusion, 51,11)

csv_file = ('Transport_Budget_1060_Extended.csv') #output excel file containing the flame solution based on cantera's adaptive grid
with open(csv_file, 'w') as outfile: #control loop for writing information in the above output file
    writer = csv.writer(outfile) #converting the data into delimited strings
    writer.writerow(['Z (mm)', 'Convection', 'Diffusion', 'Reaction'])
    for i in range(rows):
        writer.writerow([Z[i]*1e3, Convection[i], Diffusion[i], Reaction[i]])
print('Output written to {0}\n'.format(csv_file))
        
plt.clf()
plt.subplot(221)
##plt.plot(Z,-Convection,color='green',lw=3, label='Convection')
plt.plot(Z,Diffusion,color='blue',lw=3, label='Diffusion')
plt.twinx()
plt.plot(Z,Reaction,color='red',lw=3, label='Reaction')
plt.legend(loc='upper right')
plt.show()
