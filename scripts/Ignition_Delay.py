import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

def init_plotting():
    plt.rcParams['figure.figsize'] = (10, 5)
    plt.rcParams['font.size'] = 10
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['savefig.dpi'] = 600
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
    plt.rcParams['axes.linewidth'] = 4

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

init_plotting()

phi1 = 0.4
phi2 = 0.4
phi_step = 0.1
t1 = 800
t2 = 800
temp_step = 25
i = 500 #Number of iterations
pressure = 20*ct.one_atm
temp_total_steps = int((np.ceil((t2-t1)/temp_step))+1)
phi_total_steps = int((np.ceil((phi2-phi1)/phi_step))+1)
tau1 = np.zeros((temp_total_steps,phi_total_steps))
tau2 = np.zeros((temp_total_steps,phi_total_steps))
phi_idx = -1
temp_idx = -1
temp = np.zeros((temp_total_steps, phi_total_steps))
temp_array = np.zeros((i,temp_total_steps))
total_temp_array = np.zeros((i*phi_total_steps, temp_total_steps))
total_hrr = np.zeros((i*phi_total_steps, temp_total_steps))
gas = ct.Solution('dme39.cti')
io2 = gas.species_index('O2'); #Index of O2 in mix
in2 = gas.species_index('N2'); #Index of N2 in mix
ich3och3 = gas.species_index('CH3OCH3'); #Index of CH3OCH3 in mix
ich4 = gas.species_index('CH4') #Index of CH4 in mix
composition = np.zeros(gas.n_species)
times = np.zeros(i)
tms = np.zeros((i,1))
hrr = np.zeros((i,1))
time_step = 1.e-5

for phi in np.arange(phi1, phi2+0.001, phi_step):
    phi_idx +=1
    for t in np.arange(t1,t2+0.001,temp_step):
        temp_idx += 1
        composition[ich3och3] = phi
        composition[io2] = 3.0
        composition[in2] = 3*3.76
        composition[ich4] = 0.0
        gas.TPX = t, pressure, composition
        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])
        time = 0.0
        IDT1 = 0.0
        IDT2 = 0.0
        temp[temp_idx, phi_idx] = t
        print('%10s %10s %10s %14s %14s'% ('t [s]','T [K]','P [Pa]','u [J/kg]','HRR [W/m3]'))
        for n in range(i):
            time += time_step #time in s
            sim.advance(time)
            times[n] = time*1e3  # time in ms
            temp_array [n, temp_idx] = r.T
            hrr[n] = (-1)*np.dot(gas.net_production_rates, gas.partial_molar_enthalpies)
            print('%10.3e %10.3f %10.3f %14.6e %14.6e'  % (sim.time, r.T, r.thermo.P, r.thermo.u, hrr[n]))
            
        for n in range(1,i+1,1):
            IDT1 += time_step
            if (hrr[n-1] > hrr[n]):
                tau1[temp_idx,phi_idx] = IDT1*1.e3
                print('Tau-1 (ms) = ', tau1[temp_idx,phi_idx])
                break
                
        for n in range(i):
            IDT2 += time_step
            if (hrr[n] == np.amax(hrr)):
                tau2[temp_idx,phi_idx] = IDT2*1.e3
                print('Tau-2 (ms) = ', tau2[temp_idx,phi_idx])
                break   

    temp_idx = -1
    print ('End of phi = ', phi1)
    total_temp_array[phi_idx*i:(phi_idx+1)*len(temp_array),:] = temp_array[:,:]
    total_hrr[phi_idx*i:(phi_idx+1)*len(temp_array),:] = hrr[:]

#print(temp)
#print(tau)
plt.clf()
plt.figure(1)
#plt.subplot(211)
#plt.plot(temp, tau2, marker="H", color='red', lw=4)
#plt.xlabel(r'$\rm\bf{T}_u$' + ' (K)', fontsize=40, fontweight='bold')
#plt.ylabel(r'$\bf\tau$' + ' (ms)', fontsize=40, fontweight='bold')
#plt.xticks(fontsize=40, fontweight='bold')
#plt.yticks(fontsize=40, fontweight='bold')
#plt.plot((750, 750), (3.5, 4), 'k-', lw=4)
#plt.plot((950, 950), (3.5, 4), 'k-', lw=4)
#plt.annotate('',xy=(750, 3.75), xytext=(950, 3.75), arrowprops = dict(arrowstyle="<->", lw=3), fontsize=40)
#plt.text(850, 4, 'NTC Regime', horizontalalignment='center', fontweight='bold', fontsize=40)
#plt.title('Ignition Delay Time Vs. Initial Temperature for DME-Air at P=20 atm, Phi=0.4')
#plt.subplot(212)
#plt.figure(2)
plt.plot(times,total_temp_array[0:i,0], 'blue',  label = 'T=%.f'%t1) #marker = ".",
#plt.plot(times,total_temp_array[0:i,1], 'green' , label = 'T=720 K') #marker = ","
#plt.plot(times,total_temp_array[0:i,2], 'crimson', label = 'T=740 K') #marker = "o", 
#plt.plot(times,total_temp_array[0:i,3], 'yellow', label = 'T=760 K') #marker = "v", 
#plt.plot(times,total_temp_array[0:i,4], 'blue' , label = 'T=780 K') #marker = "^", 
#plt.plot(times,total_temp_array[0:i,5], 'cyan',  label = 'T=800 K') #marker = "*",
#plt.plot(times,total_temp_array[0:i,6], 'black' ,  label = 'T=820 K') #marker = "s",
#plt.plot(times,total_temp_array[0:i,7], 'cadetblue',  label = 'T=840 K') #marker = "D",
#plt.plot(times,total_temp_array[0:i,8], 'brown', label = 'T=860 K') # marker = ">",
#plt.plot(times,total_temp_array[0:i,9], 'darkmagenta',  label = 'T=880 K') #marker = "H",
#plt.plot(times,total_temp_array[0:i,10], 'darkslategray',  label = 'T=900 K') #marker = "<",
#plt.xlabel('Time(ms)')
#plt.ylabel('Temperature(K)')
#plt.title('T Vs. t for %s-Air at T=%.f K, P=%.f atm, Phi=%.1f'%(gas.species_name(ich3och3), t1, pressure/ct.one_atm, phi1))
#plt.legend()
#plt.legend(bbox_to_anchor=(1.01,1.0), loc=2, borderaxespad=0.)
#plt.subplot(212)
#plt.plot(times,total_hrr[0:i,0], 'blue',  label = 'T=%.f'%t1)
#plt.xlabel('Time(ms)')
#plt.ylabel('Heat Release Rate(W/m^3)')
#plt.title('HRR Vs. t for %s-Air at T=%.f, P=%.f, Phi=%.1f'%(gas.species_name(ich3och3), t1, pressure/ct.one_atm, phi1))
#plt.legend()
#plt.legend(bbox_to_anchor=(1.01,1.0), loc=2, borderaxespad=0.)
plt.show()

