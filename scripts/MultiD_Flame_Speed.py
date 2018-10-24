import cantera as ct
import numpy as np
from scipy.integrate import simps

NX = 1300
NY = 192
NZ = 1
delx = 7.5e-6
dely = 7.5e-6
delz = 1
starttime = 1.921e-3
endtime = 3.797e-3
timestepsize = 5e-9
checkpointfrequency = 200
gas = ct.Solution('dme30.xml')

if (starttime==0):
  l0 = int(-1)
else:
  l0 = int((-1)+((endtime-starttime)/(timestepsize*checkpointfrequency)))
file_path = '/media/sq1/Newton/7.5CAD_B8_240us_0.2_W7/data/'

xgrid = np.linspace(0.0, delx*(NX-1), NX)
ygrid = np.linspace(0.0, dely*(NY-1), NY)

MW_O2 = gas.molecular_weights[gas.species_index("O2")]

for i in np.arange(starttime,endtime,timestepsize*checkpointfrequency):
  str_time = ("%.5E"%i)
  P = np.fromfile(file_path+str_time+'/P')
  T = np.fromfile(file_path+str_time+'/T')
  Y = np.fromfile(file_path+str_time+'/Y')
  WDOT = np.fromfile(file_path+str_time+'/WDOT')

  P = P.reshape(NZ,NY,NX)
  P = np.swapaxes(P,0,2)

  T = T.reshape(NZ,NY,NX)
  T = np.swapaxes(T,0,2)

  Y = Y.reshape(gas.n_species,NZ,NY,NX)
  Y = np.swapaxes(Y,1,3)

  WDOT = WDOT.reshape(gas.n_species,NZ,NY,NX)
  WDOT = np.swapaxes(WDOT,1,3)

  rho = np.zeros((NX,NY,NZ))
  
  for j in range (NX):
    for k in range (NY):
      for l in range (NZ):
        gas.TPY = T[j,k,l], P[j,k,l], Y[:,j,k,l]
        rho[j,k,l] = gas.density

  Sc = simps(simps(WDOT[gas.species_index("O2"),:,:,0]*MW_O2,ygrid),xgrid) / ((NY-1.0)*dely * np.average(rho[0,:,:]*(Y[gas.species_index("O2"),-1,:,:]-Y[gas.species_index("O2"),0,:,:])))

  print(Sc)
