#!/usr/bin/python3

import sys
sys.path.append('/home/swapnil/Cantera/lib/python3.4/site-packages')
import cantera as ct
import numpy as np
import fileinput
import matplotlib.pyplot as plt

iterations = 490
F_S = np.zeros(iterations)
L = np.zeros(iterations)
idx = 0

def solveflame(savefilename):

  loglevel = 1  # amount of diagnostic output (0 to 8)
  refine_grid = True  # 'True' to enable refinement, 'False' to disable

  # IdealGasMix object used to compute mixture properties
  gas = ct.Solution('dme39.xml')
  
  # Flame object
  f = ct.FreeFlame(gas)

  #----------------------------------------------------------------------

  f.restore('restart.xml','energy_mix')
  f.transport_model = 'Mix'
  f.set_max_jac_age(10, 10)
  f.set_time_step(1e-5, [2, 5, 10, 20, 50])
  f.solve(loglevel, refine_grid)
  f.save(savefilename,'energy_mix', 'solution for a given z left boundary')
  #f.write_csv(csvfilename)
  #print('z = %s S-L = %.15e' %(zstr, f.u[0]))
  F_S[idx] = f.u[0]
  L[idx] = 0.005-float(zstr) 
  print('L = %.5f, S-L = %.5f' %(L[idx], F_S[idx]))  

  
#----------------------------------------------------------------------
# MAIN MAIN
#----------------------------------------------------------------------
runid='514'
zstr='0.00515'
#newfilename='baseflame'+runid+'.xml'
newfilename='flame'+runid+'-'+zstr+'.xml'
delz=0.00001

for n in range(iterations):
  runid="%03d"%(n+515)
  #----------------------------------------------------------------------
  oldfilename=newfilename
  of=open(oldfilename,'r')
  nf=open("restart.xml",'w')
  readnext=0 
  for oldline in of:
    if readnext==1 :
      words=oldline.split()
      newz=float(words[0][:-1])-delz
      zstr="%.15E" %(newz)
      words[0]=zstr+','
      newline="  ".join(words)    
      newline='          '+newline+'\n'
      zstr="%.5f"%(newz)
      readnext=0
    else:
      newline=oldline
    if 'title="z"' in oldline:
      readnext=1
    nf.write(newline)
  of.close()
  nf.close()
  newfilename='flame'+runid+zstr+".xml"
  csvfilename='soln'+runid+'_'+zstr+".csv"
  
  #----------------------------------------------------------------------
  solveflame(newfilename)
  idx +=1


plt.clf()
plt.plot(L,F_S)
plt.xlabel('Induction Length (m)')
plt.ylabel('Flame Speed (m/s)')
plt.title('Transition from Deflagration to Spontaneous propagation for DME-Air at T=800 K, P=40 atm, Phi=0.4')
plt.show()
