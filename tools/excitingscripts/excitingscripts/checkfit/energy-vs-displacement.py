import numpy as np
from typing import Any, List
import sys
import os
import checkfit

if (os.path.exists('INFO-diamond-phonon') == False ): 
    sys.exit("ERROR: file INFO-diamond-phonon not found!\n")
if (os.path.exists('energy-vs-displacement')==False): 
    sys.exit("ERROR: file energy-vs-displacement not found!\n")

input_info = open('INFO-diamond-phonon',"r")
line = input_info.readline()
lattice_parameter = float(input_info.readline().strip().split()[5])
input_info.close()

#-------------------------------------------------------------------------------
maximum_displ = input("\nEnter maximum displacement for the fit >>>> ")
order_of_derivative = input("\nEnter the order of derivative >>>> ")

try:
    maximum_displ = float(maximum_displ)
    order_of_derivative= int(order_of_derivative)
    if (order_of_derivative < 0): 
      sys.exit("ERROR: Order of derivative must be positive!\n")
except ValueError:
    print("Displacement is not float or order is not int.")
    sys.exit()

amass = input("\nEnter atomic mass in [amu] >>>> ")
try:
    amass = float(amass)
except ValueError:
    print("Atomic mass (amass) is not float.")
    sys.exit()
#-------------------------------------------------------------------------------

orderstep = 1
factor = 2./3.
if (os.path.exists('X-phonon-calculation')):
    factor = 1./4. 
    print ("_____________________________________________")
    print (" X-phonon-calculation")
if (os.path.exists('quantum-espresso')): 
    factor=factor*0.5
amu = 1.66053886
invcm2hz = 33.356409
bohr_radius = 0.529177
pi = 3.1415926535897931
joule2hartree = 4.3597482
invcm2joule = 5.03411721428
joule2rydberg = joule2hartree/2.

unitconv=invcm2hz/lattice_parameter/bohr_radius
unitconv=unitconv**2
unitconv=unitconv*factor*joule2hartree/amass/amu*10**5

#-------------------------------------------------------------------------------

energy: List[float] = []
displ: List[float] = []

#-------------------------------------------------------------------------------

input_energy = open('energy-vs-displacement',"r")

while True:
    line = input_energy.readline()
    line = line.strip()
    if len(line) == 0: break
    eta,ene = line.split() 

    if (abs(float(eta)) <= maximum_displ):  
       energy.append(float(ene))
       displ.append(float(eta))

input_energy.close()
#-------------------------------------------------------------------------------

displ, energy = checkfit.sort_entries(displ,energy)

#-------------------------------------------------------------------------------

orderlist: list[int] = []

for i in range(0,6):
    dumorder=order_of_derivative+orderstep*i
    orderlist.append(dumorder)  

#-------------------------------------------------------------------------------

output_file = open('check-energy-derivatives',"w")
nmax=len(displ)

while (len(displ) > order_of_derivative and len(displ) > 1): 
    bb: list[Any] = []
    n=len(displ)
    etam=max(displ)
    emin=min(displ)
    etam=max(abs(etam),abs(emin))
    freq = 0
    for order in orderlist: 
        cc=checkfit.fit(order,displ,energy,n,order_of_derivative)
        if (cc != 99999999):
            if (cc >= 0): freq=np.sqrt( cc*unitconv)/2./pi
            if (cc <  0): freq=-np.sqrt(-cc*unitconv)/2./pi
            bb.append(freq)
        else:
            bb.append(cc)
    if (n == nmax):
        checkfit.printnice(etam,bb,6,order_of_derivative,orderstep,nmax)
    checkfit.printderiv(etam,bb,6,output_file)

    if (abs(displ[0]+etam) < 1.e-7): 
        displ.pop(0)
        energy.pop(0)

    if (abs(displ[len(displ)-1]-etam) < 1.e-7): 
        displ.pop()
        energy.pop()

output_file.close()

output_file = open('order-of-derivative',"w")
print(order_of_derivative, file=output_file)
output_file.close()
print("#############################################\n")
os.system("export PDEFAULT=true; PLOT-checkderiv.py")
os.system("unset PDEFAULT")