"""
Python script to check energy versus strain curves are a good fit.
"""
import numpy as np
import checkfit
import os
import sys
import glob

factor=1.
if (os.path.exists('quantum-espresso') or os.path.exists('vasp')): 
    factor=0.5

if (os.path.exists('planar')): 
    alat     = float(checkfit.pick_from_file('planar',1,1))
    covera   = float(checkfit.pick_from_file('planar',1,2))
    acfactor = alat*covera*5.2917721092e-2
    factor=factor*acfactor

startorder=0
if (os.path.exists('startorder')): 
    input_startorder = open('startorder',"r")
    startorder=int(input_startorder.readline().strip().split()[0])
   
lINFO = os.path.exists('INFO-elastic-constants')
levs = os.path.exists('energy-vs-strain')

lelastic = False    
listene = sorted(glob.glob("*"+"-Energy.dat"))
if (len(listene) > 0): lelastic = True

if (not(lelastic)):
    if (not(lINFO)): sys.exit("ERROR: file INFO-elastic-constants not found!\n")
    if (not(levs)):  sys.exit("ERROR: file energy-vs-strain not found!\n")
    
    input_info = open('INFO-elastic-constants',"r")
    for i in range(3): line = input_info.readline()
    volume = float(input_info.readline().strip().split()[6])
    defcod = int(input_info.readline().strip().split()[3])
    deflab = input_info.readline().strip().split()[3]
    input_info.close()
    inputene = "energy-vs-strain"

if (lelastic):
    input_info = open('../INFO_ElaStic',"r")
    for i in range(4): line = input_info.readline()
    try:
      volume = float(input_info.readline().strip().split()[6])
    except TypeError:
        sys.exit("Could not convert input volume to float.")
    defcod = 99
    deflab = "******"
    inputene = str(listene[0])
    
#-------------------------------------------------------------------------------
order_of_derivative, maximum_strain = checkfit.ask_for_params("strain")
print("\n")
#-------------------------------------------------------------------------------

orderstep = 1
bohr_radius = 0.529177
joule2hartree = 4.3597482
joule2rydberg = joule2hartree/2.
unitconv = joule2hartree/bohr_radius**3/volume*10.**3*factor

#-------------------------------------------------------------------------------

energy: list[float] = []
strain: list[float] = []

#-------------------------------------------------------------------------------

input_energy = open(inputene,"r")

while True:
    line = input_energy.readline()
    line = line.strip()
    if len(line) == 0: break
    eta,ene = line.split() 

    if (abs(float(eta)) <= maximum_strain):  
       energy.append(float(ene))
       strain.append(float(eta))

#-------------------------------------------------------------------------------

strain,energy=checkfit.sort_entries(strain,energy)

#-------------------------------------------------------------------------------

orderlist: list[int]=[]

for i in range(startorder,startorder+6):
    dumorder=order_of_derivative+orderstep*i
    orderlist.append(dumorder)   

#-------------------------------------------------------------------------------

output_file = open('check-energy-derivatives',"w")

n_ini = len(strain)

while (len(strain) > order_of_derivative and len(strain) > 1): 
    bb = []
    n=len(strain)
    etam=max(strain)
    emin=min(strain)
    etam=max(abs(etam),abs(emin))
    for order in orderlist: 
        cc=checkfit(order,strain,energy,n,order_of_derivative)
        if (cc != 99999999):
            bb.append(cc*unitconv)
        else:
            bb.append(cc)
    printderiv(etam,bb,6,output_file)
    if (n_ini == len(strain)):
        printnice(etam,bb,6,order_of_derivative,orderstep,n_ini,defcod,deflab)

    if (abs(strain[0]+etam) < 1.e-7): 
        strain.pop(0)
        energy.pop(0)

    if (abs(strain[len(strain)-1]-etam) < 1.e-7): 
        strain.pop()
        energy.pop()

output_file.close()

output_file = open('order-of-derivative',"w")
print >> output_file, order_of_derivative
output_file.close()
print "###########################################\n"
os.system("export PDEFAULT=true; PLOT-checkderiv.py")
os.system("unset PDEFAULT")





