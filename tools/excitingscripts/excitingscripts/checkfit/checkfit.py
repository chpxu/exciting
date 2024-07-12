"""
Central file containing functions used by CHECKFIT files.
"""
from typing import IO, Any, List, Literal, Union, Tuple
import numpy as np
import os
import sys

def fit(order: int, x: list[float | int],y: list[float | int],npoints: int,nderiv: int) -> int:
    import math
    b: int = 99999999
    if (npoints > order): 
       b = np.polyfit(x,y,order)[order-nderiv]
       b = math.factorial(nderiv)*b
    return b

def printderiv(etam: float,bb: list[float],nfit: int,output_file: IO[Any]) -> None:
    print('%12.8f'%(etam), file=output_file)
    for j in range(0,nfit):
        if (bb[j] != 99999999):
            print('%14.6f'%(bb[j]), file=output_file),
    print("", file=output_file)
    return 

def printnice(etam: Union[float, int],bb: List[float],nfit:int,nderiv: int,orderstep: int, displpoints: int) -> None:
    invcm2hz = 33.356409
    print("\n")
    print(f"#############################################\n")
    print(f"Fit data------------------------------------\n")
    print(f"Maximum value of the displacement       ==>    {'%5.3f'%(etam)}")
    print(f"Number of displacement values used      ==>    {'%5i'%(displpoints)}\n")
    print(f"Fit results for the derivative of order ==>    {'%4i'%(nderiv)}\n")
    for j in range(0,nfit):
        order= nderiv+j*orderstep
        if (bb[j] != 99999999):
            print(f"Polynomial of order {'%2i'%(order)} ==> {'%8.2f'%(bb[j])} [cm^-1])")
    print("\n")
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print(f"Polynomial of order {'%2i'%(order)} ==> {'%8.4f'%(bb[j]/invcm2hz)} [THz]")
    print("\n")
    return 

def sort_entries(s: list[float],e: list[float]) -> Tuple[list[float], list[float]]:
    """
    Sort quantities.
    """
    ss: list[float]=[]
    ee: list[float]=[]
    ww: list[float]=[]
    for i in range(len(s)): ww.append(s[i])
    ww.sort()
    for i in range(len(s)):
        ss.append(s[s.index(ww[i])])
        ee.append(e[s.index(ww[i])])
    return ss, ee


def print_strain_nice(etam: Union[float, int],bb: List[float],nfit:int,nderiv: int,orderstep: int, strainpoints: int, dc: str,dl: str) -> None:
    """
    Pretty print information from the strain file.
    """
    punit = "[GPa]"
    if (os.path.exists('planar')): punit = "[N/m]"
    print("\n")
    print ("###########################################\n")
    print ("Fit data-----------------------------------\n")
    print (f"Deformation code             ==> {dc}\n")#'%2i'%(dc))
    print (f"Deformation label            ==> {dl}\n")
    print (f"Maximum value of the strain  ==> {'%10.8f'%(etam)}\n") 
    print (f"Number of strain values used ==> {'%2i'%(strainpoints)}\n") 
    print (f"Fit results for the derivative of order {'%3i'%(nderiv)}\n")  
    for j in range(0,nfit):
        order=nderiv+j*orderstep
        if (bb[j] != 99999999):
            print(f"Polynomial of order: {'%2i'%(order)} ==> {'%10.2f'%(bb[j])}, {punit}")
    print("\n")
    return 

def pick_from_file(filename: str,nrow: int,ncol: int) -> str:
    """
    Output the entry at position (nrow, ncol)
    """
    pfile = open(filename,"r")
    line = ""
    for _ in range(nrow): line = pfile.readline().strip()
    x = line.split()[ncol-1]
    pfile.close()
    return x 

def ask_for_params(mode: Literal['displ', 'strain']) -> Tuple[int, float]:
    order_of_derivative = input("\nEnter the order of derivative >>>> ")
    try: 
        order_of_derivative = int(order_of_derivative)
    except ValueError:
        sys.exit("Order of derivative is not integer.")
    if (order_of_derivative < 0): 
        sys.exit("ERROR: Order of derivative must be positive!\n")
    if (mode == 'displ'):
        maximum_displ = input("\nEnter maximum displacement for the fit >>>> ")
        try:
            maximum_displ = float(maximum_displ)
            return order_of_derivative, maximum_displ
        except ValueError:
            print("Displacement is not float.")
            sys.exit()
    if (mode == 'strain'):
        maximum_strain = input("\nEnter maximum strain for the fit >>>> ")
        try:
            maximum_strain = float(maximum_strain)
            return order_of_derivative, maximum_strain
        except ValueError:
            print("Maximum strain is not float.")
            sys.exit()


def order_and_print(mode: Literal['strain', 'displ'], 
                    quantity: list[float], 
                    energy: list[float], 
                    target: list[Any], 
                    order_of_derivative: int, 
                    orderstep: int, 
                    unitconv: float, 
                    output_file: IO[Any], 
                    defcon: int | None = None, 
                    deflab: str | None = None):
    
    nmax=len(quantity)

    while (len(quantity) > order_of_derivative and len(quantity) > 1): 
        bb: list[Any] = []
        n=len(quantity)
        etam=max(quantity)
        emin=min(quantity)
        etam=max(abs(etam),abs(emin))
        freq: float = 0.0
        for order in target: 
            cc=fit(order,quantity,energy,n,order_of_derivative)
            if (cc != 99999999):
                if (cc >= 0): freq=np.sqrt( cc*unitconv)/2./pi
                if (cc <  0): freq=-np.sqrt(-cc*unitconv)/2./pi
                bb.append(freq)
            else:
                bb.append(cc)
        if (n == nmax and mode == 'displ'):
            printnice(etam,bb,6,order_of_derivative,orderstep,nmax)
        if (n == nmax and mode == 'strain'):
            printnice(etam,bb,6,order_of_derivative,orderstep,n_ini,defcod,deflab)
        printderiv(etam,bb,6,output_file)

        if (abs(quantity[0]+etam) < 1.e-7): 
            quantity.pop(0)
            energy.pop(0)

        if (abs(quantity[len(quantity)-1]-etam) < 1.e-7): 
            quantity.pop()
            energy.pop()