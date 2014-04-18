'''
Make the ascii input file for particles for the disksplode simulations
Sam Geen, April 2014
'''

import os
import numpy as np

# HARD-CODED!!!
units_density=2.3247434e-24 # mu * 1.66e-24, where mu = 1.4
units_time = 2.5395079e15 # 1/sqrt(G * units_density)
units_length = 3.08e18 # 1 pc
Myrins = 3.15569e13
msolaring = 1.9891e33

def checkexists(fname,num=0):
    '''
    Back up the file it if exists
    '''
    if num != 0:
        full = fname+".saf"+str(num).zfill(5)
    else:
        full = fname
    if os.path.exists(full):
        checkexists(fname,num+1)
    else:
        os.system("cp "+fname+" "+full)

def calc_sfr(sigma_g):
    '''
    Star formation rate in Msolar Myr-1
    Eqn 13 in Creasey et al 2013
    '''
    area = 1.0 # in kpc^2
    # NOTE - eqn 13 is in years, not Myr
    return 2.5e2*(sigma_g**1.4)*area

def calc_nsn(sigma_g,simtime):
    '''
    Calculate the number of supernovae over the age of the simulation
    sigma_g in Msolar/pc^2, simtime in Myr
    '''
    # Number of supernovae per 100Msolar
    # Below equation 14 Creasey et al 2013
    snperm = 1.8e-2
    mpertime = calc_sfr(sigma_g)
    mtot = mpertime * simtime
    nsn = snperm * mtot
    return int(nsn+0.5)

def xdist(fg,sigma_g,nsn):
    '''
    Get the column distribution of stars
    '''
    b = 61 * (fg/0.1) / (sigma_g / 10) # in pc
    xin = np.random.rand(nsn)
    # NOTE - xin is in kpc
    rho = sigma_g / 2.0 / b * np.cosh(xin/(b/1e3))**(-2.0) # Volume density
    mcol = sigma_g
    x = rho / mcol # Distribution is the volume density / column density
    import matplotlib.pyplot as plt
    plt.plot(xin,x,".")
    plt.yscale("log")
    plt.savefig("xdisttest.png")
    return x

def run(fg,sigma_g,simtime):
    '''
    Make the input file
    sigma_g is in Msolar/pc^2
    simtime in Myr
    '''
    # Open the file and back up old files if necessary
    filename = "sn_parts.dat"
    checkexists(filename)
    f = open(filename,"w")
    # Masses stuff
    eff_sn = 1.8 # per 100Msolar of stars
    # Find mass of stars per supernova in Msolar
    mpersn = 1.0/eff_sn * 100.0 * msolaring
    mpersn /= units_density*units_length**3.0 # in code units
    # Time stuff
    toff = 10.0 # 10Myr explosion offset to remove beforehand
    
    # Set up arrays
    nsn = calc_nsn(sigma_g,simtime)
    x = xdist(fg,sigma_g,nsn)
    y = np.random.rand(nsn)
    z = np.random.rand(nsn)
    u = np.zeros(nsn)
    v = np.zeros(nsn)
    w = np.zeros(nsn)
    i = -np.arange(1,nsn+1,dtype=np.int32)
    m = np.zeros(nsn)+mpersn
    # Set supernovae off at random times inside simtime
    # Subtract offset to explode at zero-time
    t = np.random.rand(nsn)*simtime - toff
    t *= Myrins
    t /= -units_time # Has to be -ve for the GMC feedback model
    for j in range(0,nsn):
        f.write(str(x[j])+" "+str(y[j])+" "+str(z[j])+\
                str(u[j])+" "+str(v[j])+" "+str(w[j])+\
                str(m[j])+" "+str(t[j])+"\n")
    # Done!
    f.close()


if __name__=="__main__":
    run(fg=0.1,sigma_g=11.61,simtime=10.0)
