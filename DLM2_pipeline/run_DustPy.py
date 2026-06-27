"""
This package contains function to configure and run DustPy

Developped by A. Blazère in June 2026
"""

'---Intialisation---'
import numpy as np
from dustpy import constants as c
from dustpy import plot

'---Colored prints---'
#Functions from https://www.geeksforgeeks.org/python/print-colors-python-terminal/
def prRed(s): print("\033[91m {}\033[00m".format(s))
def prGreen(s): print("\033[92m {}\033[00m".format(s))
def prYellow(s): print("\033[93m {}\033[00m".format(s))
def prLightPurple(s): print("\033[94m {}\033[00m".format(s))
def prPurple(s): print("\033[95m {}\033[00m".format(s))
def prCyan(s): print("\033[96m {}\033[00m".format(s))
def prLightGray(s): print("\033[97m {}\033[00m".format(s))

# '---Configuring---'
# def configure_DustPy(parameter_array):
#     sim = Simulation() #Creates a simulation object that has DustPy default parameters
#     print("---Configuring DustPy---")
#     if parameter_array[24]==False:
#         print("#Information: User defined parameters will be used#")
#         print("Default parameters were:")
#         sim.ini
        
#         #Dust
#         sim.ini.dust.aIniMax=parameter_array[0] #Maximum particle size that will be filled initially (as in what is the max size that exits initially)
#         sim.ini.dust.allowDriftingParticles=parameter_array[1] #Radial drift on/off
#         sim.ini.dust.erosionMassRatio=parameter_array[2] #This parameter defines the threshold between full fragmentation and erosion
#         sim.ini.dust.d2gRatio=parameter_array[3] #Dust to gas ratio
#         sim.ini.dust.distExp=parameter_array[4] #nital size distribution of particles
#         sim.ini.dust.excavatedMass=parameter_array[5] #If errosive collision, how much mass is chipped of larger collision partner
#         sim.ini.dust.fragmentDistribution=parameter_array[6] #Whenever a particle collision produces fragments the fragments will have a mass distribution following a power law of index...
#         sim.ini.dust.rhoMonomer=parameter_array[7] #Monomer bulk mass density, or dust grain individual density
#         sim.ini.dust.vFrag=parameter_array[8] #Fragmentation velocity, if above erosion or fragmentation, if below sticking try increasing by an order of magnitude default is 1m/s (written 100 cm/s)
    
#         #Grid
#         sim.ini.grid.Nmbpd=parameter_array[9] #Number of mass bins per mass decade
#         sim.ini.grid.mmin=1e5 * np.pi * ((parameter_array[10]*1e-4)**3) * sim.ini.dust.rhoMonomer #Minimum particle mass
#         sim.ini.grid.mmax=1e5 * np.pi * ((parameter_array[11]*1e-4)**3) * sim.ini.dust.rhoMonomer #Maximum particle mass
#         sim.ini.grid.Nr=parameter_array[12] #Number of radial grid cells
#         sim.ini.grid.rmin=parameter_array[13]*14959787070000.0 #Location of inner radial grid boundary
#         sim.ini.grid.rmax=parameter_array[14]*14959787070000.0 #Location of outer radial grid boundary
    
#         #Star
#         sim.ini.star.M=parameter_array[15]*1.989e+33 #Star mass
#         sim.ini.star.R=parameter_array[16]*6.957e+10 #Star radius
#         sim.ini.star.T=parameter_array[17] #Star temperature
    
#         #Gas
#         sim.ini.gas.alpha=parameter_array[18] #alpha viscosity parameter
#         sim.ini.gas.Mdisk=parameter_array[19]*(1/sim.ini.dust.d2gRatio)*1.989e+33 #Inital gas disk mass
#         sim.ini.gas.mu=parameter_array[20] #Mean molecular weight of the gas
#         sim.ini.gas.SigmaExp=parameter_array[21] #Power law exponent of surface density profile
#         sim.ini.gas.SigmaRc=parameter_array[22] * (897587224200000.0/30) #Critical cut-off radius of surface density
        
#         print("Edited parameters are:")
#         sim.ini
        
#     else:
#         print("#Information: DustPy default parameters will be used#")
#         print("Parameters are:")
#         sim.ini
    
#     return(sim)

'---Refining inner grid cell---'
def refine_inner_grid_cell(DustPy_parameters):
    if DustPy_parameters.nr_in != 0 and (DustPy_parameters.nr_in/DustPy_parameters.sim.ini.grid.Nr)<=0.1:
        prYellow("#Information: The inner grid cell refinement has been activated, this may affect DustPy performance#")
        ###Code by R. Tazaki, based on MCFOST code by Christophe Pinte and contributors###
        def mcfost_make_r_walls(Rin, Rout, nr, nr_in):
            ln_delta_r = (1.0 / (nr - nr_in + 1.0)) * np.log(Rout / Rin)
            delta_r = np.exp(ln_delta_r)
        
            tab = np.empty(nr + 1, dtype=float)
            tab[0] = Rin
        
            denom = (2.0**nr_in - 1.0)
            for k in range(1, nr_in + 1):
                frac = (2.0**k - 1.0) / denom
                tab[k] = np.exp(np.log(Rin) - (np.log(Rin) - np.log(Rin * delta_r)) * frac)
        
            for k in range(nr_in + 1, nr + 1):
                tab[k] = tab[k - 1] * delta_r
        
            return tab
            
        def mcfost_make_r_centers(tab):
            tab = np.asarray(tab, dtype=float)
            return 0.5 * (tab[1:] + tab[:-1])
        
        r_i = mcfost_make_r_walls(Rin=DustPy_parameters.sim.ini.grid.rmin, Rout=DustPy_parameters.sim.ini.grid.rmax, nr=DustPy_parameters.sim.ini.grid.Nr, nr_in=DustPy_parameters.nr_in)
        r_c = mcfost_make_r_centers(r_i)
        print('cell walls:',r_i)
        print('cell centers:',r_c)
        DustPy_parameters.sim.grid.ri = r_i
        DustPy_parameters.sim.grid.r = r_c
    
    elif DustPy_parameters.nr_in != 0:
        prPurple("#ERROR: The inner grid cell refinement has been activated, but the refinement would spread the cell along more than 10% of the grid.\nThis would hinder DustPy performance too much. Reverting to the default: inner grid cell refinement deactivated#")

'---Bouncing from Dominik and Dullemond 2023---'
def activate_bouncing(sim,a_mono):
    prYellow("#Information: Bouncing has been activated in DustPy simulation#")
    ###Code and comments from Dominik and Dullemond 2023###
    
    #To implement the bouncing barrier into DustPy, we first (before calling sim.run()) computed the bouncing velocity v_{b}
    sim.dust.addfield("amonomer",a_mono,description="Monomer size (cm)")
    sim.dust.addfield("froll",1e-4,description="Force needed to roll a monomer") # Value from Heim et al. 1999
    sim.dust.addfield("mredu",sim.grid.m[:,None] * sim.grid.m[None,:] / ( sim.grid.m[:,None] + sim.grid.m[None,:] ),description="Reduced mass of two particles")
    sim.dust.v.addfield("bounce",np.sqrt(5*np.pi*sim.dust.amonomer*sim.dust.froll/sim.dust.mredu),description="Maximum sticking velocity [cm/s]") # Eq. 7 of Guettler et al. 2010

    #Next we overrode the sim.dust.p.stick() function of DustPy by defining our own
    def dd_p_stick(sim):
        dum = (sim.dust.v.bounce[None,:,:] / sim.dust.v.rel.tot)**2
        sim.dust.p.bounce = (1.5*dum + 1) * np.exp(-1.5*dum)
        pnostick = np.maximum(sim.dust.p.frag,sim.dust.p.bounce)
        p = 1. - pnostick
        p[0] = 0.
        p[-1] = 0.
        return p

    #And we linked it in through the following assignment
    sim.dust.p.stick.updater = dd_p_stick
    
'---Tentative to multithread---'
def activate_threading(sim):
    def multi_thread(sim):
        import threading
        import time
        
        def delta(frame):
            frame.dust.delta.update()
            time.sleep(1)
        def rhos(frame):
            frame.dust.rhos.update()
            time.sleep(1)
        def fill(frame):
            frame.dust.fill.update()
            time.sleep(1)
        def a(frame):
            frame.dust.a.update()
            time.sleep(1)
        def St(frame):
            frame.dust.St.update()
            time.sleep(1)
        def H(frame):
            frame.dust.H.update()
            time.sleep(1)
        def rho(frame):
            frame.dust.rho.update()
            time.sleep(1)
        def backreaction(frame):
            frame.dust.backreaction.update()
            time.sleep(1)
        def v(frame):
            frame.dust.v.update()
            time.sleep(1)
        def D(frame):
            frame.dust.D.update()
            time.sleep(1)
        def eps(frame):
            frame.dust.eps.update()
            time.sleep(1)
        def kernel(frame):
            frame.dust.kernel.update()
            time.sleep(1)
        def p(frame):
            frame.dust.p.update()
            time.sleep(1)
        def S(frame):
            frame.dust.S.update()
            time.sleep(1)
        
        t1 = threading.Thread(target=delta, args=sim)
        t2 = threading.Thread(target=rhos, args=sim)
        t3 = threading.Thread(target=fill, args=sim)
        t4 = threading.Thread(target=a, args=sim)
        t5 = threading.Thread(target=St, args=sim)
        t6 = threading.Thread(target=H, args=sim)
        t7 = threading.Thread(target=rho, args=sim)
        t8 = threading.Thread(target=backreaction, args=sim)
        t9 = threading.Thread(target=v, args=sim)
        t10 = threading.Thread(target=D, args=sim)
        t11 = threading.Thread(target=eps, args=sim)
        t12 = threading.Thread(target=kernel, args=sim)
        t13 = threading.Thread(target=p, args=sim)
        t14 = threading.Thread(target=S, args=sim)
    
        t1.start()
        t2.start()
        t3.start()
        t4.start()
        t5.start()
        t6.start()
        t7.start()
        t8.start()
        t9.start()
        t10.start()
        t11.start()
        t12.start()
        t13.start()
        t14.start()
        t1.join()
        t2.join()
        t3.join()
        t4.join()
        t5.join()
        t6.join()
        t7.join()
        t8.join()
        t9.join()
        t10.join()
        t11.join()
        t12.join()
        t13.join()
        t14.join()

    sim.dust.updater = multi_thread
    
'---Runing DustPy---'
def run_DustPy_function(DustPy_parameters):
    
    refine_inner_grid_cell(DustPy_parameters) #Extends the first cell onto other cells to match MCFOST
    DustPy_parameters.sim.initialize() #Initialiazes simulation object with chosen parameters
    
    if DustPy_parameters.bouncing==True: #Activates bouncing if chosen
        activate_bouncing(DustPy_parameters.sim,DustPy_parameters.a_monomer)
    
    # activate_threading(DustPy_parameters.sim)
        
    DustPy_parameters.sim.t.snapshots = np.hstack([DustPy_parameters.sim.t, np.geomspace(DustPy_parameters.starting_time, DustPy_parameters.stopping_time, num=DustPy_parameters.slice_number) * c.year]) #Choses starting time, stopping time, number of slices
    DustPy_parameters.sim.writer.datadir = DustPy_parameters.path #Choses output folder
    
    prGreen("======Runing DustPy======")
    DustPy_parameters.sim.run() #Runs the simulation
    
    if DustPy_parameters.plot==True:
        plot_DustPy_results(DustPy_parameters.path)
    
'---Plotting the simulation results---'
    
def plot_DustPy_results(path):
    prGreen("---Plotting DustPy results stored in "+path+" as an interactive window---")
    
    plot.ipanel(path)
