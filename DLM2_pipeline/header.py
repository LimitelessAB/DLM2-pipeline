"""
This file defines the classes used by the rest of the code

Developped by A. Blazère in June 2026
"""

'---Initialisation---'
import numpy as np
from dustpy import Simulation
from pymcfost.parameters import Params, _word_to_bool

'---Class definition---'

#A new class used for storing usefull DustPy parameters
class DustPy_parameters:
    #A function to define this new class
    def __init__(self):
            self.sim = Simulation() #Creates a simulation object that has DustPy default parameters
            self.nr_in = 0
            self.bouncing = False #Bouncing?
            self.a_monomer = 1e-4 #Monomer size in cm
            self.starting_time = 1e3 #Starting time in s
            self.stopping_time = 1e5 #Stopping time in s
            self.slice_number = 21 #Number of time slices
            self.path = "." #Path to save DustPy results, MCFOST results, and dneb measurments to
            self.plot = False #Plot DustPy results?
    
    #A function to print this new class
    def __str__(self):
        print(self.sim.ini.grid)
        print("The inner most radial cell will be extended on to "+str(self.nr_in)+" cells")
        print(self.sim.ini.star)
        print(self.sim.ini.gas)
        print(self.sim.ini.dust)
        print("Bouncing on? "+str(self.bouncing))
        print("a_monomer="+str(self.a_monomer)+" cm")
        print("Starting time="+str(self.starting_time)+ "s")
        print("Stopping time="+str(self.stopping_time)+ "s")
        print("Path="+str(self.path))
        print("Plotting results? "+str(self.plot))
        return("===End of class===")
    
    #A function to change values in this new class using a list, where None means leave the default parameter
    def set_variables(self,parameter_array):
        
        #Dust
        if parameter_array[0]!=None: self.sim.ini.dust.aIniMax=parameter_array[0] #Maximum particle size that will be filled initially (as in what is the max size that exits initially)
        if parameter_array[1]!=None: self.sim.ini.dust.allowDriftingParticles=parameter_array[1] #Radial drift on/off
        if parameter_array[2]!=None: self.sim.ini.dust.erosionMassRatio=parameter_array[2] #This parameter defines the threshold between full fragmentation and erosion
        if parameter_array[3]!=None: self.sim.ini.dust.d2gRatio=parameter_array[3] #Dust to gas ratio
        if parameter_array[4]!=None: self.sim.ini.dust.distExp=parameter_array[4] #Inital size distribution of particles
        if parameter_array[5]!=None: self.sim.ini.dust.excavatedMass=parameter_array[5] #If errosive collision, how much mass is chipped of larger collision partner
        if parameter_array[6]!=None: self.sim.ini.dust.fragmentDistribution=parameter_array[6] #Whenever a particle collision produces fragments the fragments will have a mass distribution following a power law of index...
        if parameter_array[7]!=None: self.sim.ini.dust.rhoMonomer=parameter_array[7] #Monomer bulk mass density, or dust grain individual density
        if parameter_array[8]!=None: self.sim.ini.dust.vFrag=parameter_array[8] #Fragmentation velocity, if above erosion or fragmentation, if below sticking try increasing by an order of magnitude default is 1m/s (written 100 cm/s)
    
        #Grid
        if parameter_array[9]!=None: self.sim.ini.grid.Nmbpd=parameter_array[9] #Number of mass bins per mass decade
        if parameter_array[10]!=None: self.sim.ini.grid.mmin=(4/3) * np.pi * ((parameter_array[10]*1e-4)**3) * self.sim.ini.dust.rhoMonomer #Minimum particle mass
        if parameter_array[11]!=None: self.sim.ini.grid.mmax=(4/3) * np.pi * ((parameter_array[11]*1e-4)**3) * self.sim.ini.dust.rhoMonomer #Maximum particle mass
        if parameter_array[12]!=None: self.sim.ini.grid.Nr=parameter_array[12] #Number of radial grid cells
        if parameter_array[13]!=None: self.sim.ini.grid.rmin=parameter_array[13]*14959787070000.0 #Location of inner radial grid boundary
        if parameter_array[14]!=None: self.sim.ini.grid.rmax=parameter_array[14]*14959787070000.0 #Location of outer radial grid boundary
        if parameter_array[15]!=None: self.nr_in=parameter_array[15]
    
        #Star
        if parameter_array[16]!=None: self.sim.ini.star.M=parameter_array[16]*1.989e+33 #Star mass
        if parameter_array[17]!=None: self.sim.ini.star.R=parameter_array[17]*6.957e+10 #Star radius
        if parameter_array[18]!=None: self.sim.ini.star.T=parameter_array[18] #Star temperature
    
        #Gas
        if parameter_array[19]!=None: self.sim.ini.gas.alpha=parameter_array[19] #alpha viscosity parameter
        if parameter_array[20]!=None: self.sim.ini.gas.Mdisk=parameter_array[20]*(1/self.sim.ini.dust.d2gRatio)*1.989e+33 #Inital gas disk mass
        if parameter_array[21]!=None: self.sim.ini.gas.mu=parameter_array[21] #Mean molecular weight of the gas
        if parameter_array[22]!=None: self.sim.ini.gas.SigmaExp=parameter_array[22] #Power law exponent of surface density profile
        if parameter_array[23]!=None: self.sim.ini.gas.SigmaRc=parameter_array[23] * (897587224200000.0/30) #Critical cut-off radius of surface density
        
        #Bouncing
        if parameter_array[24]!=None: self.bouncing = parameter_array[24] #Bouncing ?
        if parameter_array[25]!=None: self.a_monomer = parameter_array[25] #Monomer size in cm
        
        #Sim para
        if parameter_array[26]!=None: self.starting_time = parameter_array[26]
        if parameter_array[27]!=None: self.stopping_time = parameter_array[27]
        if parameter_array[28]!='': self.path = parameter_array[28]
        if parameter_array[29]!=None: self.plot = parameter_array[29]
    
    #A function to change values in the parameters txt file
    def update_txt_file(self):
        labels=np.array([": Maximum particle size that will be filled initially (as in what is the max size that exits initially) in cm",": Radial drift on/off (1/0)",": This parameter defines the threshold between full fragmentation and erosion",
                    ": Dust to gas ratio",": Inital size distribution of particles",": If errosive collision, how much mass is chipped of larger collision partner",": Whenever a particle collision produces fragments the fragments will have a mass distribution following a power law of index",
                    ": Monomer bulk mass density, or dust grain individual density",": Fragmentation velocity (If above this velocity erosion or fragmentation. If below sticking.), default is 1m/s (written 100 cm/s)",
                    ": Number of mass bins per mass decade",": Minimum particle size in microns",": Maximum particle size in microns",": Number of radial grid cells",": Location of inner radial grid boundary in AU",": Location of outer radial grid boundary in AU", ": The inner most radial cell will be extended on to (cells) (0 will turn this off)",
                    ": Star mass in solar mass",": Star radius in solar radius",": Star temperature in K",": alpha viscosity parameter",": Inital dust disk mass in solar mass",": Mean molecular weight of the gas",": Power law exponent of surface density profile",
                    ": Critical cut-off radius of surface density in AU",": Bouncing on/off (1/0) (introduces a state between sticking and erosion or fragmentation)",": : Monomer size in cm (usually the smallest grain you have)",": Starting time in s",": Stopping time in s",
                    ": Path",": Plotting results on/off (1/0):"])
        parameter_array=np.array([self.sim.ini.dust.aIniMax, self.sim.ini.dust.allowDriftingParticles, self.sim.ini.dust.erosionMassRatio, self.sim.ini.dust.d2gRatio, self.sim.ini.dust.distExp, self.sim.ini.dust.excavatedMass, 
                self.sim.ini.dust.fragmentDistribution, self.sim.ini.dust.rhoMonomer, self.sim.ini.dust.vFrag, self.sim.ini.grid.Nmbpd, (((3/4) * (self.sim.ini.grid.mmin) / (np.pi*self.sim.ini.dust.rhoMonomer))**(1/3))*1e4,
                (((3/4) * (self.sim.ini.grid.mmax) / (np.pi*self.sim.ini.dust.rhoMonomer))**(1/3))*1e4,self.sim.ini.grid.Nr,self.sim.ini.grid.rmin/14959787070000.0,self.sim.ini.grid.rmax/14959787070000.0,
                self.nr_in,self.sim.ini.star.M/1.989e+33,self.sim.ini.star.R/6.957e+10,self.sim.ini.star.T,self.sim.ini.gas.alpha,self.sim.ini.gas.Mdisk/((1/self.sim.ini.dust.d2gRatio)*1.989e+33),
                self.sim.ini.gas.mu,self.sim.ini.gas.SigmaExp,self.sim.ini.gas.SigmaRc/(897587224200000.0/30),self.bouncing,self.a_monomer,self.starting_time,self.stopping_time,self.path,self.plot]).astype(str)
        np.savetxt(self.path+"/DustPy_params.txt",np.swapaxes([parameter_array,labels],0,1),fmt="%s")
        
        return(parameter_array)
            
    #A function to change values in this new class using a txt file
    def read_txt_file(self):
        try:
            data = np.loadtxt(self.path+"/DustPy_params.txt",dtype=str)
            parameter_array = float(data[:,0])
            self.set_variables(parameter_array)
        except FileNotFoundError: #If the txt file does not exist it is created using default parameters
            parameter_array=self.update_txt_file()
        
        return(parameter_array)      
    
    #A function to change values in this new class using an interface
    def choose_variables(self):  
        
        default_parameter_array=self.read_txt_file()
        
        from tkinter import Tk,Label,Button,W,Entry, Scrollbar,RIGHT,BOTH,X,Y,Canvas,Frame,NW,ALL,TOP

        def onFrameConfigure(canvas): #Allows scrolling on the GUI
            canvas.configure(scrollregion=canvas.bbox(ALL))
        
        root = Tk() #Initiates GUI
        root.title("2DMLpipeline DustPy parameters selector window") #Sets name
        root.geometry('1000x1000') #Sets default size
        
        scrollbarx = Scrollbar(root) #Creates a horizontal scrollbar
        scrollbarx.pack(side=TOP, fill=X)
        
        scrollbary = Scrollbar(root) #Creates a vertical scrollbar
        scrollbary.pack(side=RIGHT, fill=Y)
        
        canvas = Canvas(root) #Creates the canvas on will our GUI will be displayed. It is separate from the root to allow scrolling
        canvas.configure(yscrollcommand=scrollbary.set)
        
        #Attaches scrollabar to canvas
        canvas.configure(xscrollcommand=scrollbarx.set)
        canvas.pack(side=RIGHT, fill=BOTH,expand = True)
        
        scrollbarx.configure(command=canvas.xview)
        scrollbary.configure(command=canvas.yview)
        
        #Updates what is shown when scolling
        content = Frame(canvas) 
        content.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
        
        canvas.create_window((0,0), window=content, anchor=NW)
        
        #Create fields and confirmation button that when pressed collects data
        labels=["Maximum particle size that will be filled initially (as in what is the max size that exits initially) in cm:","Radial drift on/off (1/0):","This parameter defines the threshold between full fragmentation and erosion:",
                    "Dust to gas ratio:","Inital size distribution of particles:","If errosive collision, how much mass is chipped of larger collision partner:","Whenever a particle collision produces fragments the fragments will have a mass distribution following a power law of index:",
                    "Monomer bulk mass density, or dust grain individual density:","Fragmentation velocity (If above this velocity: erosion or fragmentation. If below: sticking.), default is 1m/s (written 100 cm/s):",
                    "Number of mass bins per mass decade:","Minimum particle size in microns:","Maximum particle size in microns:","Number of radial grid cells:","Location of inner radial grid boundary in AU:","Location of outer radial grid boundary in AU:", "The inner most radial cell will be extended on to (cells) (0 will turn this off):",
                    "Star mass in solar mass:","Star radius in solar radius:","Star temperature in K:","alpha viscosity parameter:","Inital dust disk mass in solar mass:","Mean molecular weight of the gas:","Power law exponent of surface density profile:",
                    "Critical cut-off radius of surface density in AU:","Bouncing on/off (1/0) (introduces a state between sticking and erosion or fragmentation):","Monomer size in microns (usually the smallest grain you have):","Starting time in s:","Stopping time in s:",
                    "Path:","Plotting results on/off (1/0):","Parameters:","Your inputs, blancks will be left as DustPy defaults:"]
        textboxes=[None] * 30
        values=[None] * 30
        lbl = Label(content, text = labels[30], justify="left", font='Helvetica 18 bold')
        lbl.grid(sticky = W,column =0, row =0)
        lbl = Label(content, text = labels[31], justify="left", font='Helvetica 18 bold')
        lbl.grid(sticky = W,column =1, row =0)
        for i in range(0,30):
            lbl = Label(content, text = labels[i], justify="left")
            lbl.grid(sticky = W,column =0, row =i+1)
            textboxes[i] = Entry(content, width=10)
            if i==1 or i==24 or i==29: textboxes[i].insert(0, str(int(default_parameter_array[i])))
            else: textboxes[i].insert(0, str(default_parameter_array[i]))
            textboxes[i].grid(column =1, row =i+1)
        def clicked():
            for i in range(0,30):
                if i==1 or i==24 or i==29: #We want these fields to be converted to bool
                    try:
                        if textboxes[i].get()=="":
                            continue
                        else:
                            values[i]=bool(float(textboxes[i].get()))
                            continue
                    except ValueError:
                        print("#ERROR: {"+labels[i]+"} value was given in an incompatible format. It is supposed to be a float. Keeping DustPy default value.")
                    except TypeError:
                        print("#ERROR: {"+labels[i]+"} value was given in an incompatible format. It is supposed to be a float. Keeping DustPy default value.")
                elif i==27: values[i]=textboxes[i].get() #We want this field to be converted to str
                else: #We want these field to be converted to float
                    try:
                        if textboxes[i].get()=="":
                            continue
                        else:
                            values[i]=float(textboxes[i].get())
                            continue
                    except ValueError:
                        print("#ERROR: {"+labels[i]+"} value was given in an incompatible format. It is supposed to be a float. Keeping DustPy default value.")
                    except TypeError:
                        print("#ERROR: {"+labels[i]+"} value was given in an incompatible format. It is supposed to be a float. Keeping DustPy default value.")
            root.destroy()
        btn = Button(content,text = "Confirm parameters",font='Helvetica 18 bold',fg = "red", command=clicked)
        btn.grid(column=1, row=30)
        root.mainloop()
        
        self.set_variables(values)
        self.update_txt_file()
        
#A new class used for storing usefull DustPy parameters
class MCFOST_parameters:
    #A function to define this new class
    def __init__(self,path=None):
            self.path=path #Where to find MCFOST para file. 
            if self.path==None: #If there's no MCFOST para file, than this cannot work. 
                print("###FATAL ERROR: User cannot initialise this class with no MCFOST para file, if you want a default para file pls run $mcfost -get_para, in your terminal.\nFYI: The files are usually downloaded to the root.###")
            else:
                self.parameters=Params(path)
     
    #A function to print this new class
    def __str__(self):
        print("The following parameters come from a parameter file stored in :"+self.path)
        print(self.parameters)
        return("===End of class===")
    
    #A function to update the class parameters from the DustPy outputs, and chosen r_in
    def update_parameters(self,parameter_array,n_rin):
        cm_to_AU=6.68459e-14
        g_to_solar_mass=1/1.989e+33
        # -- Grid --
        self.parameters.grid.n_rad = int(parameter_array[12])
        self.parameters.grid.nz = int(parameter_array[13])
        self.parameters.grid.n_rad_in = int(n_rin)
        # -- Density structure --
        for k in range(0, self.parameters.simu.n_zones):
            self.parameters.zones[k].dust_mass = parameter_array[7] * g_to_solar_mass
            self.parameters.zones[k].gas_to_dust_ratio = 1/parameter_array[8]
            self.parameters.zones[k].h0 = parameter_array[5] * cm_to_AU
            self.parameters.zones[k].Rref = parameter_array[4] * cm_to_AU
            self.parameters.zones[k].Rin = parameter_array[3] * cm_to_AU
            self.parameters.zones[k].Rout = parameter_array[4] * cm_to_AU
            self.parameters.zones[k].flaring_exp = parameter_array[6]
        # -- Grain properties --
        for k in range(0, self.parameters.simu.n_zones):
            for j in range(self.parameters.zones[k].n_species):       
                self.parameters.zones[k].dust[j].amin = parameter_array[9]
                self.parameters.zones[k].dust[j].amax = parameter_array[10]
                self.parameters.zones[k].dust[j].n_grains = int(parameter_array[11])
        # -- Star properties --
        for k in range(0, self.parameters.simu.n_zones):
            self.parameters.stars[k].Teff = parameter_array[0]
            self.parameters.stars[k].R = parameter_array[1] * cm_to_AU
            self.parameters.stars[k].M = parameter_array[2] * g_to_solar_mass
            self.parameters.stars[k].file = "lte"+str(int(round(self.parameters.stars[k].Teff*1e-2)*1e2))+"-4.0.NextGen.fits.gz"
            
    #A function to update the para file.
    def update_para_file(self):
        self.parameters.writeto(self.path)
        
# =============================================================================
#     def update_para_file(self,parameters_array):
#         Parameters.phot.nphot_T = parameters_array[0]
#         Parameters.phot.nphot_SED = parameters_array[1]
#         Parameters.phot.nphot_image = parameters_array[2]
#         # -- Wavelengths --
#         Parameters.wavelengths.n_wl = parameters_array[3]
#         Parameters.wavelengths.wl_min = parameters_array[4]
#         Parameters.wavelengths.wl_max = parameters_array[5]
#         Parameters.simu.compute_T = _word_to_bool(parameters_array[6])
#         Parameters.simu.compute_SED = _word_to_bool(parameters_array[7])
#         Parameters.simu.use_default_wl = _word_to_bool(parameters_array[8])
#         Parameters.wavelengths.file = parameters_array[9]
#         Parameters.simu.separate_contrib = _word_to_bool(parameters_array[9])
#         Parameters.simu.separate_pola = _word_to_bool(parameters_array[10])
#         # -- Grid --
#         Parameters.grid.type = parameters_array[11]
#         Parameters.grid.n_rad = parameters_array[12]
#         Parameters.grid.nz = parameters_array[13]
#         Parameters.grid.n_az = parameters_array[14]
#         Parameters.grid.n_rad_in = parameters_array[15]
#         # -- Maps --
#         Parameters.map.nx = parameters_array[16]
#         Parameters.map.ny = parameters_array[17]
#         Parameters.map.size = parameters_array[18]
#         Parameters.map.RT_imin = parameters_array[19]
#         Parameters.map.RT_imax = parameters_array[20]
#         Parameters.map.RT_ntheta = parameters_array[21]
#         Parameters.map.lRT_centered = _word_to_bool(parameters_array[22])
#         Parameters.map.RT_az_min = parameters_array[23]
#         Parameters.map.RT_az_max = parameters_array[24]
#         Parameters.map.RT_n_az = parameters_array[25]
#         Parameters.map.distance = parameters_array[26]
#         Parameters.map.PA = parameters_array[27]
#         # -- Scattering method --
#         Parameters.simu.scattering_method = parameters_array[28]
#         # -- Symetries --
#         Parameters.simu.image_symmetry = _word_to_bool(parameters_array[29])
#         Parameters.simu.central_symmetry = _word_to_bool(parameters_array[30])
#         Parameters.simu.axial_symmetry = _word_to_bool(parameters_array[31])
#         # -- Disk physics --
#         Parameters.simu.dust_settling_type = parameters_array[32]
#         Parameters.simu.dust_settling_exp = parameters_array[33]
#         Parameters.simu.a_settling = parameters_array[34]
#         Parameters.simu.radial_migration = _word_to_bool(35)
#         Parameters.simu.dust_sublimation = _word_to_bool(parameters_array[36])
#         Parameters.simu.hydrostatic_eq = _word_to_bool(parameters_array[37])
#         Parameters.simu.viscous_heating = _word_to_bool(parameters_array[38])
#         Parameters.simu.viscosity = parameters_array[39]
#         # -- Number of zones --
#         n_zones = parameters_array[40]
#         Parameters.simu.n_zones = n_zones
#         # -- Density structure --
#         for k in range(n_zones):  
#             Parameters.zones[k].geometry = parameters_array[41]
#             Parameters.zones[k].dust_mass = parameters_array[42]
#             Parameters.zones[k].gas_to_dust_ratio = parameters_array[43]
#             Parameters.zones[k].h0 = parameters_array[44]
#             Parameters.zones[k].Rref = parameters_array[45]
#             Parameters.zones[k].vertical_exp = parameters_array[46]
#             Parameters.zones[k].Rin = parameters_array[47]
#             Parameters.zones[k].edge = parameters_array[48]
#             Parameters.zones[k].Rout = parameters_array[49]
#             Parameters.zones[k].Rc = parameters_array[50]
#             Parameters.zones[k].flaring_exp = parameters_array[51]
#             Parameters.zones[k].surface_density_exp = parameters_array[52]
#             Parameters.zones[k].m_gamma_exp = parameters_array[53]
#         # -- Grain properties --
#         for k in range(n_zones):
#             n_species = parameters_array[54]
#             Parameters.zones[k].n_species = n_species
#             for j in range(n_species):       
#                 Parameters.zones[k].dust[j].type = parameters_array[56]
#                 n_components = parameters_array[57]
#                 Parameters.zones[k].dust[j].n_components = n_components
#                 Parameters.zones[k].dust[j].mixing_rule = parameters_array[58]
#                 Parameters.zones[k].dust[j].porosity = parameters_array[59]
#                 Parameters.zones[k].dust[j].mass_fraction = parameters_array[60]
#                 Parameters.zones[k].dust[j].DHS_Vmax = parameters_array[61]
#                 for l in range(n_components):
#                     Parameters.zones[k].dust[j].component[l].file = parameters_array[62]
#                     Parameters.zones[k].dust[j].component[l].volume_fraction = parameters_array[63]
#                 Parameters.zones[k].dust[j].heating_method = parameters_array[64]
#                 Parameters.zones[k].dust[j].amin = parameters_array[65]
#                 Parameters.zones[k].dust[j].amax = parameters_array[66]
#                 Parameters.zones[k].dust[j].aexp = parameters_array[67]
#                 Parameters.zones[k].dust[j].n_grains = parameters_array[68]
#         -- Molecular settings --
#         Parameters.mol.compute_pop = _word_to_bool(parameters_array[69])
#         Parameters.mol.compute_pop_accurate = _word_to_bool(parameters_array[70])
#         Parameters.mol.LTE = _word_to_bool(parameters_array[71])
#         Parameters.mol.profile_width = parameters_array[72]
#         Parameters.mol.v_turb = parameters_array[73]
#         Parameters.mol.v_turb_unit = parameters_array[74]
#         n_mol = parameters_array[75]
#         Parameters.mol.n_mol = n_mol
#         for k in range(n_mol):
#             Parameters.mol.molecule.append(Molecule())
#             Parameters.mol.molecule[k].file = line[0]
#             Parameters.mol.molecule[k].level_max = int(line[1])
#             if (Parameters.simu.version < 4.1):
#                 Parameters.mol.molecule[k].v_max = float(line[0])
#                 Parameters.mol.molecule[k].v_min = - Parameters.mol.molecule[k].v_max
#                 Parameters.mol.molecule[k].nv = 2*int(line[1])+1
#             Parameters.mol.molecule[k].cst_abundance = _word_to_bool(line[0])
#             Parameters.mol.molecule[k].abundance = line[1]
#             Parameters.mol.molecule[k].abundance_file = line[2]
#             Parameters.mol.molecule[k].ray_tracing = _word_to_bool(line[0])
#             nTrans = int(line[1])
#             Parameters.mol.molecule[k].n_trans = nTrans
#             Parameters.mol.molecule[k].transitions = list(
#                 map(int, line[0:nTrans])
#             )  # convert list of str to int
#             if (Parameters.simu.version > 4.0):
#                 Parameters.mol.molecule[k].v_min = float(line[0])
#                 Parameters.mol.molecule[k].v_max = float(line[1])
#                 Parameters.mol.molecule[k].nv = int(line[2])
#         if (Parameters.simu.version > 3.0):
#             # -- Atom settings --
#             n_atoms = int(line[0])
#             Parameters.atomic.n_atoms = n_atoms
#             for k in range(n_atoms):
#                 Parameters.atomic.atom.append(Atom())
#                 Parameters.atomic.atom[k].file = line[0]
#                 Parameters.atomic.atom[k].nLTE = _word_to_bool(line[0])
#                 Parameters.atomic.atom[k].initial_solution = int(line[0])
#                 Parameters.atomic.atom[k].v_max = float(line[0])
#                 Parameters.atomic.atom[k].nv = int(line[1])
#                 Parameters.atomic.atom[k].images = _word_to_bool(line[0])
#                 n_trans = int(line[1])
#                 Parameters.atomic.atom[k].n_trans = n_trans
#                 Parameters.atomic.atom[k].lower = np.zeros(n_trans, dtype=int)
#                 Parameters.atomic.atom[k].upper = np.zeros(n_trans, dtype=int)
#                 for l in range(n_trans):
#                     Parameters.atomic.atom[k].lower[l] = int(line[0])
#                     Parameters.atomic.atom[k].upper[l] = int(line[1])
#         else:
#             Parameters.atomic.n_atoms = 0
#         -- Star properties --
#         n_stars = 1
#         Parameters.simu.n_stars = n_stars
#         for k in range(n_stars):
#             Parameters.stars[k].Teff = output[i,0]
#             Parameters.stars[k].R = output[i,1] * cm_to_AU
#             Parameters.stars[k].M = output[i,2] * g_to_solar_mass
#             Parameters.stars[k].x = 0
#             Parameters.stars[k].y = 0
#             Parameters.stars[k].z = 0
#             Parameters.stars[k].is_bb = _word_to_bool('F')
#             Parameters.stars[k].file = "lte"+str(np.round(Parameters.stars[k].Teff*1e-2)*1e2)+"-4.0.NextGen.fits.gz"
#             Parameters.stars[k].fUV = 0.04
#             Parameters.stars[k].slope_UV = 2.0
#             
#         Params.writeto(Parameters, path)
# =============================================================================
        
#A new class used for storing usefull PSF parameters
class PSF_parameters:
    #A function to define this new class
    def __init__(self,Nx_list,Ny_list,Pxscl_list,path):
        self.path=path #Folder where to find the fits of the PSFs
        self.Nx_list = Nx_list
        self.Ny_list = Ny_list
        self.Pxscl_list = Pxscl_list #arsec/pixel
    
    #A function to print this new class
    def __str__(self):
        print("The PSFs have "+str(self.Nx_list)+" pixels on the X axis")
        print("The PSFs have "+str(self.Ny_list)+" on the Y axis")
        print("The PSFs have a pixel scale of "+str(self.Pxscl_list)+" arsec/pixel")
        return("===End of class===")

#A new class used for storing usefull MCFOST images parameters
class MCFOST_images_parameters:
    #A function to define this new class
    def __init__(self,wavelength_list,Nx_list,Ny_list,Grid_size,distance):
        self.wavelength_list = wavelength_list #microns
        self.Nx_list = np.array(Nx_list, dtype=float)
        self.Ny_list = np.array(Ny_list, dtype=float)
        self.Grid_size = Grid_size #AU
        self.distance = distance #Pc
        self.Pxscl_list = ( self.Grid_size*1.496e+11 / (self.Nx_list * self.distance * 3.086e+16) ) * (180.0/np.pi)*(60.0**2)  #arsec  
        
    #A function to print this new class
    def __str__(self):
        print("The chosen material is file: "+self.material)
        print("The wavelegths to be explored are (microns): "+str(self.wavelength_list))
        print("The desired : "+str(self.wavelength_list))
        print("The chosen default sizes for the X axis are (pixels): "+str(self.Nx_list))
        print("The chosen default sizes for the Y axis are (pixels): "+str(self.Ny_list))
        print("The physical size of the image grid is: "+str(self.Grid_size)+" AU")
        print("The distance to the disk is: "+str(self.distace)+" Pc")
        return("===End of class===")