# =============================================================================
# Dustpy Linked to MCOFST for Dark Lane Measurment (DLM2) pipeline
# =============================================================================

"""
This file is the main file of DLM2 pipeline

Developped by A. Blazère in June 2026
"""

'---Initialisation---'
from .header import *
from .run_DustPy import *
from .run_MCFOST import *
from .measure_dark_lanes import *
from pymcfost import make_density_no_settling, make_density_parametric, make_density_Dubrulle, make_density_Fromang
import os

'---Colored prints---'
#Functions from https://www.geeksforgeeks.org/python/print-colors-python-terminal/
def prRed(s): print("\033[91m {}\033[00m".format(s))
def prGreen(s): print("\033[92m {}\033[00m".format(s))
def prYellow(s): print("\033[93m {}\033[00m".format(s))
def prLightPurple(s): print("\033[94m {}\033[00m".format(s))
def prPurple(s): print("\033[95m {}\033[00m".format(s))
def prCyan(s): print("\033[96m {}\033[00m".format(s))
def prLightGray(s): print("\033[97m {}\033[00m".format(s))

'---Runing all the components of 2DML pipeline---'
def master_function(DustPy_parameters,number_of_vertical_cells,Settling_type,MCFOST_parameters,DustPy_slices_list_to_run,MCFOST_image_parameters,PSF_parameters=None,convolution_triggered_list=None,parametric_parameters=None,skip_dustpy=None):
    prLightPurple("="*68)
    prLightPurple("||    DML2 pipeline code developped by A. Blazère in June 2026    ||")
    prLightPurple("="*68)
    print("")
    
    #Setting the new root so that MCFOST writes files in the correct place
    os.chdir(DustPy_parameters.path)
    
    #DustPy run
    if skip_dustpy==None: run_DustPy_function(DustPy_parameters)
    
    #DustPy file conversion to MCFOST density file format
    para_file_changes=[]
    match Settling_type:
        
        case "No settling":
            prGreen("======Converting DustPy linear density to 2D density map with no settling for MCFOST======")
            print("")
            
            para_file_changes=make_density_no_settling(DustPy_parameters.path,DustPy_parameters.slice_number,number_of_vertical_cells,1)
            
        case "Parametric":
            prGreen("======Converting DustPy linear density to 2D density map with parametric settling for MCFOST======")
            print("")
            
            if parametric_parameters==None: 
                return("###FATAL ERROR: User chose parametric settling but parametric settling parameters are missing. If you do not know what to chose, take a look at Duchêne et al. 2024 (https://doi.org/10.3847/1538-3881/acf9a7)###")
            else:
                para_file_changes=make_density_parametric(DustPy_parameters.path,DustPy_parameters.slice_number+1,number_of_vertical_cells,parametric_parameters,1) #With the data file 0 there is slice_number + 1 files
        
        case "Dubrulle":
            prGreen("======Converting DustPy linear density to 2D density map with Dubrulle settling for MCFOST======")
            print("")
            
            para_file_changes=make_density_Dubrulle(DustPy_parameters.path,DustPy_parameters.slice_number+1,number_of_vertical_cells,1)
        
        case "Fromang":
            prGreen("======Converting DustPy linear density to 2D density map with Fromang settling for MCFOST======")
            print("")
            
            para_file_changes=make_density_Fromang(DustPy_parameters.path,DustPy_parameters.slice_number+1,number_of_vertical_cells,1)
            
        case _:
            prRed("The settling type entered is incorrect. Possible options are: No settling, Parametric, Dubrulle, Fromang.")
            return("###Fatal error: Unable to compute vertical density###")
    
    if convolution_triggered_list==None: convolution_triggered_list=[[False] * len(MCFOST_image_parameters.wavelength_list)] * len(DustPy_slices_list_to_run)
    var=-1
    for slice in DustPy_slices_list_to_run:
        var+=1
        #Creating the MCFOST para file
        MCFOST_parameters.update_parameters(para_file_changes[slice],DustPy_parameters.nr_in)
        MCFOST_parameters.path=DustPy_parameters.path+"/data"+str(slice)+".para"
        MCFOST_parameters.update_para_file()
        
        #Runing MCFOST on the last para file
        complete_run_MCFOST_function(slice,DustPy_parameters.path,Settling_type,MCFOST_parameters.path,convolution_triggered_list[var],DustPy_parameters.sim.ini.grid.rmax*6.68459e-14,MCFOST_image_parameters,PSF_parameters)
    
    prLightPurple("="*46)
    prLightPurple("||    All processes have been completed!    ||")
    prLightPurple("="*46)
    
            
    