"""
This package contains function to configure and run MCFOST

Developped by A. Blazère in June 2026
"""

'---Initialisation---'
import pymcfost as mcfost
from pymcfost import Image
from pymcfost import run
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import scipy.ndimage
# from tp_utils import *

'---Colored prints---'
#Functions from https://www.geeksforgeeks.org/python/print-colors-python-terminal/
def prRed(s): print("\033[91m {}\033[00m".format(s))
def prGreen(s): print("\033[92m {}\033[00m".format(s))
def prYellow(s): print("\033[93m {}\033[00m".format(s))
def prLightPurple(s): print("\033[94m {}\033[00m".format(s))
def prPurple(s): print("\033[95m {}\033[00m".format(s))
def prCyan(s): print("\033[96m {}\033[00m".format(s))
def prLightGray(s): print("\033[97m {}\033[00m".format(s))
      
'---MCFOST image calculations---'
#A funtion that calculates how much to zoom, and reshape MCFOST images to match PSF pxscl and keep enough resolution for dneb measurments
def find_MCFOST_options(MCFOST_image_parameters,Rout,PSF_parameters,f=3): #f is how much space should the disk take in the image
    options_list=np.zeros((len(PSF_parameters.Pxscl_list),3))
    for i in range (0,len(PSF_parameters.Pxscl_list)):
        if PSF_parameters.Nx_list[i]!=None: #PSF list should match MCFOST image list, but if there is no PSF or no desired convolution, then just place None in the PSF list slots
            zoom=0
            reshape=0
            option=1
            C1 = MCFOST_image_parameters.Pxscl_list[i]/PSF_parameters.Pxscl_list[i]
            # C2 = PSF_parameters.Nx_list[i] / MCFOST_image_parameters.Nx_list[i]
            # C3 = f * Rout / (MCFOST_image_parameters.Grid_size * C1)
            
            # if C3>=1: #This case means that the resulting image needs to be enlarged to match PSF number of pixels
            #     zoom = 1/C1
            #     reshape = MCFOST_image_parameters.Nx_list[i] * C2
            # elif C3<1: #This case means that the resulting image needs to be shrinked to match PSF number of pixels, this might hinder resolution so we do not apply C3 and later enlarge the PSF
            #     zoom = 1 / (C1 * C3) 
            #     reshape = MCFOST_image_parameters.Nx_list[i] * C2 * C3
            zoom = C1
            reshape = MCFOST_image_parameters.Nx_list[i]
            if reshape > PSF_parameters.Nx_list[i]: option=2
            options_list[i,:]=[zoom,reshape,option]
        else:
            options_list[i,:]=[None,MCFOST_image_parameters.Nx_list[i],None]
    return(options_list)

'---Convolution util---'
###Code of the following functions by Éric Thiébaut###
def zeropad(arr, shape):
    """Zero-pad array ARR to given shape.

    The contents of ARR is approximately centered in the result."""
    rank = arr.ndim
    if len(shape) != rank:
        raise ValueError("bad number of dimensions")
    diff = np.asarray(shape) - np.asarray(arr.shape)
    if diff.min() < 0:
        raise ValueError("output dimensions must be larger or equal input dimensions")
    offset = diff//2
    z = np.zeros(shape, dtype=arr.dtype)
    if rank == 1:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        z[i0:n0] = arr
    elif rank == 2:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        i1 = offset[1]; n1 = i1 + arr.shape[1]
        z[i0:n0,i1:n1] = arr
    elif rank == 3:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        i1 = offset[1]; n1 = i1 + arr.shape[1]
        i2 = offset[2]; n2 = i2 + arr.shape[2]
        z[i0:n0,i1:n1,i2:n2] = arr
    elif rank == 4:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        i1 = offset[1]; n1 = i1 + arr.shape[1]
        i2 = offset[2]; n2 = i2 + arr.shape[2]
        i3 = offset[3]; n3 = i3 + arr.shape[3]
        z[i0:n0,i1:n1,i2:n2,i3:n3] = arr
    elif rank == 5:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        i1 = offset[1]; n1 = i1 + arr.shape[1]
        i2 = offset[2]; n2 = i2 + arr.shape[2]
        i3 = offset[3]; n3 = i3 + arr.shape[3]
        i4 = offset[4]; n4 = i4 + arr.shape[4]
        z[i0:n0,i1:n1,i2:n2,i3:n3,i4:n4] = arr
    elif rank == 6:
        i0 = offset[0]; n0 = i0 + arr.shape[0]
        i1 = offset[1]; n1 = i1 + arr.shape[1]
        i2 = offset[2]; n2 = i2 + arr.shape[2]
        i3 = offset[3]; n3 = i3 + arr.shape[3]
        i4 = offset[4]; n4 = i4 + arr.shape[4]
        i5 = offset[5]; n5 = i5 + arr.shape[5]
        z[i0:n0,i1:n1,i2:n2,i3:n3,i4:n4,i5:n5] = arr
    else:
        raise ValueError("too many dimensions")
    return z

def ifft(x): return np.real(np.fft.ifftn(x))

'---Convolution of MCFOST images with PSF---'
def convolution(file_number,wavelength,options,PSF_path,Image_path):
    prGreen("======Convolving "+str(wavelength)+" microns image with PSF======")
    print("")
    PSF=fits.open(PSF_path)[0].data
    MCFOST_image=fits.open(Image_path)[0].data
    if options[2]==1:
        MTF=np.fft.fftn(np.fft.ifftshift(PSF)) #If we don't ifftshift our convolution will make for corner images on image corners.
        Convolved=ifft(MTF*np.fft.fftn(MCFOST_image[0,0,0,:,:]))
            
    else:
        Resized_shifted_psf=np.fft.ifftshift(zeropad(PSF,(int(options[1]),int(options[1])))) #To preserver our disk resolution, we have made the image larger than the PSF. We therefore increase it's size.
        MTF=np.fft.fftn(Resized_shifted_psf)
        Convolved=ifft(MTF*np.fft.fftn(MCFOST_image[0,0,0,:,:]))
    
    return(Convolved)

'---Image scales from pixel to arsec---'
def px_to_arsec(px,Pxscl,Image_center):
    return (px-Image_center) * Pxscl

'---Image plotting---'
def plot_MCFOST_images(convolution_triggered,wavelength_list,Nb_pixels_list,Pxscl_list,MCFOST_image_path):
    plt.figure(figsize=(80,80))
    for i in range (0,len(wavelength_list)):
        ax=plt.subplot(round(len(wavelength_list)/10)+1,10,i+1)
        if convolution_triggered[i]==True:
            img=fits.open(MCFOST_image_path+"/data_"+str(wavelength_list[i])+"/convolved.fits")[0].data
        else:
            img=(fits.open(MCFOST_image_path+"/data_"+str(wavelength_list[i])+"/RT.fits.gz")[0].data)[0,0,0,:,:]
        
        plt.xlabel(r'$\Delta$X (arcsec)')
        plt.ylabel(r'$\Delta$Y (arcsec)')
        plt.imshow(np.sqrt(img),origin='lower',cmap="cividis",extent=[px_to_arsec(0,Pxscl_list[i],round(Nb_pixels_list[i]/2)),
                                                                      px_to_arsec(Nb_pixels_list[i],Pxscl_list[i],round(Nb_pixels_list[i]/2)),
                                                                      px_to_arsec(0,Pxscl_list[i],round(Nb_pixels_list[i]/2)),
                                                                      px_to_arsec(Nb_pixels_list[i],Pxscl_list[i],round(Nb_pixels_list[i]/2))])
        plt.title('Image at '+str(wavelength_list[i])+' micrometers',fontweight="bold")
        plt.colorbar()
    plt.savefig(MCFOST_image_path+"/MCFOST_images.png")   
    
    prYellow("#Information: Saving a plot of all MCFOST images#")
    print("")
    
    plt.show()

'---Dark lane fitting---'
###Code of the following block by Gaspard Duchêne, with edits from A. Blazère to measure radius###
def get_spines_poly(im,pxscl,off_as,width_as,cen_x,cen_y,length_as,maxoff_as,order,save_plot=0):    

    length = int(round(length_as / pxscl))
    width = int(round(width_as / pxscl))
    maxoff = int(round(maxoff_as / pxscl))
    offset = int(round(off_as / pxscl))

    x = np.arange(0-length,length+1,1)
    oversampl = 10
    xx = np.arange(0-length,length+1/oversampl,1/oversampl)

    if save_plot == 1:
        length_plot = int(3*length)
        x_plot = (np.arange(0-length_plot,length_plot+1,1)) * pxscl
        y_plot = np.mean(im[cen_y-length_plot:cen_y+length_plot+1,cen_x+offset-width:cen_x+offset+width+1],axis=1)
        plt.plot(x_plot,y_plot,'.',color='darkred',label='Observed')
        plt.xlabel('$\Delta$Y (arcsec)')
        plt.ylabel('Surface Brightness (arb. unit)')

    # top nebula
    peak_approx = np.argmax(im[cen_y:cen_y+maxoff,cen_x+offset])
    y = np.mean(im[cen_y+peak_approx-length:cen_y+peak_approx+length+1,cen_x+offset-width:cen_x+offset+width+1],axis=1)
    pars = np.polyfit(x,y,order)
    fit = np.zeros(len(xx))
    for i in range(order+1):
        fit += pars[i] * xx**(order-i)
    top_peak_loc, top_peak_bright = cen_y+peak_approx+xx[np.argmax(fit)], max(fit)

    if save_plot == 1:
        plt.plot((xx+peak_approx)*pxscl,fit,'k',label='Polynomial Fit')        

    # bottom nebula
    peak_approx = np.argmax(im[cen_y-maxoff:cen_y,cen_x+offset])
    y = np.mean(im[cen_y-maxoff+peak_approx-length:cen_y-maxoff+peak_approx+length+1,cen_x+offset-width:cen_x+offset+width+1],axis=1)
    pars = np.polyfit(x,y,order)
    fit = np.zeros(len(xx))
    for i in range(order+1):
        fit += pars[i] * xx**(order-i)
    bot_peak_loc, bot_peak_bright = cen_y-maxoff+peak_approx+xx[np.argmax(fit)], max(fit)

    if save_plot == 1:
        plt.plot((0-maxoff+peak_approx+xx)*pxscl,fit,'k')
        plt.legend()
        plt.savefig('vertprof_polyfit.eps')

    return top_peak_loc,bot_peak_loc,top_peak_bright, bot_peak_bright

def measure_dneb_and_radius(convolution_triggered,wavelengths,nb_px_list,physical_size,pxscl_list,path,f=3):
    radius_list=[]
    Dark_lanes=[]
    prGreen("======Measuring dneb and radius======")
    print("")
    fig1=plt.figure(1,figsize=(80,80))
    fig2=plt.figure(2,figsize=(80,80))
    fig3=plt.figure(3,figsize=(80,80))
    # Load images from a directory containing MCFOST outputs
    for wavelength_number in range(0,len(pxscl_list)):
        ax1=fig1.add_subplot(round(len(pxscl_list)/10)+1,10,wavelength_number+1)
        ax2=fig2.add_subplot(round(len(pxscl_list)/10)+1,10,wavelength_number+1)
        ax3=fig3.add_subplot(round(len(pxscl_list)/10)+1,10,wavelength_number+1)
        if convolution_triggered==True:
            img_to_fit = (fits.open(path+"/data_"+str(wavelengths[wavelength_number])+'/convolved.fits')[0].data)
        else:
            img_to_fit = (fits.open(path+"/data_"+str(wavelengths[wavelength_number])+'/RT.fits.gz')[0].data)[0,0,0,:,:]
        
        # Input image properties
        pxscl_fit=pxscl_list[wavelength_number]
        xcen_init=round(nb_px_list[wavelength_number]/2)
        ycen_init=xcen_init

        # Parameters to vary for fit to evaluate uncertainties and extract only the spine
        order = 8 # Order of polynomial function to fit
        length_arcsec = 0.1 # Half width of vertical profile to fit, centered on local maximum [arcsec]
        window_arcsec = 0.2 # Half width of horizontal window to average the profile over [arcsec]
        maxoffset_arcsec = (1/f * nb_px_list[wavelength_number] * pxscl_fit) # Maximum distance along vertical direction to find the local maximum [arcsec]
        span = (1/f * nb_px_list[wavelength_number] * pxscl_fit)/2 # Maximum distance along horizontal direction to cover [arcsec]
        step = pxscl_fit # Step in horizontal direction [arcsec]

        xrange_2ndfit = (1/f * nb_px_list[wavelength_number] * pxscl_fit)/2 # Range to use in fitting polynomial to spine [arcsec]
        xavoid_bottom = 0 # Range to avoid around the symmetry axis for the bottom nebula [arcsec]
        xavoid_top = 0 # Range to avoid around the symmetry axis for the top nebula [arcsec] 
        order_2ndfit = 2 # Polynomial order for spine fitting
        sampl = pxscl_fit/10 # Sampling for finer polynomial evaluation and finding minimum distance [arcsec]
        
        offset_arcsec = np.arange(0-span,span+step,step)
        noff = len(offset_arcsec)
        ytop = np.zeros(noff)
        ybot = np.zeros(noff)
        ftop = np.zeros(noff)
        fbot = np.zeros(noff)

        for i in range(noff):
            yt, yb, ft, fb = get_spines_poly(img_to_fit,pxscl_fit,offset_arcsec[i],window_arcsec,xcen_init,ycen_init,length_arcsec,maxoffset_arcsec,order)
            ytop[i], ybot[i], ftop[i], fbot[i] = yt, yb, ft, fb

        ax1.imshow(np.log10(img_to_fit-1.01*np.min(img_to_fit)),origin='lower',vmin=-20,vmax=-16,cmap='gray_r',extent=[px_to_arsec(0,pxscl_fit,xcen_init),
                                                                                                                        px_to_arsec(nb_px_list[wavelength_number],pxscl_fit,xcen_init),
                                                                                                                        px_to_arsec(0,pxscl_fit,xcen_init),
                                                                                                                        px_to_arsec(nb_px_list[wavelength_number],pxscl_fit,xcen_init)])
        ax1.plot(offset_arcsec,px_to_arsec(ytop, pxscl_fit, ycen_init),'+',color='indianred')
        ax1.plot(offset_arcsec,px_to_arsec(ybot, pxscl_fit, ycen_init),'+',color='indianred')
        #plt.xlim([xcen_init-fov_fig,xcen_init+fov_fig])
        #plt.ylim([ycen_init-fov_fig,ycen_init+fov_fig])
        
        ax1.set_xlabel(r'$\Delta$X (arcsec)')
        ax1.set_ylabel(r'$\Delta$Y (arcsec)')
        plt.tight_layout()

        # plt.savefig('spine_illustr_f444w.eps')
        # plt.show()

        res = np.array([offset_arcsec,ytop*pxscl_fit,ybot*pxscl_fit,ftop,fbot])

        xx = np.arange(0.-xrange_2ndfit,xrange_2ndfit+sampl,sampl)
        fit_top = np.zeros(len(xx))
        fit_bot = np.zeros(len(xx))


        x = res[0,:]
        top = res[1,:]
        bot = res[2,:]
        ax2.plot(x,ftop)
        ax2.set_xlabel(r'$\Delta$X (arcsec)')
        ax2.set_ylabel(r'$Intensity (W\cdot m_{-2}\cdot pixel_{-1})$')
        plt.tight_layout()
        temp=[]
        nan_purged_ftop=np.nan_to_num(ftop)
        for n in range (0,len(nan_purged_ftop)):
            if ftop[n]>0.1*np.max(nan_purged_ftop):
                temp.append(x[n])
        radius_list.append(abs(np.max(temp)-np.min(temp))/2)

        subset = np.where((np.abs(x) <= xrange_2ndfit) & (np.abs(x) >= xavoid_top))

        ax3.plot(x,top-np.min(top),'+')
        pars = np.polyfit(x[subset],top[subset],order_2ndfit)
        for j in range(order_2ndfit+1):
            fit_top += pars[j] * xx**(order_2ndfit-j)
        ax3.plot(xx,fit_top-np.min(top))

        subset = np.where((np.abs(x) <= xrange_2ndfit) & (np.abs(x) >= xavoid_bottom))
            
        ax3.plot(x,bot-np.min(top),'+')
        pars = np.polyfit(x[subset],bot[subset],order_2ndfit)
        for j in range(order_2ndfit+1):
            fit_bot += pars[j] * xx**(order_2ndfit-j)
        ax3.plot(xx,fit_bot-np.min(top))
        ax3.set_xlabel(r'$\Delta$X (arcsec)')
        ax3.set_ylabel(r'$\Delta$Y (arcsec) shifted so as to have the top spine minimum at zero')
        plt.tight_layout()
        dist_fit = fit_top - fit_bot
            
        Dark_lanes.append(np.min(dist_fit))
    plt.show()
    fig1.savefig(path+"/Top_and_bottom_fits.png")
    fig2.savefig(path+"/Radius_fits.png")
    fig3.savefig(path+"/Spine_fits.png")
    prYellow("#Information: Saved all three figures as png files#")
    print("")
    return(Dark_lanes,radius_list)

'---Runing MCFOST and measuring dark lanes---'
def complete_run_MCFOST_function(file_number,DustPy_path,settling_type,para_path,convolution_triggered,Rout,MCFOST_image_parameters,PSF_parameters=None):
    
    prGreen("======Runing MCFOST on DustPy slice "+str(file_number)+"======")
    print("")
    
    if convolution_triggered != ([False]*len(MCFOST_image_parameters.wavelength_list)): #Convolution implies same pixel scale in MCFOST_image and PSF
        options_list=find_MCFOST_options(MCFOST_image_parameters, Rout, PSF_parameters)
        
    prGreen("---Rendering temperature map---")
    print("")
    
    run(para_path,options="-root_dir data"+str(file_number)+" -density_file "+DustPy_path+"/data"+str(file_number)+"_density_"+str(settling_type)+"_settling_file.fits") # Run MCFOST with a parameter file to create temperature file
    for wavelength_number in range (0, len(MCFOST_image_parameters.wavelength_list)):
        
        if PSF_parameters.Nx_list[wavelength_number]==None: #Because we need a pxscl for every image, even the ones that the user doesn't want convolved, we steal for those images the MCFOST_image_parameters pxscl and use it to replace the None in the PSF_parameters class.
            PSF_parameters.Pxscl_list[wavelength_number]=MCFOST_image_parameters.Pxscl_list[wavelength_number]
            
        prGreen("---Rendering "+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+" microns image---")
        print("")
        
        if convolution_triggered[wavelength_number]==True:
            #Run MCFOST with a parameter file to create image
            # prYellow("#User asked for a resolution of "+str(MCFOST_image_parameters.Nx_list[wavelength_number])+" by " +str(MCFOST_image_parameters.Nx_list[wavelength_number])+
            #          ".\nBut in order to convolve with PSF, switching to a resolution of "+str(options_list[wavelength_number,1])+" by " +str(options_list[wavelength_number,1])+".\nAnd applying a zoom of "+str(options_list[wavelength_number,0])+"#")
            # print("")   
            prYellow("#Information: User asked for a resolution of "+str(MCFOST_image_parameters.Nx_list[wavelength_number])+" by " +str(MCFOST_image_parameters.Nx_list[wavelength_number])+
                     ".To match PSF pxscl, applying a zoom of "+str(options_list[wavelength_number,0])+"#")
            print("") 
            
            if options_list[wavelength_number,2]==1 or options_list[wavelength_number,2]==2:
                #Run an MCFOST image applying a zoom to match PSF pxscl and changing the image size from MCFOST_Image class
                print(para_path)
                print("-root_dir data"+str(file_number)+" -density_file "+DustPy_path+"/data"+str(file_number)+"_density_"+str(settling_type)+
                    "_settling_file.fits -img "+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+" -resol "+
                    str(int(options_list[wavelength_number,1]))+" "+str(int(options_list[wavelength_number,1]))+" -zoom "+
                    str(options_list[wavelength_number,0]))
                run(para_path,options="-root_dir data"+str(file_number)+" -density_file "+DustPy_path+"/data"+str(file_number)+"_density_"+str(settling_type)+
                    "_settling_file.fits -img "+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+" -resol "+
                    str(int(options_list[wavelength_number,1]))+" "+str(int(options_list[wavelength_number,1]))+" -zoom "+
                    str(options_list[wavelength_number,0]))
                
            else:
                #Run an MCFOST image just applying a zoom to match PSF pxscl 
                run(para_path,options="-root_dir data"+str(file_number)+" -density_file "+DustPy_path+"/data"+str(file_number)+"_density_"+str(settling_type)+
                    "_settling_file.fits -img "+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+" -zoom "+str(options_list[wavelength_number,0]))
                
            #Convolve MCFOST image
            Convolved_image=convolution(file_number,MCFOST_image_parameters.wavelength_list[wavelength_number],options_list[wavelength_number,:],
                                        PSF_parameters.path+"/PSF-"+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+".fits",DustPy_path+"/data"+
                                        str(file_number)+"/data_"+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+'/RT.fits.gz')
            
            #Save convolved image
            convolved_hdr = fits.Header()
            convolved_hdu = fits.PrimaryHDU(Convolved_image, header=convolved_hdr)
            final = fits.HDUList([convolved_hdu])
            final.writeto(DustPy_path+"/data"+str(file_number)+"/data_"+str(MCFOST_image_parameters.wavelength_list[wavelength_number])+"/convolved.fits")
            
        else:
            #Run an MCFOST image
            run(para_path,options="-root_dir data"+str(file_number)+" -density_file "+DustPy_path+"/data"+str(file_number)+"_density_"+str(settling_type)+"_settling_file.fits -img "+str(MCFOST_image_parameters.wavelength_list[wavelength_number])) # Run MCFOST with a parameter file to create image
            
    if convolution_triggered != ([False]*len(MCFOST_image_parameters.wavelength_list)): #Plot and measure dneb and radius in the case where convolution is applied to at least one image
        plot_MCFOST_images(convolution_triggered,MCFOST_image_parameters.wavelength_list,options_list[:,1],PSF_parameters.Pxscl_list,DustPy_path+"/data"+str(file_number))
        Dark_lanes, radius_list = measure_dneb_and_radius(convolution_triggered,MCFOST_image_parameters.wavelength_list,options_list[:,1],(1/options_list[:,0])*MCFOST_image_parameters.Grid_size,PSF_parameters.Pxscl_list,DustPy_path+"/data"+str(file_number))
    else: #Plot and measure dneb and radius in the case where no convolution is applied to any image
        plot_MCFOST_images(convolution_triggered,MCFOST_image_parameters.wavelength_list,MCFOST_image_parameters.Nx_list,MCFOST_image_parameters.Pxscl_list,DustPy_path+"/data"+str(file_number))
        Dark_lanes, radius_list = measure_dneb_and_radius(convolution_triggered,MCFOST_image_parameters.wavelength_list,MCFOST_image_parameters.Nx_list,MCFOST_image_parameters.Grid_size,MCFOST_image_parameters.Pxscl_list,DustPy_path+"/data"+str(file_number))
    
    prGreen("======Plotting dneb and radius======")
    print("")
    
    #Plot the dneb values to a text file with space as the delimiter
    plt.plot(MCFOST_image_parameters.wavelength_list,Dark_lanes,label="dneb in arsec as a function of wavelength")
    plt.scatter(MCFOST_image_parameters.wavelength_list,Dark_lanes)
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("dneb (arsec)")
    plt.xscale('log')
    plt.legend()
    plt.show()
    #Save the dneb values to a text file with space as the delimiter
    
    #Plot the radius values to a text file with space as the delimiter
    np.savetxt(DustPy_path+"/data"+str(file_number)+"/dneb.txt",np.swapaxes([MCFOST_image_parameters.wavelength_list,Dark_lanes],0,1), delimiter=' ')
    plt.plot(MCFOST_image_parameters.wavelength_list,radius_list,label="Radius in arsec as a function of wavelength")
    plt.scatter(MCFOST_image_parameters.wavelength_list,radius_list)
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Radius (arsec)")
    plt.xscale('log')
    plt.legend()
    plt.show()
    #Save the radius values to a text file with space as the delimiter
    np.savetxt(DustPy_path+"/data"+str(file_number)+"/radius.txt",np.swapaxes([MCFOST_image_parameters.wavelength_list,radius_list],0,1), delimiter=' ')
    
    prYellow("#Information: Saved dneb and radius in a txt file#")