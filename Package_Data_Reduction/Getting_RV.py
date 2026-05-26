#--------------------------------------------------------------
# Packages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits

import os
import sys
from CrossCorrelation import CCF

#--------------------------------------------------------------
# Folder with the spectra
path_files = path_folder+EXTRA_FOLDER_RV

# List of Spectra (full path required)
list_file  = spectrums


#----------
# During my PhD I was testing a few things to get a good RV and I have left this option here that might work weel in some cases but I have not included the options to run it from INPUT.py
# Therefore just the default (using a syntetic spectra for all stars) is included here

# The method refers to if you want to use a synthetic spectra for all the stars ('Synth')
# Or if you just want to use it for the synthetic spectra for a template star and then use this rest frame star to get the RV of the others
method='Synthetic' # 'Synthetic', 'Template'



# Depending On the method used (not need on semiautomatic mode)

# Synthetic spectra
path_syn_spec  = complete_path_template # just to save this info

# Template Star
# path_temp_star = 'path_to_template_star'
# you have to include a complete path for this star (you can use any star, is not need to be of the same folder)


#-----------

# Just one work at the time
# PLOT_EACH_ONE= True # Defined on INPUT.py #Plot each one of the cross correlation and fit (it will display for 1 sec, you can changue the time in the function that calls the CCF)
PLOT_sumary  = not PLOT_EACH_ONE  # Plot a summary of the fits and a overploted spectra (to see if the features match in all the spectra)    

SAVE=True
SAVE_NAME='OUTPUT_RV.txt' # It will be stored in the same forlder with all the data 
#--------------------------------------------------------------

# MIN and MAX are values to remove the extreme of your spectra (aberrations on the edge) if None will be used the whole spectra
MIN=MIN_WAVE_CCF
MAX=MAX_WAVE_CCF

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

# Variables given on INPUT.py

# READ FILES HERE IF YOU WANT TO USE THIS MANUALLY

# OBJ_Template='Synth_Spectra'
# Template_wave=syn_spec[:,0]
# Template_flux=syn_spec[:,1]

Template_wave=WAVE_TEMP
Template_flux=FLUX_TEMP


#--------------------------------------------------------------

# Useful Functions

# Shift wavelenght according to RV (Doppler Shift)
def corr_rv(wave,rv):
    return wave/(1+rv/299792.458)

# Function that reads fits Output of Merging ans substracting sky 
def get_w_f(path):
    with fits.open(path) as hdul:
        OBJ  = hdul[0].header['OBJECT']
        RVC  = hdul[0].header['RVBARY']
        RA   = hdul[0].header['RA']
        DEC  = hdul[0].header['DEC']
        MAG  = float(hdul[0].header['MAG'])
        wave = hdul[1].data['wave']

        if (MIN_WAVE_CCF is None) or (MAX_WAVE_CCF is None): # If values are not given we trim 1/20 per extreme
            Ext    =[min(wave),max(wave)]
            delta  =(Ext[1]-Ext[0])/20
            MIN,MAX=Ext[0]+delta,Ext[1]-delta

        mask = (MIN<wave)&(wave<MAX)     # Cut the edges os spectra
        wave = wave[mask]
        flux = hdul[1].data['flux'][mask]

        MOflux= np.median(flux)

        if not (0.5<MOflux<1.5):
            flux/=MOflux  # Scaled to 1
    
    return wave,flux,OBJ,MAG,RVC,RA,DEC

# Set frame of figure for Summary plot
def MYFigure():
    Fig= plt.figure(777,figsize=(13,6))
    gs = gridspec.GridSpec(
        3,2,                    # Number of axis y,x
        height_ratios = [1,1,1],# relatives ratios of heigh
        width_ratios  = [1,1],  # relatives ratio of with
        left  = 0.1,            # Space to the edge of left   from the nearest axis
        right = 0.95,           # Space to the edge of right  from the nearest axis
        bottom= 0.08,           # Space to the edge of bottom from the nearest axis
        top   = 0.95,           # Space to the edge of top    from the nearest axis
        wspace= 0.2,            # Space horizontal between each of the axis
        hspace= 0.2)            # Space vertical   between each of the axis

    ax0=Fig.add_subplot(gs[:,0])
    ax1=Fig.add_subplot(gs[0,1])
    ax2=Fig.add_subplot(gs[1,1])
    ax3=Fig.add_subplot(gs[2,1])

    return ax0,ax1,ax2,ax3
    

# Function to add the calculated RV to the original fits file
def UPDATE_header(path,RV,RVERR,RVTRUE):
    with fits.open(path) as F:
        F[0].header['RV']    = (RV,'Radial Velocity Obtained with CCF')
        F[0].header['RVERR'] = (RVERR,'Estimated Error in RV measure on CCF')
        F[0].header['RVTRUE']= (RVTRUE,'True Radial Velocity (Measured RV+RVBARY)')
        F.writeto(path, overwrite=True)

# Function that creates a New Fits file with the flux on restframe
def new_fits(path,name,RV):
    with fits.open(path+name) as F:
        header=F[0].header
        try:
            wave  =F[1].data.wave
            flux  =F[1].data.flux
        except:
            flux =hdul[0].data
            CRVAL=hdul[0].header['CRVAL1']
            CDELT=hdul[0].header['CDELT1']
            wave =CRVAL+np.arange(len(flux))*CDELT

        try: # No longer 1D file since wavelength is not evenly spaced
            header.pop('CRVAL1')
            header.pop('CDELT1')
        except:pass

        phdu=fits.PrimaryHDU(header=header)

        # Create the columns of data
        columns=[
            fits.Column(name='wave', format='D', unit='Angstrom',array= corr_rv(wave,RV) ),
            fits.Column(name='flux',       format='D', unit='adu',     array= flux )
            ]
        # Generate the secondary fits with the data
        NEWFITS = fits.BinTableHDU.from_columns(columns)

        # Set the EXTNAME keyword to label the binary table HDU
        NEWFITS.header['EXTNAME'] = 'SPECTRUM (restframe)'
        # Comments about the columns
        NEWFITS.header['TCOMM1']  = 'wavelenght'
        NEWFITS.header['TCOMM2']  = 'Flux' 
        # Save the fits file
        hdul = fits.HDUList([phdu, NEWFITS])
        hdul.writeto(path+'RESTFRAME/'+name, overwrite=True)

# Function that stick all function togethers (Calls the CCF)
# Also if you decided to plot each one of the plots the info is defined here
def get_rv(path):
    wave ,flux ,OBJ ,MAG ,RVC ,RA ,DEC =get_w_f(path)

    myCCF=CCF(template_wavelength   = Template_wave,  
                template_flux       = Template_flux, 
                observed_wavelength = wave, 
                observed_flux       = flux,
                IGNORE              = False,
                PRINT               = False,
                RV_range            = RADVEL_RANGE,
                )

    RV_sample,Corr,RV,ERR=myCCF.find_rv(steps3=500)

    print(f'{OBJ:<17} {RV-float(RVC):9.3f}\t {ERR:5.3f}\t {RVC:9.3f}\t {RV:9.3f}\trelative_to_{OBJ_Template}')

    if PLOT_EACH_ONE and not PLOT_sumary:
        myCCF.PLOT(PAUSE=True,time=1)

    return RV_sample,Corr,RV,ERR, wave ,flux ,OBJ ,MAG ,RVC ,RA ,DEC

#--------------------------------------------------------------


if SAVE:
    try:
         os.makedirs(path_files+'RESTFRAME', exist_ok=True)
    except:
        pass
    my_file=open(path_files+SAVE_NAME,'w')
    my_file.write(f'# Method used {method}\n')
    my_file.write(f'# Path to sythetic spectra: {path_syn_spec}\n')
    

if PLOT_sumary:
    ax0,ax1,ax2,ax3= MYFigure()



###############################################################################################
# Checking method used
#-----------------------
if method=='Synthetic':
    pass
#-----------------------
elif method=='Template':# Use template for one star , an then that star to match the others
    print('\n#--------------------------------------------------------\n')
    print('# Template')
    RV_sample,Corr,RV,ERR, wave ,flux ,OBJ ,MAG ,RVC ,RA ,DEC= get_rv(path_temp_star)

   
    if SAVE:
        my_file.write(f'# Path to Template star: {path_temp_star}\n')
        my_file.write(f'# RV of Template star: {RV} [km/s]\n')

    OBJ_Template ='Temp_Star_'+str(OBJ)
    Template_wave= corr_rv(wave,RV)
    Template_flux= flux.copy()
    
    MIN=MIN+10
    MAX=MAX-10
    print('\n#--------------------------------------------------------\n')
#-----------------------
else:
    raise Exception('Method inputed is not recognised, please set method to "Synthe" or "Template" ')
#-----------------------
###############################################################################################


if SAVE:
    my_file.write('# Radial velocities are in km/s \n')
    my_file.write(f'#Object\tRA\tDEC\tRV\tRV_ERR\tBarycentric_RV\tRV_corrected\ttemplate_used_in_CCF\n')


print(f'{"OBJ":<17} {"RV":<9}\t {"RV_ERR":<5}\t {"BarycenRV":<9}\t {"RV_corr":<9}\trelative_to_What?')


for File in list_file:
    # Fits file that we will read
    fits_file= path_files+File

    # Get relevant info
    RV_sample,Corr,RV,ERR, wave ,flux ,OBJ ,MAG ,RVC ,RA ,DEC= get_rv(fits_file)
    
    # Save in a file relevant info
    if SAVE: 
        # On txt 
        my_file.write(f'{OBJ}\t{RA}\t{DEC}\t{RV-float(RVC)}\t{ERR}\t{RVC}\t{RV}\trelative_to_{OBJ_Template}\n')
        # On original fits
        UPDATE_header(fits_file,RV-float(RVC),ERR,RV)
        # And we create a new Fits
        new_fits(path_files,File,RV)
    
    # Show Summary plot
    if PLOT_sumary: 
        ax0.plot(RV_sample,Corr)
        new_wave=corr_rv(wave,RV)
        ax1.plot(new_wave,flux)   
        ax2.plot(new_wave,flux)
        ax3.plot(new_wave,flux)  
        
  

if PLOT_sumary:

    MTflux= np.median(Template_flux)
    
    if not (0.25<MTflux<1.5):
        Template_flux/=MTflux
    
    # just reference points just in case it was not defined
    Npoints   = len(new_wave)
    Npoint_ref= (int(Npoints*0.2),int(Npoints*0.55),int(Npoints*0.6),int(Npoints*0.725))
    
    try:
        if None in WINDOW_1:
            WINDOW_1=(new_wave[Npoint_ref[0]],new_wave[Npoint_ref[1]])
    except:
        WINDOW_1=(new_wave[Npoint_ref[0]],new_wave[Npoint_ref[1]])
    else:
        print('Window 1 was not identified, it will show just a generic window')
        WINDOW_1=(new_wave[Npoint_ref[0]],new_wave[Npoint_ref[1]])
    
    try:
        if None in WINDOW_2:
            WINDOW_2=(new_wave[Npoint_ref[2]],new_wave[Npoint_ref[3]])
    except:
        WINDOW_2=(new_wave[Npoint_ref[2]],new_wave[Npoint_ref[3]])
    else:
        print('Window 2 was not identified, it will show just a generic window')
        WINDOW_2=(new_wave[Npoint_ref[2]],new_wave[Npoint_ref[3]])
    

    print('\nTemplate spectra diplayed on black\n')
    ax1.plot(Template_wave,Template_flux,color='black',lw=2)
    ax1.set_ylim(0,np.nanpercentile(flux,86)*1.35)


    ax2.plot(Template_wave,Template_flux,color='black',lw=2)
    ax2.set_ylim(0,np.nanpercentile(flux,86)*1.35)
    ax2.set_xlim(WINDOW_1[0],WINDOW_1[1])

    ax3.plot(Template_wave,Template_flux,color='black',lw=2)
    ax3.set_ylim(0,np.nanpercentile(flux,86)*1.35)
    ax3.set_xlim(WINDOW_2[0],WINDOW_2[1])

    ax0.set_ylabel('Correlation')
    ax1.set_xlabel('RV [km/s]')
    ax3.set_xlabel('Wavelength')

    if SAVE_SUMMARY_PLOT:
        plt.savefig(path_files+'CCFS.png',dpi=100)
        plt.close(777)
    else:
        plt.show()
        plt.close(777)

if SAVE:
    my_file.close()