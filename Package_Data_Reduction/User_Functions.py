#-----------------------------
# Packages
import numpy as np 
from astropy.io import fits
import glob
import os
#-----------------------------

# Keys that you have to provide in a function if you create one for your specific data
# Rules associated with the data are informed after the comments
KEYS_FOR_TABLES_READER_FUNCTIONS=['MJD',   # DATE of observation must be MODIFIED JULIAN DATE (MJD), if you dont have it use astropy utilities to convert to MJD
                                  'RAS',    # Right Ascention associated with the specific measure (DATE) len(RAS)==len(DATE)
                                  'DECS',   # Declinations    associated with the specific measure (DATE) len(RAS[n])==len(DECS[n])
                                  'GLAT',   # Coord: Latitude of the telescope
                                  'GLON',   # Coord: Longitud of the telescope
                                  'HTEL',   # Coord: Eleveation above sea of the telescope
                                  'WAVE',   # Wavelength associated with each one of the positions (RA, DEC) len(RAS [n]==len(WAVE[n]))
                                  'FLUX',   # Flux       associated with each one of the positions (RA, DEC) len(DECS[n]==len(FLUX[n]))
                                  'NOT_SKY', # (Optional) If among your data you have sky fibers or files this will be a way to provide that info (must be strictly associated with len(RAS[n]==len(DECS[n])==len(NOT_SKY[n]) and len(NOT_SKY)==len(DATE)
                                  'OBJS',    # (Optional) name of each of the sources must be len(RAS[n]==len(DECS[n])==len(OBJS[n]))
                                  'SHAPE',   # (Optional) refers the shape of the data this will be printed if given
                                  'MAGS',    # (Optional) refers to the magnitude of the targets, will be used to plot results otherwise generic values will be used instead
                                  ] 

# You must convert the flux to a np.array in order to properly work

#____________________________________________________________________________________________________

# Function that given a epoch (folder) will look retrieve relevant info for listing and merging files (this was done accordingly to the output given by ESOREFLEX)
def TABLES_UVES_MOS_REDUCED(epoch,ARM='U'): # ARM could be  'U' or 'L'

    DATA_FROM_FILE= dict()

    table = glob.glob(epoch+f'*FIB_SCI_INFO_TAB_RED{ARM}.fits')[0]   # Look for the bintable with the info

    with fits.open(table) as hdul:
        mjd    = hdul[0].header['MJD-OBS']
        types  = hdul[1].data.TYPE
        mask   = (types!='')  # The pipeline generates some lines empty
        types  = types[mask]
        RA     = hdul[1].data.RA[mask]
        DEC    = hdul[1].data.DEC[mask]
        OBJ    = hdul[1].data.OBJECT[mask]

        # Is a sky fiber? True or False
        not_sky= (types!='s')&(types!='S') # Indentify the sky fibers in the table
        DATA_FROM_FILE['NOT_SKY']=not_sky  # Corresponding to the Files containing actual data, which are not sky fibers
        
        DATA_FROM_FILE['RAS']  = RA  # Positions on Rigth Ascention
        DATA_FROM_FILE['DECS'] = DEC # Positions on Declination
        DATA_FROM_FILE['OBJS'] = OBJ
        
        # Data of Observation
        DATA_FROM_FILE['MJD']= mjd # Dates of each set of data
        DATA_FROM_FILE['GLAT']= hdul[0].header['HIERARCH ESO TEL GEOLAT' ]
        DATA_FROM_FILE['GLON']= hdul[0].header['HIERARCH ESO TEL GEOLON' ]
        DATA_FROM_FILE['HTEL']= hdul[0].header['HIERARCH ESO TEL GEOELEV']

    # Inside of each folder of epoch we look the files that we want 
    files = list(filter(lambda x:f'_MWXB_SCI_RED{ARM}' in x and not 'ERR' in x,glob.glob(epoch+'*.fits')))
    files.sort()

    DATA_FROM_FILE['SHAPE']= len(files) # This is optional (will be used to print alongside in the list files)

    FLUXES=[]
    WAVES =[]
    for file in files:
        with fits.open(file) as F:
            # Store Data
            Flux=F[0].data       

            CRVAL=F[0].header['CRVAL1']
            CDELT=F[0].header['CDELT1']
            N_pix=F[0].header['NAXIS1']

            # # Print info just to check that we are talking about the same kind of files
            # print(file,Flux.shape,end='\t')
            # print(round(CRVAL,3),round(CDELT,3),N_pix)

            wave = CRVAL+np.arange(N_pix)*CDELT

            FLUXES.append(Flux)
            WAVES.append(wave)
         
    DATA_FROM_FILE['FLUX'] = np.array(FLUXES) # Saving Flux data # All have the same shape for a single exposure

    DATA_FROM_FILE['WAVE'] = wave #np.array(WAVES) # All WAVES are equal for a single exposure
    
    return DATA_FROM_FILE
    # Just the first 5 are strictly required but the other parameters will be used to complete headers if not given just generic info will be added

#____________________________________________________________________________________________________

# Function that given list of files and return the info required for Giraffe in a dictionary with variables names that will be need
# In order to properly work define the same key names of variables for the dictionary that is output
def TABLES_GIRAFFE_MOS_REDUCED(file):
    # We will store the data in a dictionary
    DATA_FROM_FILE= dict()
    
    with fits.open(file) as F: 
        # EXPT=F[0].header['EXPT']
        # GAIN=F[0].header['HIERARCH ESO DET OUT1 GAIN']
        # Flux=GAIN*F[0].data/EXPT 
        Flux=F[0].data      
        DATA_FROM_FILE['FLUX'] = Flux # Saving Flux data in e-/s

        CRVAL=F[0].header['CRVAL2'] # Number 2 is due it refers to the second axis
        CDELT=F[0].header['CDELT2']
        N_pix=F[0].header['NAXIS2']

        DATA_FROM_FILE['WAVE'] = CRVAL+np.arange(N_pix)*CDELT

        # Print info just to check that we are talking about the same kind of files
        # print(file,Flux.shape,end='\t')
        # print(round(CRVAL,3),round(CDELT,3),N_pix)


        # Store relevant info
        DATA_FROM_FILE['OBJS'] = F[1].data['OBJECT'   ]
        DATA_FROM_FILE['RAS']  = F[1].data['RA'       ]
        DATA_FROM_FILE['DECS'] = F[1].data['DEC'      ]
        DATA_FROM_FILE['MAGS'] = F[1].data['MAGNITUDE'] # This is optional (will be used to plot with corresponding colors)
        DATA_FROM_FILE['SHAPE']= Flux.shape             # This is optional (will be used to print alongside in the list files)
        DATA_FROM_FILE['FIBER']= F[1].data['RP']        # This is optional (may help in to track any issue with specific fibers)
        # Is a sky fiber? True or False
        TYPE= F[1].data['TYPE'] 
        DATA_FROM_FILE['NOT_SKY'] = (TYPE!='S')&(TYPE!='s') # This is Optional, in case of not being provided it will assume that all your files are science targets

        # Data of Observation
        DATA_FROM_FILE['MJD'] = F[0].header['MJD-OBS'                 ]
        DATA_FROM_FILE['GLAT']= F[0].header['HIERARCH ESO TEL GEOLAT' ]
        DATA_FROM_FILE['GLON']= F[0].header['HIERARCH ESO TEL GEOLON' ]
        DATA_FROM_FILE['HTEL']= F[0].header['HIERARCH ESO TEL GEOELEV']
        # FIBN=F[1].data.FPD

    return DATA_FROM_FILE
    #-----------------------------------------------------------


#____________________________________________________________________________________________________
# EXAMPLE function for continuum fititng
# giving an hdu list (that will be read with astropy)
def READ_FITS_CONTINUUM(hdul):
    try:# Multidimensional data format
        wave=hdul[1].data['wave']
        flux=hdul[1].data['flux']
    except:
            try:
                # We scale the values of teh spectra around one for easier computation
                wave=hdul[1].data['wavelength']
                flux=hdul[1].data['Flux']
            except:# 1D format standard IRAF
                flux =hdul[0].data
                CRVAL=hdul[0].header['CRVAL1']
                CDELT=hdul[0].header['CDELT1']
                wave =CRVAL+np.arange(len(flux))*CDELT
    return wave,flux

#____________________________________________________________________________________________________
# EXAMPLE function for reading syntehtic spectra
# works for outputs of AUTOKUR spectra

def TEMPLATE_READER(path):

    syn_spec=np.genfromtxt(path,skip_header=47)
    Template_wave=syn_spec[:,0]
    Template_flux=syn_spec[:,1]
    # Template_flux/=np.median(Template_flux)

    return Template_wave,Template_flux
