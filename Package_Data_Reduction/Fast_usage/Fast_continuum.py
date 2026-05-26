# Add path where is the location of the files
import sys
path_package=__file__.replace('Fast_usage/Fast_continuum.py','')
sys.path.append(path_package)

# Packages
from astropy.io import fits
from Continuum_fitting import cont_fit
from glob import glob
import os


# Path where your spectras are stored
PATH=path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_data_ian/OB-0/RESTFRAME/'

# Name of your files in a list
# Files=['MED_RGB_107']

# Just a lazy way to put all my files started with MED for the example case
Files=list(filter(lambda x:'MED' in x,glob(PATH+'*')))
Files=list(map(lambda x: x.split('/')[-1],Files))


# Range of wavelength that you want to asjust the Continuum (usually extremes values in the CCD are messed so this is an easy way to get rid of those extremes)
# Also you can just simply set a value between 0 and 1 which will represent the percentage of the lenght of your spectra that you just want to get rid of the edges of a lot of different spectrums
# WAVE_RANGE: float|tuple[float] =0.04 
WAVE_RANGE: float|tuple[float] =(6391.4,6613) # HR14b


# Provide a Function that can retrieve your flux and wavelength from a hdul
# read_fits=None
def read_fits(hdul):
    try:
        wave=hdul[-1].data.wave
        flux=hdul[-1].data.flux
    except:
        wave=hdul[-2].data.wave
        flux=hdul[-2].data.flux

    return wave,flux


################################################################

CONT=cont_fit(Files,                   # a list with the files that we will work with (not full path)
              PATH,                    # The path where these files exist 
              wave_range=WAVE_RANGE,   # Wavelength range that will be used to perform the fit
              Continue=True,           # Look for the file that contains previous adjustments in the folder
              READ_FUNCTION=read_fits  # Function used to read files (just fits files) if None it will use defaul values
            )

# Once the complete "OB" is done it will save the corrected by continuum in a subfolder
CONT.store_cont_sub(Force=False,) 