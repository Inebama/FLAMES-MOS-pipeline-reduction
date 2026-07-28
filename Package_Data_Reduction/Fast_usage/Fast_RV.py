# Add path where is the location of the files
import sys
path_package=__file__.replace('Fast_usage/Fast_RV.py','')
sys.path.append(path_package)

# Easy way to have a CCF for a single/multiple stars
from CrossCorrelation import CCF
from astropy.io import fits
from User_Functions import TEMPLATE_READER
import numpy as np
from glob import glob

# Path to stars
path_star=path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_data_ian/OB-0/'
path_stars=glob(path_star+'*.fits')[3:5]
print(path_star)



# Get a template
complete_path_template:str =path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR14b_R37000.dat'
Template_wave,Template_flux=TEMPLATE_READER(complete_path_template)

RADVEL_RANGE:tuple[float]=(-120,120)


# Modify as your needs (this is for tipycal output of the software)
def get_w_f(path:str,             # Path to the fits file
            MIN :float =6391.4,   # Trim the wavelength range to cut edges were aberrations may appear
            MAX :float =6613,
            ):
    with fits.open(path) as hdul:
        OBJ  = hdul[0].header['OBJECT']
        RVC  = hdul[0].header['RVBARY']
        RA   = hdul[0].header['RA']
        DEC  = hdul[0].header['DEC']
        MAG  = float(hdul[0].header['MAG'])
        wave = hdul[1].data['wave']

        mask = (MIN<wave)&(wave<MAX)     # Cut the edges of spectra
        wave = wave[mask]
        flux = hdul[1].data['flux'][mask]

        # MOflux= np.median(flux)

        # if not (0.5<MOflux<1.5):
        #     flux/=MOflux  # Scaled to 1
    
    return wave,flux,OBJ,MAG,RVC,RA,DEC

# General way to use the CCF
def get_rv(path):
    wave ,flux ,OBJ ,MAG ,RVC ,RA ,DEC =get_w_f(path)

    myCCF=CCF(template_wavelength   = Template_wave,  
                template_flux       = Template_flux, 
                observed_wavelength = wave, 
                observed_flux       = flux,
                IGNORE              = False,
                PRINT               = True,
                RV_range            = RADVEL_RANGE,
                )

    RV_sample,Corr,RV,ERR=myCCF.find_rv(steps3=500)

    print(f'{OBJ:<17} {RV+float(RVC):9.3f}\t {ERR:5.3f}\t {RVC:9.3f}\t {RV:9.3f}')
    myCCF.get_errs_info(PRINT=True)
    myCCF.PLOT()

#-------------------------------------

for path in path_stars:
    get_rv(path)

