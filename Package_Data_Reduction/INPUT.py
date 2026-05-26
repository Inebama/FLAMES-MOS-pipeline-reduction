#############################################################################
# Instructions 
#
# Here you have to set the values in these specific variables
# DO NOT CHANGE the name of the variables or functions, just edit their content  
# Make sure that the outputs of the function deliver all the information required 
#
# Run this code and see the outputs untill you are sure that the functions are properly printing the required information
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
import numpy as np 
from User_Functions import TABLES_UVES_MOS_REDUCED
from User_Functions import TABLES_GIRAFFE_MOS_REDUCED
from User_Functions import TEMPLATE_READER
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
# General 

# This line will read were was placed this file, for personalized runs prefer to give full paths 
path_package=__file__.replace('INPUT.py','')

# Paths (make sure that you are using full paths and ends in '/')

## Testing software
# path_input :str =path_package+'TESTING_FOLDER/DATA_GIRAFFE/'
# path_output:str =path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_F/'
path_input :str ='/scratch/home/ibaeza/Desktop/DATA/GIRAFFE_MOSSPECTRA/GIRAFFE/'
path_output:str ='/scratch/home/ibaeza/Desktop/MY_OUTPUT/GIRAFFE_median_24h_corrected/'


## Maca DATA test (all data treated separately)
# path_input :str=path_package+'TESTING_FOLDER/DATA_MACA_GIRAFFE_6273/'
# path_output:str=path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_2/'

## After first experiment I did splitted the data into 2 groups and I have merged them
# path_input :str=path_package+'TESTING_FOLDER/DATA_MACA_GIRAFFE_6273/BLOCK1/'
# path_output:str=path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_Split1/'
# path_input :str=path_package+'TESTING_FOLDER/DATA_MACA_GIRAFFE_6273/BLOCK2/'
# path_output:str=path_package+'TESTING_FOLDER/TEST_OUTPUT/Test_Split2/'

# Here you can explicitly enter what task do you want to perform 
# Also you can give this information whe running on command line (see MAIN.py for details)

LIST      :bool = True  # Looks for files that could be treated as the same Observing block
MERGE     :bool = False  # Make the sky substraction and merge of files
RADIAL_VEL:bool = False  # Make a cross correlation to determine Radial velocities and store them
CONTINUUM :bool = False  # Start a continuum fitting for files and store them 


###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________


# Specific methods setup

###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________

# Checking Observing Blocks

# This is asking to you if the files that must list is in files or folders

FITS_FILE:bool = True  # implies that must search for *.fit* files
#FITS_FILE:bool = False # implies that will look for */ 

# (Specific search can be modified below in fine adjustments)  

Range_in_hours_to_be_consider_same_Observing_Block :float =24 # 12 hours imply basically one night, and less than 0.25 likely will just separate individual exposures

# Stablish this function, this should give the information associated to the fits files considering an input complete path to a single file 

# read_table=TABLES_UVES_MOS_REDUCED
read_table=TABLES_GIRAFFE_MOS_REDUCED

# Do you want in addition to store the OB block in a human readable txt?
SAVE_TXT= True

###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________

# Merging and Substracting sky

# Choose function to merge
function_used_to_merge= np.nanmedian       # Strongest against outliers
# function_used_to_merge= np.nanmean         # Converges faster to "true" values


# Range to display in figures to check good merge
# you should check that are in the units of your data

# GIRAFFE HR14b
# min_wave: float|None  = 655
# max_wave: float|None  = 657.5 

# UVES 580 L
# min_wave : float|None = 5150
# max_wave : float|None = 5300 

# UVES 580 U
# min_wave : float|None = 6550
# max_wave : float|None = 6575  

# GENERIC
min_wave: float|None  = None
max_wave: float|None  = None  

# Factor to convert your wavelength to angstrom (you must know this from the header of your data)

wave_to_angstrom: float=10   # For GIRAFFE (pipeline retrieves wavelenght in Nanometers)
# wave_to_angstrom: float=1  # For UVES    (pipeline retrieves wavelenght in Angstrom)


###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________

# Radial Velocity 

# The semiautomatic mode will look in the output folder for files to perform the CCF
# Depending on how you have performed the tasks (first getting RV or fit the continuum) you will have to provid ea subfolder were the actual files that you want to adjust are placed 

EXTRA_FOLDER_RV  :str = ''               # Letting empty means that will be looking directly in the output folder
# EXTRA_FOLDER_RV:str = 'Continuum_sub/' # Continuum sub is the folder that it is created in the output directory when running the Continuum


# Path to synthetic spectra or template star to perform CCF (Just for the Radial Velocity)

complete_path_template:str =path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR14b_R37000.dat'
# complete_path_template:str =path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/GIRAFFE_TEMPLATE_HR13_R27000.dat'
# complete_path_template:str =path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/UVES_L_ARM_580_R47000.dat'
# complete_path_template:str =path_package+'TESTING_FOLDER/SYNTH_TEMPLATE/UVES_R_ARM_580_R47000.dat'


# Function that given the path reads the template
reads_template=TEMPLATE_READER


# Radial velocity Range that you want to explore in the CCF
RADVEL_RANGE:tuple[float]=(-250,250)


# MIN and MAX are values to remove the extreme of your spectra (aberrations on the edge) if None will be used the whole spectra
# Generic
MIN_WAVE_CCF: float =None
MAX_WAVE_CCF: float =None

# GIRAFFE HR14b
# MIN_WAVE_CCF: float =6391.4
# MAX_WAVE_CCF: float =6613

# UVES 580 L
# MIN_WAVE_CCF: float =4800
# MAX_WAVE_CCF: float =5750

# UVES 580 U
# MIN_WAVE_CCF: float =5850
# MAX_WAVE_CCF: float =6800


# Name to identify to respect what we will be doing the CCF 
OBJ_Template:str ='Synth_Spectra' # This will appear on the output file


# Alternatively you could just use the following variables and directly gave the template flux and wavelength reading the template on this file
WAVE_TEMP: float|None =None
FLUX_TEMP: float|None =None


# Do you want to see each one of the fits done with the CCF (for all the spectra given)
PLOT_EACH_ONE:bool= False # If False either way you will get a summary plot to check that everything is working well

# For summary plot (PLOT_EACH_ONE= False)
SAVE_SUMMARY_PLOT:bool= True

# Windows that will be displayed 

# WINDOW_1: tuple =None
# WINDOW_2: tuple =None

# HR14b
# WINDOW_1: tuple =(6488, 6500.5)
# WINDOW_2: tuple =(6552, 6576  )

# UVES L arm
# WINDOW_1: tuple =(5150, 5300)
# WINDOW_2: tuple =(5400, 5450)

# UVES U arm
# WINDOW_1: tuple =(5890, 5905)
# WINDOW_2: tuple =(6552, 6576)

# Compute all RV for each observation at once
Across_OB: bool = True

###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________

# Continuum Fitting 

# The semiautomatic mode will look for the files on the output directory and for the files there will perform the continuum fitting
# Depending in the order that you want ot perfomr the task (first fit continuum and then get RV) you will have to rpovide the path top the subfolder created
EXTRA_FOLDER_CONT : str = 'RESTFRAME/' # Restframe is the folder that it is created in th eoutput directory when running the radial velocity
# EXTRA_FOLDER_CONT:str = ''           # Letting empty means that will be looking directly in the output folder


# WARNING: For now the step of continuum fiting just works with fits files
# Because of how is constructed this package if you want to edit this function it should take as an input the HDU-list and from that point create the specific way of read
# you can check a an example in the file of User_Functions.py called READ_FITS_CONTINUUM (this is the default and is alredy inegrated)

read_fits=None
# read_fits=None # Setting this variable to None will use a generic reader that is compatible with previous outputs of the code and iraf like fits

# Range of wavelength that you want to asjust the Continuum (usually extremes values in the CCD are messed so this is an easy way to get rid of those extremes)
# Also you can just simply set a value between 0 and 1 which will represent the percentage of the lenght of your spectra that you just want to get rid of the edges of a lot of different spectrums
WAVE_RANGE: float|tuple[float] =0.04 
# WAVE_RANGE: float|tuple[float] =(6391.4,6613)


# If you want to automatically fit a continuum you can do it by setting this keyword (however we highly recommend to do it intercatively)
# While this creates a good approaximation be cautios of the edges of your spectra and in the vicinity of strong lines (this last point might be more mitigated using larger shifts)
AUTOMATED_CONT:bool=False

###################################################################################################################
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------

# Fine adjustments for semi-automatic mode

# Which kind of files should Look for?
key: str ='*.fit*' if FITS_FILE else '*/'

# When merging fluxes on which axis we should use?
AXIS_TO_MERGE:int=0



#############################################################################
#############################################################################
# HERE AFTER NOT TOUCH
#############################################################################
#############################################################################

# Do not change this, The code below will allow you to test that you have anything that you need if running this script in python does not arise you an error you should probably be fine
# Any way if somenthign happens I'm providing after the import the version of the packages that I have used to do these codes
if __name__=='__main__': # This tell to python that this piece of code only will be run if you are directly compiling this code and not when importing it
    # All Packages that we use on this software  
    #import numpy as np                            # numpy      version '1.26.4'
    import matplotlib as mpl                       # matplotlib version '3.9.0'
    import matplotlib.pyplot as plt                
    import matplotlib.gridspec as gridspec
    from matplotlib import lines 
    from matplotlib.widgets import Slider, Button
    from scipy.optimize import curve_fit
    from scipy.interpolate import interp1d  
    from scipy.interpolate import UnivariateSpline # scipy      version '1.11.4'
    from astropy import units as u                 # astropy    version '5.3.4'
    from astropy.io import fits                    
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, Angle

    #----------------------------------
    from prettytable import PrettyTable           # prettytable version '3.10.0' (pip install prettytable)
    #----------------------------------

    import os
    import sys
    import glob 
    # import typing #Deprecated
    try:
        from tkinter import filedialog as fd # Will ask for folders 
    except:
        pass
    #----------------------------------
    # Start Testing

    # We will test if the functions provided can be used to read the files
    if LIST or MERGE:
        files = glob.glob(path_input+key)
        files.sort()
        OUTPUT= read_table(files[0])
        print('\n\n mean RA\t mean DEC\t DATE\t How Many Skyes?\t First_File_Found\n')
        print(np.mean(OUTPUT['RAS']),np.mean(OUTPUT['DECS']),OUTPUT['MJD'],sum(~OUTPUT['NOT_SKY']),files[0])
        print('\nIf this match what you expect (check please) the LIST and MERGE method should work fine')

    if MERGE:
        FLUX=OUTPUT['FLUX']
        shape1=FLUX[0].shape
        if len(shape1)>1:
            if shape1[0]>shape1[1] and AXIS_TO_MERGE!=0:
                sys.exit('WARNING:; I think that AXIS_TO_MERGE must be 1 for your data')
            elif shape1[1]>shape1[0] and AXIS_TO_MERGE!=1:
                sys.exit('WARNING:; I think that AXIS_TO_MERGE must be 0 for your data')
        elif AXIS_TO_MERGE!=0:
            sys.exit('WARNING:; I think that AXIS_TO_MERGE must be 0 for your data')


    # we will test that the reader or the variables are given 
    if RADIAL_VEL:
        if not WAVE_TEMP is None and not FLUX_TEMP is None and len(WAVE_TEMP)==len(FLUX_TEMP):
            print('\nYou have provided directly the template WAVE and FLUX, so I will asume that you know what you are doing \n')
        else:
            temp_wave,temp_flux=TEMPLATE_READER(complete_path_template)
            if len(temp_wave)!=len(temp_flux):
                sys.exit(' WARNING:; While checking the template wave and flux I have found that they missmatch please check your function of template reader')
    
    
    # We will check that the data could be read with the function provided
    if CONTINUUM:
        if not read_fits is None:
            with fits.open(path_output) as hdul:
                wave,flux=read_fits(hdul)
            if len(wave)!=len(flux):
                sys.exit(' WARNING:; While checking the continuum function I have found that the flux and wavelength does not match please check this')
    

    
    print(f'\n\nDONE tests, as long as you have checked the printed data above, for the methods selected at the beggining of this file:\n\n   LIST      :{LIST}\n\n   MERGE     :{MERGE}\n\n   RADIAL_VEL:{RADIAL_VEL}\n\n   CONTINUUM :{CONTINUUM}\n\n everything should work well\n')
    
    THIS_SHOULD_NOT_BE_DEFINED_ON_RUNNING_MAIN_FILE=170499