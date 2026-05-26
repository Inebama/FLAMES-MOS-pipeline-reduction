#----------------------------------------------------------------
#Importing the inputs defined in the input file
from INPUT import *
import numpy as np

#----------------------------------------------------------------
### To run options from command line Uncomment the following block of code

# import sys

# # you can customize how do you want to gave this input
# OPTIONS=['-list','-merge','-continuum','-ravel']
# # e.g
# # python MAIN.py -list

# # If some of the methods were given as a command line they would have priority 
# # Otherwise we will use the values defined in INPUT file
# try: 
#     LIST= OPTIONS[0] in sys.argv
# except:
#     LIST
# try: 
#     MERGE= OPTIONS[1] in sys.argv
# except:
#     MERGE
# try: 
#     CONTINUUM= OPTIONS[2] in sys.argv
# except:
#     CONTINUUM
# try: 
#     Radial_Vel= OPTIONS[3] in sys.argv
# except:
#     Radial_Vel


#----------------------------------------------------------------
# To Handle the output of the functions that reads tables
def merge_dicts(dicts:list[dict]) -> dict:
    merged    = {}
    keys_order= ['MJD','RAS','DECS','NOT_SKY','GLAT','GLON','HTEL','WAVE','FLUX','OBJS','SHAPE','MAGS']
    # Get all unique keys from the dictionaries
    all_keys = set(key for d in dicts for key in d.keys())
    for key in keys_order:
        # Collect values from all dictionaries for the current key
        if key in all_keys:
            data        = [d.get(key) for d in dicts]
            array       = np.empty(len(data), dtype=object) # Must be done in this way to avoid numpy trying to reshape array to have the same size
            array[:]    = data
            merged[key] = array
        else:
            merged[key] = None  # Use None for keys that are not in all the dictionaries
    
    return merged


#----------------------------------------------------------------



# Here we will execute the OBSERVING BLOCK code 
# You can use this file as a reference of how this code must be used manually

if LIST or MERGE:
    
    import glob # To look files (like an ls in unix)
    import os   # To make sure that the output directory exist or create it in case of neccesary
    
    # We look for all the files with the key defined in INPUT in the input directory 
    Files=[i.replace(path_input,'') for i in glob.glob(path_input+key)]
    Files.sort() # Ordered by name

    print('\nReading data ...\n')
    # We obtain the relevant info of the files using the USER FUNCTIONS defined for you specific kind of File
    DATA=[read_table(path_input+file) for file in Files]

    # We unify the informations of the different files in an order way (each file data will correspond to an index of all variables in the dicitonary)
    DATA=merge_dicts(DATA)

    print('Ready\n\nChecking Output Directory ...\n')
    # We create the output file if does not exist already
    os.makedirs(path_output,exist_ok=True)

    if LIST:
        from List_files_related import OBSERVING_BLOCK
        print('Ready\n\nGetting Observing Blocks ...\n')

#____________________________________________________________________

        # We create the "observing blocks"
        OBSERVING_BLOCK(OBS_BLOCK_TIME=Range_in_hours_to_be_consider_same_Observing_Block, # Range in hours that you want to consider as the same observing block (12 would work for separate by night) 
                        path   = path_output,        # path where the file with the info of the OB will be stored 
                        epochs = Files,              # Name of files/folders for easy identification of data that is associated with the OB
                        time   = DATA['MJD'],        # time in mJD or similar (must be in days)
                        RA     = DATA['RAS' ],       # List of Right Ascension associated with the specific times of measure (Variable time)
                        DEC    = DATA['DECS'],       # List of Declination     associated with the specific times of measure (Variable time)
                        not_sky= DATA.get('NOT_SKY'),# A List with a mask of the size of RA and DEC (len(RA)==len(DEC)==len(not_sky) and ) / if set just true imply that you don't have any sky in the given
                        shapes = DATA.get('SHAPE'),  # (Optional) A variable with the shapes of the data (might be useful to check)
                        SAVE_TXT= SAVE_TXT,          # Do you want to store data in human readable txt?
                        )
        
#____________________________________________________________________
    if MERGE:
        # Here I will asume that the observation blocks already have been created
        from Merging_and_substract_sky import MERGING_MULTIPLE_SPECTRA

        OBS_BLOCKS=np.load(path_output+'Observing_blocks.npy',allow_pickle=True).item()

        for i in OBS_BLOCKS.keys():
            # We identify which data will be associated with the OBS block
            mask=np.array([Files.index(file) for file in OBS_BLOCKS[i]],dtype='int')

            os.makedirs(path_output+i,exist_ok=True)

            print('Merging and substracting sky ',i)

#____________________________________________________________________

            MERGING_MULTIPLE_SPECTRA(path  = path_output+i+'/',
                                     OFiles= OBS_BLOCKS[i],
                                     MJD   = DATA['MJD'][mask],
                                     Wave  = DATA['WAVE'][mask],
                                     Flux  = DATA['FLUX'][mask],
                                     RAS   = DATA['RAS'][mask],
                                     DECS  = DATA['DECS'][mask],
                                     GLAT  = DATA['GLAT'][mask],
                                     GLON  = DATA['GLON'][mask],
                                     HTEL  = DATA['HTEL'][mask],
                                     NOT_SKY   = None if DATA.get('NOT_SKY') is None else DATA.get('NOT_SKY')[mask],
                                     OBJS      = None if DATA.get('OBJS')    is None else DATA.get('OBJS')   [mask],
                                     MAGS      = None if DATA.get('MAGS')    is None else DATA.get('MAGS')   [mask],
                                     function_used_to_merge=function_used_to_merge,
                                     SAVE      = True,
                                     SHOW      = False,
                                     SHOW_SKIES= False,
                                     window    = (min_wave,max_wave),
                                     AXIS_TO_MERGE=AXIS_TO_MERGE,
                                     wave_to_angstrom=wave_to_angstrom,
                                     )
        
#____________________________________________________________________
# If you are lazy as me you probably would check that RV works for one or two OB and then you will want to run all at once
if RADIAL_VEL and Across_OB:
    RADIAL_VEL=False

    OBS_BLOCKS=np.load(path_output+'Observing_blocks.npy',allow_pickle=True).item()

    listOB=list(OBS_BLOCKS.keys())
    print('\nGetting RV for each OB')

    for i,OB in enumerate(listOB):
        print(i,OB)
        path_folder=path_output+listOB[i]+'/'
    
        spectrums=np.genfromtxt(path_folder+'Oficial_list.list',dtype='U100')
        
        WAVE_TEMP,FLUX_TEMP=TEMPLATE_READER(complete_path_template)

        with open(path_package+'Getting_RV.py') as f:
            exec(f.read())# This will execute this code (look at it for details and tunning)

#____________________________________________________________________
# Whenever we want to run one of these methods I will ask you to pick one of the observing blocks
if (RADIAL_VEL and not Across_OB) or (CONTINUUM and not AUTOMATED_CONT):
        
    try:
        from tkinter import filedialog as fd
        path_folder:str=fd.askdirectory(title='Select the Observing Block that you want to work',
                                       initialdir=path_output) + '/'
    except:

        OBS_BLOCKS=np.load(path_output+'Observing_blocks.npy',allow_pickle=True).item()

        listOB=list(OBS_BLOCKS.keys())
        print('\nAccordingly to the following list:')

        for i,OB in enumerate(listOB):
            print(i,OB)
            
        OB=int(input('\n Please select which Observing block do you want to work with (e.g 0): '))
        print()
        path_folder=path_output+listOB[OB]+'/'
    
    spectrums=np.genfromtxt(path_folder+'Oficial_list.list',dtype='U100')

#____________________________________________________________________

    if RADIAL_VEL: 

        if WAVE_TEMP is None or FLUX_TEMP is None:
            WAVE_TEMP,FLUX_TEMP=TEMPLATE_READER(complete_path_template)

        with open(path_package+'Getting_RV.py') as f:
            exec(f.read())# This will execute this code (look at it for details and tunning)
#____________________________________________________________________

    if CONTINUUM and not AUTOMATED_CONT:
        from Continuum_fitting import cont_fit
        CONT=cont_fit(spectrums,                  # a list with the files that we will work with
                    path_folder+EXTRA_FOLDER_CONT,# The path where thes efiles 
                    wave_range=WAVE_RANGE,   # Wavelength range that will be used to perform the fit
                    Continue=True,           # Look for the file that contains previous adjustments in the folder
                    READ_FUNCTION=read_fits  # Function used to read files (just fits files) if None it will use defaul values
                    )
        
        # Once the complete OB is done it will save the corrected by continuum in a subfolder
        CONT.store_cont_sub(Force=False,)     # If you want to save it any time you run despite that you have not end the obs block change to Force=True
                                
#____________________________________________________________________

if CONTINUUM and AUTOMATED_CONT:
    from automated_cont_fit import automated_continuum,multiple_spectra_cont,default_read

    OBS_BLOCKS=np.load(path_output+'Observing_blocks.npy',allow_pickle=True).item()

    listOB=list(OBS_BLOCKS.keys())
    print('\nGetting an automatic continuum for each OB')

    for i,OB in enumerate(listOB):
        print(i,OB)
        path_folder=path_output+listOB[i]+'/'
    
        spectrums=np.genfromtxt(path_folder+'Oficial_list.list',dtype='U100')
    
        multiple_spectra_cont(spectra      = spectrums,                       # Iterable of name files containing spectra
                                path       = path_folder+EXTRA_FOLDER_CONT,   # path where the spectrums are allocated
                                wave_range = WAVE_RANGE,                      # wavelenth range (in case you want to remove the extremes)\n
                                )