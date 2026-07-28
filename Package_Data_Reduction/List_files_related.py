#----------------------------------------
# #Packages
import numpy as np
from prettytable import PrettyTable
from astropy.coordinates import SkyCoord
from astropy import units as u

#----------------------------

# This function looks for concatenated measures wihting a time range
# The idea behind this function is to identify measures that you want to combine  

def OBSERVING_BLOCK(path,              # path where the file with the info of the OB will be stored 
                    time,              # time in mJD or similar (must be in days)
                    epochs,            # Name of files/folders for easy identification of data that is associated with the OB
                    RA,                # List of Right Ascention associated with the specific times of measure (Variable time)
                    DEC,               # List of Declination     associated with the specific times of measure (Variable time)
                    not_sky=None,      # A List with a mask of the size of RA and DEC (len(RA)==len(DEC)==len(not_sky) and ) / if set just true imply that you dont have any sky in the given
                    OBS_BLOCK_TIME=12, # Range in hours that you want to consider as the same observing block (12 would work for separate by night) 
                    shapes=None,       # (Optional) A variable with the shapes of the data (might be useful to check)
                    Fname ='Observing_blocks.npy', # Name of output file (numpy kind of file)
                    SAVE_TXT= True     # It will create a defaul TXT file with the same infor printed
                    ) -> None:         # The function does not return nothing Just creates a file with the Observing blocks in a Dictionary stored with numpy likie format 
    # We sort the files by time
    idx   = np.argsort(time)
    time  = np.array(time)[idx]
    epochs = np.array(epochs)[idx]

    RA     = np.array(RA     ,dtype='object')[idx]
    DEC    = np.array(DEC    ,dtype='object')[idx]
    
    if not_sky is None:
        not_sky=True
    else:
        not_sky= np.array(not_sky,dtype='object')[idx]

    if not shapes is None: 
        shapes= np.array(shapes)[idx]   

    # Auxiliar variables
    index= np.arange(len(time))
    temp = True
    AUX  = [time[0]]

    DATA = dict()

    eliminate_space=lambda x: [str(i).replace(' ','') for i in x]
    # For TXT
    if SAVE_TXT:
        F=open(path+Fname.replace('.npy','')+'.txt','w')
    
    # We will look for start points in observing blocks
    for j,i in enumerate(AUX):
        # We look for files that were taken within the same nigth
        mask = temp & (abs(time-i)<OBS_BLOCK_TIME/24)  # 12 hours should be more than enough to separate one night in the observatory
        
        # We also will check that the configurations are the same (same RA and DEC)
        current_index=index[mask]
        for k in range(len(current_index)-1):
            # We will check that every RA/DEC not associated with sky measure must be contained in the next set of RA and DEC
            
            current_measure= SkyCoord(ra=RA[current_index[k]]  [not_sky[current_index[k]]]  *u.degree, dec=DEC[current_index[k]]  [not_sky[current_index[k]]]  *u.degree)
            next_measure   = SkyCoord(ra=RA[current_index[k+1]][not_sky[current_index[k+1]]]*u.degree, dec=DEC[current_index[k+1]][not_sky[current_index[k+1]]]*u.degree)

            indexes,d2d,d3d=current_measure.match_to_catalog_sky(next_measure)

            sep= d2d>0.5*u.arcsec # Threshold of difference in position 

            # If the lenght of the arrays (without considering the sky) and non of the targets was 
            if (not sum(sep)) and (len(current_measure)==len(next_measure)):
                ras,decs = True,True
        
            else: # In case that the arrays of RA DEC does not match or not the same targets
                ras,decs=False,False

            mask[current_index[k+1]]= (ras&decs) & mask[current_index[k]]

        # Storing Observing Block
        DATA[f'OB-{j}']= np.ravel(epochs[mask]) 

        # Display files in an easy way to look them
        table = PrettyTable()
        table.add_column("File_Name"       , epochs[mask] )
        table.add_column("Date"            , time[mask]  )
        table.add_column("Diff_in_Hour"    , abs(time[mask]-i)*24   )
        table.add_column("Diff_in_Minutes" , abs(time[mask]-i)*24*60)
        if not shapes is None: 
            table.add_column("Shapes_of_Data"  , eliminate_space(shapes[mask]))
        
        # Aligment of columns
        table.align["File_name"]      = "l"
        table.align["Date"]           = "l"
        table.align["Diff_in_Hour"]   = "l"
        table.align["Diff_in_Minutes"]= "l"

        # We update the times that has been checked
        temp= (temp & ~mask)
        if sum(temp)>0:
            AUX.append(time[temp][0])
        
        print('OBS ',j)            
        print(table)
        if SAVE_TXT:
            F.write('# OBS'+str(j)+'\n')
            F.write('#'+table.get_string(border=False)[1:])
            F.write('\n\n')
        print()
        print()



    np.save(path+Fname,DATA)  
    if SAVE_TXT:
        F.close()

    print(f'\n\nData per OB had been stored on\n {path+Fname}\n\n')
    print(f'This can be read with:\n np.load({path+Fname} ,allow_pickle=True).item()\n')


