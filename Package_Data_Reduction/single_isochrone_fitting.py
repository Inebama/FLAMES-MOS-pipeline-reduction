# Packages 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from tkinter import filedialog as fd
from prettytable import PrettyTable
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button

# Just for typing annotation
try:
    from typing import Iterable
except:
    from collections.abc import Iterable


class ISO_FIT:
    """
    Class thought to easyly fit parameters need on an Isochrone and retrieve atmospheric parameters
    Make sure your sample of Targets is a relatively clean sample 

    #################################
    # PARAMETERS : data type = default value # Meaning
    

        MAG       : Iterable ,               # Magnitude of stars to be used
        COLOR     : Iterable ,               # Color of stras to be used
        iso_path  : str      | None = None,  # If you want you can provide directly a path otherwise it will ask you to point one with an emergent window
        OBJ       : Iterable | None = None,  # [Optional] an iterable with some identifier to be add to the output file          
        RA        : Iterable | None = None,  # [Optional] an iterable with Right Ascention of the stars to be add in the output file
        DEC       : Iterable | None = None,  # [Optional] an iterable with the Declination of the stars to be add on the output file
        photometry: str  = 'Gaia',           # gaia or ubvir are allowed by default (does not matter upper or lower case)
        Header    : int  = 13,               # Number of columns of header in the isochrone (old version 13, PARSEC V2 ->14)
        Stages    : list = [                 # List Of stages that you want to see and fit according to the PARSEC isochrones 'label' keyword 
                            #0, # Pre Main Sequence # by default not used
                            1,  # Main Sequence
                            2,  # Sub Giant
                            3,  # RGB
                            4,  # CHEB
                            #5, # Red-CHEB          # by default not used
                            #6, # Blue-CHEB         # by default not used
                            7,  # E-AGB
                            8,  # TP-AGB
                            #9, # Post AGB          # by default not used
                            ],
                 ):
    
    """

    def __init__(self,
                 MAG       : Iterable ,               # Magnitude of stars to be used
                 COLOR     : Iterable ,               # Color of stars to be used
                 iso_path  : str      | None = None,  # If you want you can provide directly a path otherwise it will ask you to point one with an emergent window
                 OBJ       : Iterable | None = None,  # [Optional] an iterable with some identifier to be add to the output file          
                 RA        : Iterable | None = None,  # [Optional] an iterable with Right Ascention of the stars to be add in the output file
                 DEC       : Iterable | None = None,  # [Optional] an iterable with the Declination of the stars to be add on the output file
                 photometry: str  = 'Gaia',           # gaia or ubvir are allowed by default (does not matter upper or lower case)
                 Header    : int  = 13,               # Number of columns of header in the isochrone (old version 13, PARSEC V2 ->14)
                 Stages    : list = [                 # List Of stages that you want to see and fit according to the PARSEC isochrones 'label' keyword 
                                     #0, # Pre Main Sequence # by default not used
                                     1,  # Main Sequence
                                     2,  # Sub Giant
                                     3,  # RGB
                                     4,  # CHEB
                                     #5, # Red-CHEB          # by default not used
                                     #6, # Blue-CHEB         # by default not used
                                     7,  # E-AGB
                                     8,  # TP-AGB
                                     #9, # Post AGB          # by default not used
                                     ],
                 ):
        

        # Storing Parameters
        self.MAG      = MAG    # Target Magnitude
        self.COLOR    = COLOR  # Target Color

        self.OBJ      = OBJ      # Array of objject if given
        self.RA       = RA       # rigth ascention to identify targets
        self.DEC      = DEC      # declination to identify targets

        self.Stages = Stages # Evolutionary stages to be used
        self.HEADER = Header
        # Stablishes the photometric system and their absorption coeficients
        self.photometric_system(photometry)
        
        self.run_app(path=iso_path)

    #____________________________________________________________
    def run_app(self,path=None):
        try:
            plt.close()
        except:
            print('Not figure')
        # We read the isochrone and define the variables
        self.isochrone_read(path,False)

        # Create the figure 
        self.Figure_setup()
        plt.show()
    #____________________________________________________________
    def photometric_system(self,
                           system:str,
                           #
                           )-> None:
        """
        Function to set photometric system used 

        You must define it according to your photometric system 
        the Magnitude used in plot and interpolation will be Mag1 and the color will be defined as Mag2-Mag3
        
        You must also provide a absorption coefficient relative to visual absorption (Av) that will be used to estimate the reddening
        
        Name of variables on absorption coeficcient defined must be the same as the defined in system look in the function for clear example
        """
        
        # Set you photometric system and their combination of filters
        if system.upper()=='GAIA': # DR3 and EDR3
            self.system={'Mag1':'G_fSBmag',
                         'Mag2':'G_BP_fSBmag',
                         'Mag3':'G_RP_fSBmag'}
            self.absorption_coeff={'G_fSBmag':0.83627, 	 	
                                   'G_BP_fSBmag':1.08337,
                                   'G_RP_fSBmag':0.63439,
                                  }
        if system.upper()=='GAIA_OLD': # DR3 and EDR3 V1.2
            self.system={'Mag1':'Gmag',
                         'Mag2':'G_BPmag',
                         'Mag3':'G_RPmag'}
            self.absorption_coeff={'Gmag':0.83627, 	 	
                                   'G_BPmag':1.08337,
                                   'G_RPmag':0.63439,
                                  }
        # I use this for Stetson photometry (see Stetson et.al. 2019)
        # Old version Parsec
        if system.upper()=='UBVIR_OLD':
            self.system={'Mag1':'Imag',
                         'Mag2':'Vmag',
                         'Mag3':'Imag'}
            self.absorption_coeff={'Umag':1.55814,   
                                   'Bmag':1.32616,     
                                   'Vmag':1.00096,    
                                   'Rmag':0.80815,    
                                   'Imag':0.59893,    
                                   'Jmag':0.28688,    
                                   'Hmag':0.18103,    
                                   'Kmag':0.11265,
                                  }
        # New photometric system V2
        if system.upper()=='UBVIR':
            self.system={'Mag1':'I_fSBmag',
                         'Mag2':'V_fSBmag',
                         'Mag3':'I_fSBmag'}
            self.absorption_coeff={'U_fSBmag':1.55814,   
                                   'B_fSBmag':1.32616,     
                                   'V_fSBmag':1.00096,    
                                   'R_fSBmag':0.80815,    
                                   'I_fSBmag':0.59893,    
                                   'J_fSBmag':0.28688,    
                                   'H_fSBmag':0.18103,    
                                   'K_fSBmag':0.11265,
                                  }
        # Set your own photometric system    
        if system.upper()=='system_random':
            self.system={'Mag1':'magnitude1',
                         'Mag2':'magnitude2',
                         'Mag3':'magnitude3'}
            
            self.absorption_coeff={'magnitude1':1000,
                                   'magnitude2':10,
                                   'magnitude3':1,
                                   }
        
        ###########################
        self.abs_mag=self.absorption_coeff[self.system['Mag1']]
        self.abs_col=self.absorption_coeff[self.system['Mag2']]-self.absorption_coeff[self.system['Mag3']]
    #___________________________________________________________________________
    
    # Functions that defines the frame of the figure the interactive actions and the first plot
    def Figure_setup(self,id=170499):
        # We create the figure where we will be seeing the spectra and the continuum
        # We create figure and axes using GridSpec
        fig = plt.figure(id,figsize=(8, 8))
        gs  = GridSpec(9, 11,
                       height_ratios=[0.1,0.075,1,0.1, 0.1, 0.075, 0.15 ,0.075, 0.15], 
                       width_ratios =[0.075,1,0.15,0.075,0.075,0.1,0.075,0.075,0.1,0.075,0.075],
                       left   = 0.07, 
                       right  = 0.97, 
                       bottom = 0.08, 
                       top    = 0.95, 
                       wspace = 0.0, 
                       hspace = 0.0)

        # Axes for spectrum plot, sliders, and spectra list
        ax_CMD = fig.add_subplot(gs[0:4, 0:5])    # For CMD
       

        # Avail MH 
        ax_MH      = fig.add_subplot(gs[2, 6:8])        
        ax_prev_mh = fig.add_subplot(gs[4, 6])
        ax_next_mh = fig.add_subplot(gs[4, 7])

        # Avail AGE 
        ax_AGE      = fig.add_subplot(gs[2, 9:11])        
        ax_prev_age = fig.add_subplot(gs[4, 9])
        ax_next_age = fig.add_subplot(gs[4, 10])

        # Button Open New isochrone

        ax_new_iso = fig.add_subplot(gs[0, 6:])

        # Button run infere stellar parameters

        ax_kiel = fig.add_subplot(gs[6:, 6:])

        # Slider
        ax_slider1 = fig.add_subplot(gs[6, 1])      # AV
        ax_slider2 = fig.add_subplot(gs[8, 1])      # UM
        # Buttons for incrementing and decrementing sliders
        ax_inc1 = fig.add_subplot(gs[6, 4])
        ax_dec1 = fig.add_subplot(gs[6, 3])
        ax_inc2 = fig.add_subplot(gs[8, 4])
        ax_dec2 = fig.add_subplot(gs[8, 3])


        # display of intial plot and labels
        self.data_scatter = ax_CMD.scatter(self.COLOR, self.MAG,label='Targets',marker='o',color='gray',s=5,alpha=0.7)
        self.iso_scatter  = ax_CMD.scatter(self.Fit_Color_iso,self.Fit_Mag_iso,color='red',marker='+',label='Isochrone',s=10)
        
        ax_CMD.set_xlabel(f'{self.system["Mag2"]}-{self.system["Mag3"]}')
        ax_CMD.set_ylabel(f'{self.system["Mag1"]}')

        # just for display purposes
        ax_CMD.set_xlim(min(self.COLOR),max(self.COLOR))
        ax_CMD.set_ylim(max(self.MAG),min(self.MAG))

        ax_CMD.legend(loc='lower right',bbox_to_anchor=(1.01,0.97),framealpha=0)

        # Set up sliders
        self.AV_slider  = Slider(ax_slider1, 'Visual\nAbsorption', 0 , 3 , valinit=self.AV, valstep=0.01)
        self.UM_slider  = Slider(ax_slider2, 'Distance\nModulus' , 0 , 25, valinit=round(self.UM,1), valstep=0.02)
        self.MH_slider  = Slider(ax_MH ,     'MH'                , self.MH[0] ,self.MH[-1] , valinit=self.current_MH , valstep=self.MH,orientation='vertical')
        self.AGE_slider = Slider(ax_AGE,     'logAge'            , self.AGE[0],self.AGE[-1], valinit=self.current_AGE, valstep=self.AGE,orientation='vertical')
    
        # Set up buttons 
        # Buttons for each slider (increment and decrement)
        self.btn_inc1   = Button(ax_inc1, '+')
        self.btn_dec1   = Button(ax_dec1, '-')
        self.btn_inc2   = Button(ax_inc2, '+')
        self.btn_dec2   = Button(ax_dec2, '-')
        
        self.btn_mh_inc = Button(ax_next_mh,'+')
        self.btn_mh_dec = Button(ax_prev_mh,'-')
        self.btn_age_inc= Button(ax_next_age,'+')
        self.btn_age_dec= Button(ax_prev_age,'-')


        self.btn_inc1.on_clicked(lambda event: self.increment_slider(self.AV_slider,  0.01))
        self.btn_dec1.on_clicked(lambda event: self.increment_slider(self.AV_slider, -0.01))
        self.btn_inc2.on_clicked(lambda event: self.increment_slider(self.UM_slider,  0.02))
        self.btn_dec2.on_clicked(lambda event: self.increment_slider(self.UM_slider, -0.02))
        
        self.btn_mh_inc.on_clicked(lambda event:  self.increment_slider_mh_age(self.MH_slider,  self.MH,  1))
        self.btn_mh_dec.on_clicked(lambda event:  self.increment_slider_mh_age(self.MH_slider,  self.MH, -1))
        self.btn_age_inc.on_clicked(lambda event: self.increment_slider_mh_age(self.AGE_slider, self.AGE, 1))
        self.btn_age_dec.on_clicked(lambda event: self.increment_slider_mh_age(self.AGE_slider, self.AGE,-1))

  
        # Figure name will indicate the current isochrone
        ax_CMD.set_title(f'MH={self.current_MH:.2f} dex, Age={10**self.current_AGE:.2e} yr', fontsize=12)

        self.fig   =fig
        self.ax_CMD=ax_CMD

        # Connect sliders to the update function
        self.AV_slider.on_changed(self.update)
        self.UM_slider.on_changed(self.update)
        self.MH_slider.on_changed(self.update_mh_age)
        self.AGE_slider.on_changed(self.update_mh_age)

        # Special Buttons
        self.btn_new_iso=Button(ax_new_iso,'Open other Isochrone')
        self.btn_kiel_dg=Button(ax_kiel,   'Obtain\n stellar parameters')

        self.btn_new_iso.on_clicked(lambda event: self.run_app())
        self.btn_kiel_dg.on_clicked(lambda event: self.match_nearest())
        
    #___________________________________________________________________
    # Function to update the plot
    def update(self,event=None):

        self.UM=self.UM_slider.val
        self.AV=self.AV_slider.val

        self.fitted_isochrone()

        self.iso_scatter.set_offsets(np.c_[self.Fit_Color_iso,self.Fit_Mag_iso])

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
    
    # Function that update the values of the isochrone for the current input
    def fitted_isochrone(self):
        self.Fit_Mag_iso  =self.Mag_iso  +self.UM + self.AV*self.abs_mag
        self.Fit_Color_iso=self.Color_iso+ self.AV*self.abs_col

    # Function that updates the mh and age in case of have it
    def update_mh_age(self,event=None):
        self.current_AGE=self.AGE_slider.val
        self.current_MH =self.MH_slider.val
        self.change_iso(update=False)
        self.renew_points_iso()

    # Update the plto of the isochrone from zero (need if the number of points change)
    def renew_points_iso(self):
        # Change title
        self.ax_CMD.set_title(f'MH={self.current_MH:.2f} dex, Age={10**self.current_AGE:.2e} yr', fontsize=12)

        self.iso_scatter.remove()
        self.iso_scatter  = self.ax_CMD.scatter(self.Color_iso,self.Mag_iso,color='red',marker='+',label='Isochrone',s=10)
        self.update()

    # Increment/Decrement button logic (keep values within slider limits)
    def increment_slider(self,slider:Slider, amount:float):
        new_val = np.clip(slider.val + amount, slider.valmin, slider.valmax)
        slider.set_val(new_val)
    
    # Increment/Decrement button for MH and AGE
    def increment_slider_mh_age(self,slider:Slider,LIST:list, amount:int):
        n=LIST.index(slider.val)
        try:
            slider.set_val(LIST[n+amount])
        except IndexError:
            slider.set_val(LIST[0])

    #___________________________________________________________________________
    def isochrone_read(self,
                       path  : str  = None,
                       update: bool = True,
                       #
                       ) -> None:
        """
        Function thougth to work with PARSEC ISOCHRONES

        It will use the stages defined on __init__  

        it has the option update for the firt iteration 
        
        """
        # If none path is give you will have to provide one
        if path is None:
            self.iso_path = self.ask_path()
        else:
            self.iso_path=path
        
        # Reading Isochrone
        Iso=np.genfromtxt(self.iso_path,names=True,skip_header=self.HEADER,comments='#')

        # We just keep the values of the stages of our interest
        aux=True
        for stage in self.Stages:
            aux= (aux & (Iso['label']!=stage))

        Iso=Iso[~aux]

        # Storing the isochrone parameters
        self.Iso=Iso # The isochrone but just the stages that we care about

        # Metallicity and AGE
        self.MH =list(set(Iso['MH']))
        self.AGE=list(set(Iso['logAge']))

        self.MH.sort()
        self.AGE.sort()

        # Stablish initial value to display
        self.current_MH =self.MH[0]
        self.current_AGE=self.AGE[0]

        # set up variables used in the plots and interpolations
        self.change_iso(update)

        # Init values for display purposes
        AV= (np.mean(self.COLOR)-np.mean(self.Color_iso))/self.abs_col
        if AV<0: # you can not have negative absorption
            AV=0

        self.AV= AV 
        self.UM= (np.mean(self.MAG  )-np.mean(self.Mag_iso)  )

        # Shift the isochrone near the data
        self.fitted_isochrone()
    
   
            

    #___________________________________________________________________
    # Function to update the values of the isochrone
    def change_iso(self,update=True):
        current_iso= ((self.Iso['MH']    == self.current_MH )&
                      (self.Iso['logAge']== self.current_AGE))

        self.current_iso=self.Iso.copy()[current_iso]

        self.Mag_iso  =  self.current_iso[self.system['Mag1']]
        self.Color_iso= (self.current_iso[self.system['Mag2']]-self.current_iso[self.system['Mag3']])
        
        if update:
            self.update()

    #______________________________________________________________________
    # Funtion to ask for path if you desire to open another isochrone
    def ask_path(self) -> str:
        path: str =fd.askopenfilename(title='Select the file of PARSEC isochrone') 
        return path

    #______________________________________________________________________
    # Function to add points in between of isochrone to improve representation 
    def increment_points(self,
                         n_points=10000,
                         max_jump_allowed_sigmas=4,
                         #
                         ) -> None:

        # We will work with the current Isochrone
        Iso=self.current_iso
        NNN=len(Iso)
        AUX1=np.arange(NNN)
        
        Mag_iso  = self.Fit_Mag_iso   
        Color_iso= self.Fit_Color_iso

        distance= ((Mag_iso[1:]-Mag_iso[:-1])**2+(Color_iso[1:]-Color_iso[:-1])**2)**0.5
        mean,std=np.nanmean(distance),np.nanstd(distance)

        # At least we will always have the entire isochrone
        mask = distance>mean+std*max_jump_allowed_sigmas
        idx  = [-1]+list(AUX1[:-1][mask])+[NNN]

        distance[mask]=0
        distance=[0]+list(distance)
        total_d=np.nansum(distance)

        store_iso=dict()
        for i in range(len(idx)-1):

            lower,upper=idx[i]+1,idx[i+1]+1

            distributed_points=int(n_points*(np.nansum(distance[lower:upper])/total_d))
            AUX2=np.linspace(AUX1[lower],AUX1[:upper][-1],distributed_points)

            for key in Iso.dtype.names:
                if i==0:
                    store_iso[key]=[]
                var=interp1d(AUX1[lower:upper],Iso[key][lower:upper])(AUX2)
                store_iso[key]+=list(var)

                if key=='label':
                    Iso[key]=np.round(Iso[key],0).astype(int)
        
        for key in Iso.dtype.names:
            store_iso[key]=np.array(store_iso[key])

        self.current_iso= store_iso
        self.Mag_iso    =  self.current_iso[self.system['Mag1']]
        self.Color_iso  = (self.current_iso[self.system['Mag2']]-self.current_iso[self.system['Mag3']])
        self.fitted_isochrone() # This updates the variable used in the next function
        

    #____________________________________________________________________________________
    # Function that match the nearest isochrone point to the data
    # It will store the data in a file and after it will plot the interpolated results that will pop up in another window
    def match_nearest(self,Output_file:str='Isochrone_values.txt') -> None:
            # We increase the number of point in the isochrone
            self.increment_points()
            
            self.waiting_figure()


            # Target data
            MAG   = self.MAG
            COLOR = self.COLOR

            # Increased data of isochrone
            iso_mag=self.Fit_Mag_iso
            iso_col=self.Fit_Color_iso

            with open(Output_file,'w') as file:

                file.write(f'# Isochrone Used ={self.iso_path}\n')
                file.write(f'# MH             ={self.current_MH}\n')
                file.write(f'# logAGE         ={self.current_AGE}\n')
                file.write(f'# AV             ={self.AV}\n')
                file.write(f'# UM             ={self.UM}\n')
                file.write(f'# abs_coeff_Mag  ={self.abs_mag}\n')
                file.write(f'# abs_coeff_Color={self.abs_col}\n')
          

                # Display files in an easy way to look them
                table = PrettyTable()

                if not self.OBJ is None:
                    table.add_column('Obj',self.OBJ)
                if not self.RA is None:

                    table.add_column('RA',self.RA)
                if not self.DEC is None:
                    table.add_column('DEC',self.DEC)

                table.add_column('MAG'  ,MAG)
                table.add_column('COLOR',COLOR)


                indexes=[]
                for i in range(len(MAG)):
                    mag,color=MAG[i],COLOR[i]
                    dist=((iso_mag-mag)**2+(iso_col-color)**2)**0.5
                    indexes.append(np.nanargmin(dist))
                
                indexes=np.array(indexes)

                for key in self.current_iso.keys():
                    table.add_column(key,self.current_iso[key][indexes])

                file.write('#'+table.get_string(border=False)[1:])
            #-------------------
            self.kiel_diagram(indexes)
    #_____________________________________________________________________
    # warning plot to explicitly say that the computer is doing somenthing
    def waiting_figure(self,id=1704) -> None:
        plt.close()
        FIG= plt.figure(id,figsize=(12,4))
        ax =FIG.add_subplot(xticks=[], yticks=[])
        # The first word, created with text().
        text = ax.text(-.1, .5, 'The stellar parameters are in process of estimate them,\n please wait till another window pop up,\n A warninn might appear that matplotlib is not responding but just wait', color="k",fontsize=20)
        ax.spines[:].set_visible(False)
        plt.show(block=True)
        plt.pause(0.1)

        print("\n\nRunning the stellar parameters ...\n\n")

    # Figure to create the Kiel diagramm must give the info of the warning plot and the indexes according to the current iso
    def kiel_diagram(self,indexes,id=1704):
        plt.close()

        self.Figure_setup()
        plt.figure(id)
        plt.scatter(10**self.current_iso['logTe'][indexes],self.current_iso['logg'][indexes])
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.xlabel('Teff')
        plt.ylabel('Log g')
        plt.title('Kiel Diagram')
        plt.show(block=False)


if __name__=='__main__':
    from astropy.io import fits

    path_data='/home/ian/Desktop/Data_base/'
    diff_redd='Diferrential_red_Pancino/Pancino_2024_diferential_reddening.fit'
    catalogue='FullCatalogue.fits'
    # My data
    with fits.open(path_data+catalogue) as hdul:
        Prob= hdul[1].data.Prob
        mask= Prob>0.99
        I   = hdul[1].data.I
        V   = hdul[1].data.V
        mask= mask & (abs(V-I)<5)
        V   = V[mask]
        I   = I[mask]
        id  = hdul[1].data.Id[mask]
        ra  = hdul[1].data.ra[mask]
        dec = hdul[1].data.dec[mask]
        gmag= hdul[1].data.Gmag[mask]
        bprp= hdul[1].data.BpRp[mask]

    from scipy.interpolate import LinearNDInterpolator
    # Reddening map
    my_cluster='NGC104'
    with fits.open(path_data+diff_redd) as hdul:
        # Mean of the cluster
        Cluster=hdul[1].data.Cluster
        mask   = (Cluster==my_cluster)
        # mean_C = hdul[1].data['E_B-V_'][mask]
        mean_C = 0.03 #Brogaard 2017
        
        # Differential
        Cluster = hdul[2].data.Cluster
        flag    = hdul[2].data.Flag
        mask    = (Cluster==my_cluster)&(~((flag!=0)&(flag!=1)))
        # mask    = (Cluster==my_cluster)&(flag==0)
        RA      = hdul[2].data.RAICRS[mask]
        DEC     = hdul[2].data.DEICRS[mask]
        EBV     = hdul[2].data["dE_B-V_"][mask]#+mean_C
        EBV-=EBV.min()

    # Create the reddening map
    intp=LinearNDInterpolator(list(zip(RA, DEC)),EBV)

    newerr=intp(ra,dec)

    
    B_abs_av  = 1.32616 
    V_abs_av  = 1.00096
    I_abs_av  = 0.59893 
    EBV_abs_av=(B_abs_av-V_abs_av)
    print(EBV_abs_av)

    g_abs_av  = 0.83627 	 	
    bp_abs_av = 1.08337
    rp_abs_av = 0.63439

    # # plt.figure()
    # # plt.scatter(V-I,I)
    # # plt.gca().invert_yaxis()

    V_or=V.copy()
    I_or=I.copy()

    #Correction redeninig
    V=V-((V_abs_av/EBV_abs_av)*newerr) 	
    I=I-((I_abs_av/EBV_abs_av)*newerr) 	
    	
    NAN_mask=~np.isnan(V+I)    
    V=V[NAN_mask]
    I=I[NAN_mask]

    V_or=V_or[NAN_mask]
    I_or=I_or[NAN_mask]

    gmag=gmag-((g_abs_av/EBV_abs_av)*newerr) 	
    bprp=bprp-(((bp_abs_av-rp_abs_av)/EBV_abs_av)*newerr) 	
    
    mask=~np.isnan(gmag+bprp)

    print(np.nanmin(newerr),np.nanmax(newerr))
    plt.figure()

    IDX=np.argsort(newerr[~np.isnan(newerr)])[::-1]

    SC=plt.scatter(ra[IDX],dec[IDX],c=newerr[IDX],vmin=np.nanmin(newerr),vmax=np.nanmax(newerr),s=3)
    plt.colorbar(SC)
    plt.gca().invert_yaxis()

    plt.figure()
    SC=plt.scatter(RA,DEC,c=EBV,vmin=EBV.min(),vmax=EBV.max())
    plt.colorbar(SC)
    plt.gca().invert_yaxis()

    plt.show()
    import sys
    sys.exit()
    
    # plt.figure()
    # plt.scatter(V-I,I)
    # plt.gca().invert_yaxis()

    # plt.figure()
    # plt.scatter(I_or-I,I)
    # plt.gca().invert_yaxis()

    # plt.figure()
    # plt.scatter(V-I,(V-I)-(V_or-I_or))
    # plt.gca().invert_yaxis()

    # plt.show()

    # import sys
    # sys.exit()


    mask=()

    iso1='isochrone.txt'
    iso2='ParseV2-10-15Gyr__-0.89_+0.13MH_n0.4.txt'
    iso3='ParsecV2_12.86-13.26Gyr_-0.89_-0.1MH_n0.4.txt'
    iso3='ParsecV2_12.86-13.26Gyr_-0.89_-0.1MH_n0.65.txt'
    # iso3='ParsecV2_10-13.25Gyr_-0.8_-0.45MH_n0.65.txt'
    iso3='ParsecV2_10-15Gyr_-0.8_-0.45MH_n0.65.txt'
    # iso3='ParsecV2_9-14Gyr_-0.75_-0.4MH_n0.46_GAIA.txt'
    # iso3='ParsecV2_9-14Gyr_-0.75_-0.4MH_n0.46.txt'
    # ISO_FIT(MAG=gmag[mask],COLOR=bprp[mask],OBJ=id[mask],RA=ra[mask],DEC=dec[mask],iso_path='/home/ian/Downloads/'+iso3,photometry="gaia",Header=14)

    ISO_FIT(MAG=I,COLOR=V-I,OBJ=id,RA=ra,DEC=dec,iso_path='/home/ian/Downloads/'+iso3,photometry="ubvir",Header=14)

# MH =-0.59
# AGE=10.16879
# AV =0.01
# UM =13.219999999999997
# abs coeff Mag  =0.59893
# abs coeff Color=0.4020300000000001