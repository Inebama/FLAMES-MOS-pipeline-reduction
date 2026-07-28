# Code developed on Python 3.11.7 by Ian Baeza 
# Packages 
import matplotlib.pyplot as plt         # Plot results
import matplotlib.gridspec as gridspec  # Plot results
from matplotlib import lines            # Plot results
import numpy as np                      # Management of data and useful biult in functions
from scipy.interpolate import interp1d  # Interpolate the template spectrum
from scipy.optimize import curve_fit    # To get errors of the fit

# plt.style.use('dark_background') # My preferences in color

# To create Gif (not necessary)
# from PIL import Image                  
# import os

######################################################################################

# It is assumed that you will provide an Observed spectrum that has been treated before, this means to remove edges of low efficiency and if you provide an spectrum with gaps with no values of wavelenght associated within the gaps 

class CCF:
    """
    ____________________________________________________________________
    #Variables Description and usage

    ---------------------
    template_wavelength: list or array of wavelenght of template spectrum  (must be the same size as template_flux, and same units of observed_wavelength)
    
    #WARNING: It is need to have a larger range than of the wavelength than the range in observed_wavelenght
    
    ---------------------
    template_flux:       list or array of flux of template spectrum        (must be the same size as template_wavelenght)
    
    #WARNING: It is recomended to have the flux normalized, or at least scaled to values similar to 1 (dividing by median, etc.)
    
    ---------------------
    observed_wavelength: list or array of wavelenght of Observed spectrum  (must be the same size as observed_flux, and same units of template_wavelength)
    
    #WARNING: It is required to have a smaller range in wavelength than the range in template_wavelenght 
             (min(template_wavelength)<min(observed_wavelength) and max(observed_wavelength)<max(template_wavelength))
             This will limit the range of RV in which you are able to shift the spectra and propely interpolate
    ---------------------
    observed_flux: list or array of wavelenght of Observed spectrum  (must be the same size as observed_wavelength)
    
    #WARNING: It is recomended to have the flux normalized, or at least scaled to values similar to 1 (dividing by median, etc.)
    
    ---------------------
    RV_range:            list or an array, with the minimum and maximum value of the Radial velocity to be tested if not defined if will test fom the limits relative to the template and observed spectrum
    
    ---------------------
    IGNORE:              bool (True or False) , if you want to avoid the check of relative wavelength of the template and the observed spectrum and flux scaled to 1
    
    ---------------------
    PRINT:               bool, if you want to see printed on terminal information during the checks 
    
    ---------------------
    Cosmic_rays:         bool, In case that your data has cosmic rays (likely) set this command to true, it will search for relevant cosmic rays and will eliminate them for the computation of the CCF (this will also cut emission lines)
    ---------------------
    Check_template_relative_to_obs: Bool ask if you want to plot the max red/blue shifted of the templated
                                    Based on the given or calculated RV range plot the max and min RV shift
                                    of the template in order to check that will alway cover the obs spectra
    _______________________________________________________________________
    """

    def __init__(self,template_wavelength,  
                      template_flux, 
                      observed_wavelength, 
                      observed_flux,
                      RV_range=None,
                      IGNORE  = False,
                      PRINT   = False,
                      Cosmic_rays=True,
                      Check_template_relative_to_obs=False,
                      ):
        
        
        # Store the variables (these migth be modified below)
        self.t_wave   = np.array(template_wavelength )
        self.t_flux   = np.array(template_flux       ) 
        self.obs_wave = np.array(observed_wavelength ) 
        self.obs_flux = np.array(observed_flux       )
        self.RV_range = RV_range

        # Minimum error that we will have given the step in wavelength provided (half of the characteristical step in wavelength with respect to central wavelength converted to RV)
        self.ERR_init = np.median(self.obs_wave[1:]-self.obs_wave[:-1])/np.median(self.obs_wave)*299_792.458/2

        #------------------------------------------------------------
        self.Computed = False    # This means that RV_range was given
                                 # if not it will be re-defined below  
        #------------------------------------------------------------
        #______________
        # I do not recomend to ignore unless you are sure 
        # that your range of your template could be shifted 
        # in the ranges of RV provided, unexpected behaviour 
        # when interpolate could lead to bad measures of RV
            
        if not IGNORE: 
            #--------------------------------------------------------
            # Check if the fluxes are scaled to 1 (or normalized)

            MTflux= np.median(self.t_flux  )
            MOflux= np.median(self.obs_flux)

            if not (0.5<MOflux<1.5):
                self.obs_flux/=MOflux
                print('The Observed Flux has been divided by the median value')

            if not (0.5<MTflux<1.5):
                self.t_flux/=MTflux
                print('The Template Flux has been divided by the median value')

            #--------------------------------------------------------
            # Check the limits to test the correlation
            Check_RV  = lambda obs,temp: ((obs/temp)-1)*299792.458

            Check_wave= lambda temp,RV: temp*(1+RV/299792.458)

            MinT,MinO = min(template_wavelength),min(observed_wavelength)
            MaxT,MaxO = max(template_wavelength),max(observed_wavelength)
            

            if RV_range is None:
                # If your template is smaller in range of wavelength than your obs spectra 
                # I will stop unless you provide a range of RV so I can cut properly the obs spectra
                if MinT>MinO or MaxT<MaxO:
                    raise Exception('The min or max of your Obs spectra exceeds you template, unless you provide another Template or a range of RV we can not continue (and in this late case I will have to cut your Obs Spectra accordilgly)')
                #-------
                self.Computed=True

                RV_min= Check_RV(MaxO,MaxT) # Blue shift need to reach the border
                RV_max= Check_RV(MinO,MinT) # Red  shift need to reach the border 

                RV_range=[RV_min,RV_max]

                if PRINT:
                    print(f'Given the template and the Observed wavelength\n you will be able to test RV in the range of  \n\t[{RV_min:.2f},{RV_max:.2f}] km/s')
                
                self.RV_range=RV_range

            else:
                # Making sure that you give min and max RV ordered
                if len(RV_range)!=2 or RV_range[0]>=RV_range[1]:
                    raise Exception('RV_range Variable is expected to have a min and max value of Radial velocities that will be tested ordered like [min,max], please enter correct values or leave in blank')
                
     
                minWave=Check_wave(MinT,RV_range[-1]) # max +shift of the shortest wavelength
                maxWave=Check_wave(MaxT,RV_range[0])  # max -shift of the longest  wavelength

                if MinO<minWave or maxWave<MaxO:
                    print(f'In order to test RV=[{RV_range[0]:.2f},{RV_range[-1]:.2f}] I\'m gona cut the obs spectra since your template is to short in range of wavelength')
                    mask= (minWave<self.obs_wave)&(self.obs_wave<maxWave)

                    self.obs_wave=self.obs_wave[mask]
                    self.obs_flux=self.obs_flux[mask]

                # If you provided a much larger template we will trim it in order to make more efficient the calculation
                if MinT<minWave*0.97 or maxWave*1.05<MaxT:
                    mask= (minWave*0.99<self.t_wave)&(self.t_wave<maxWave*1.01)

                    self.t_wave=self.t_wave[mask]
                    self.t_flux=self.t_flux[mask]
                elif PRINT:
                    print(f'The range you provide of RV=[{RV_range[0]},{RV_range[1]}]\n is within the range relative to your Observed Spectra ')
                    print(f'On maximum Min Temp ={minWave:.2f}< Min Obs {MinO:.2f}')
            
            if Check_template_relative_to_obs:
                plt.plot(self.t_wave*(1+self.RV_range[0]/299792.458),self.t_flux,color='blue',label='Temp Blue Shift')
                plt.plot(self.t_wave*(1+self.RV_range[1]/299792.458),self.t_flux,color='red', label='Temp Red Shift')
                plt.plot(self.obs_wave,self.obs_flux,color='black',label='Obs Spectra')
                plt.legend()
                plt.show()
            #--------------------------------------------------------
        #______________
        elif RV_range is None:
            raise Exception('If you want to avoid the checks with IGNORE you have to provide the RV_range')
        
        #------------------------------------------------------------
        if Cosmic_rays:
            p99 =np.percentile(self.obs_flux,99)
   
            # All the flux that it is over p99 will be removed
            mask  = (self.obs_flux)>p99
            self.obs_flux[mask]=p99
            
          


    #################################################################
    # Auxiliar functions

    # Gaussian function
    def gaussian(self,x, x0, sigma, a=1):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
    
    # Parabola function (concave downward)
    def parabola(self,x,x_center,a,c=1):
        return -a* (x - x_center)**2 + c
    
    # Scale array of values between 0 and 1
    def scale_to_1(self,x):
        return (x-min(x))/(max(x)-min(x))

    # Check values of fits and get error on the measure (you can print the relevant info with this)
    def get_errs_info(self,PRINT=False,ALL_info=False):

        RV=self.RV
        center1,err1,popt1,cov1=self.params_F_gauss
        center2,err2,popt2,cov2=self.params_F_parab
        center3,err3,popt3,cov3=self.params_parabola
        center4,err4,popt4,cov4=self.params_gauss

        CENTERS= [RV,center1,center2,center3,center4]
        ERRS   = [   err1   ,err2   ,err3   ,err4   ]

        M_cent   = np.nanmedian(CENTERS)
        std_cent = np.nanstd(CENTERS)
        M_errs   = np.nanmedian(ERRS)

        ERR = np.sqrt(std_cent**2+M_errs**2+self.ERR_init**2)

        if PRINT:
            print('\n\n-----------------------------------------------------------------------------------------\n')
            print('If any value is in "nan" this means that that method didn\'t converge \n')
            print( 'Fixed shapes (This means that we forced to perform the fit with max value equal 1)\n')
            print(f'Fixed Gaussian got center on RV = {center1:.3f}+/-{err1:.3f}')
            print(f'Fixed Prabola  got center on RV = {center2:.3f}+/-{err2:.3f}')
            print('\n Non fixed Shapes\n')
            print(f'Gaussian got center on RV = {center3:.3f}+/-{err3:.3f}\t Height = {cov3[-1,-1]:.3f}')
            print(f'Gaussian got center on RV = {center4:.3f}+/-{err4:.3f}\t Height = {cov4[-1,-1]:.3f}')
            print(f'\nStatistics\n')
            print(f'Median of RVs          = {M_cent:.3f}')
            print(f'Standard Deviation RVs = {std_cent:.3f}')
            print(f'Median Errors          = {M_errs:.3f}')
            print('\n\nFound in the CCF and Final ERR\n')
    
            print(f'Central RV = {RV:.3f} +/-{ERR:.3f}')

            print('\n\nFinal Error= ((Standard Deviation RVs)**2+(Median Errors)**2)**0.5\n')
            print('\n-----------------------------------------------------------------------------------------\n\n')
        #----
        if ALL_info:
            return RV,ERR,M_cent,std_cent,M_errs
        #----   
        return ERR



    #################################################################
    
    # Main function to compute correlation
    # Shift the template interpolate it and evalueate in the wavelength of obs spectra
    def cross_correlation(self,RV):

        template_wavelength  = self.t_wave
        template_flux        = self.t_flux
        observed_wavelength  = self.obs_wave
        observed_flux        = self.obs_flux
          
        # Shift the observed wavelength according to the radial velocity
        shifted_template_wavelength = template_wavelength * (1 + RV / 299792.458)
        
        # Interpolate the shifted template flux
        interp_flux = interp1d(shifted_template_wavelength, template_flux, bounds_error=False, fill_value=0)

        # Calculate the correlation 
        synth_flux= interp_flux(observed_wavelength)
        X=observed_flux * synth_flux
        correlation   = np.sum(X)   # /np.sum(observed_flux) # Normalize?

        # Template match method
        # template_match= np.sum(abs(observed_flux-synth_flux))
        # return correlation/template_match 

        return correlation

    #######################################################################
    
    def find_rv(self,steps3    = 300, # Steps of to test in the range given when Manual=True or used in the third iteration
                     steps1    = 100, # Steps used in first iteration 
                     steps2    = 100, # Steps used in first iteration
                     Tolerance = 90,  # Min Relative height % between min and max of values of correlation to be considered for the fit of gaussian an polinomy, to get errors
                     Manual    = False# Indicates that you don't go through the iterations and directly calculates the steps with just one reoslution (steps3)
                ):
        """
        Main method to derive radial velocities
        if the range was given you could set the number of steps 
        to search in the range given setting steps3  
        If the range was automatically computed by the init function 
        It will look efficiently in 3 iterations increasing precision
        The number of each iteration could be directly modified when calling this function
    
        Arguments of Function:data type = Default

        steps3   :int  = 300, # Steps of to test in the range given when Manual=True or used in the third iteration
        steps1   :int  = 100, # Steps used in first iteration 
        steps2   :int  = 100, # Steps used in first iteration
        Tolerance:int  = 90,  # Min Relative height % between min and max of values of correlation to be considered for the fit of gaussian an polinomy, to get errors
        Manual   :bool = False# Indicates that you don't go through the iterations and directly calculates the steps with just one reoslution (steps3)

        This function will return 4 variables 

        Sampled RV
        Cross Correlation mesured at the sampled RV

        Best RV based on CCF
        Error of RV measured

        """
        #--------------------
        # Range given or computed
        RV_r = self.RV_range

        #--------------------
        # If an Range of velocities was given
        if not self.Computed and Manual:
            RV_sample  = np.linspace(RV_r[0],RV_r[-1],steps3) # Start in the whole range
            Correlation= np.array(list(map(self.cross_correlation,  RV_sample)))
            
            RV = RV_sample[Correlation.argmax()]
        #--------------------
        # If the range of velocities are the availables 
        elif self.Computed or not Manual:
            # Mapping the value of correlation iter 1

            RV_sample1  = np.linspace(RV_r[0],RV_r[-1],steps1) # Start in the whole range
            result1     = np.array(list(map(self.cross_correlation,  RV_sample1)))
            peak,width1 = self.find_peaks(result1,SHOW=False)

            try:
                # Mapping the value of correlation iter 2
                RV_sample2  = np.linspace(RV_sample1[width1[0][0]],RV_sample1[width1[0][1]],steps2)
                result2     = np.array(list(map(self.cross_correlation,  RV_sample2)))
                peak,width2 = self.find_peaks(result2,SHOW=False)


                # Mapping the value of correlation iter 3
                RV_sample3  = np.linspace(RV_sample2[width2[0][0]],RV_sample2[width2[0][1]],steps3)
                result3     = np.array(list(map(self.cross_correlation,  RV_sample3)))

                # Best Value
                RV = RV_sample3[result3.argmax()]

                # Concadenate the results 
                RV_sample = np.concatenate((RV_sample1[:width1[0][0]],
                                            RV_sample2[:width2[0][0]],
                                            RV_sample3,
                                            RV_sample2[ width2[0][1]+1:],
                                            RV_sample1[ width1[0][1]+1:]))

                Correlation= np.concatenate((result1[:width1[0][0]],
                                            result2[:width2[0][0]],
                                            result3,
                                            result2[ width2[0][1]+1:],
                                            result1[ width1[0][1]+1:]))

            except:
                RV_sample  = np.linspace(RV_r[0],RV_r[-1],steps1+steps2+steps3) #Look in the whole range if somenthing went wrong
                Correlation= np.array(list(map(self.cross_correlation,  RV_sample)))
                peak,width = self.find_peaks(result1,SHOW=False)
                
                RV = RV_sample[Correlation.argmax()]

            

            #############################

            
        #--------------------
        # Scaling Values of Correlation between 0 and 1
        # This will make easier to fit a function in order to obtain the error in RV

        Correlation=self.scale_to_1(Correlation)

        # Store the calculations made
        self.RV_sample   = RV_sample.copy()
        self.Correlation = Correlation.copy()
        self.RV          = RV    
        #--------------------
        try:
            # We look the peak of the correlation
            peak,width = self.find_peaks(Correlation,tolerance=Tolerance,SHOW=False)
            
            # We will save these values for plot purposes and easier usage
            idx = [width[0][0],peak,width[0][1]]
            self.idx_peak = idx 
            #--------------------
            # Sub sample of just the peak
            rvs  = RV_sample[idx[0]:idx[-1]]
            corrs= (self.Correlation[idx[0]:idx[-1]])

            # Guess initial parameters for functions fit
            p0 = [RV, rvs[-1]-rvs[0], 1]  # Center , width , height 

        except IndexError:
            rvs  = RV_sample
            corrs= self.Correlation

            # Guess initial parameters for functions fit
            p0 = [RV, rvs[-1]-rvs[0], 1]  # Center , width , height 
        
        # Attemp to fit functions to get errors of the peak center

        # Fixed will be refer to make those function to have at the peak a value of 1 (which we enforce to be also our peak)
        try: # Fixed Gaussian
            popt, pcov = curve_fit(self.gaussian, rvs, corrs, p0=p0[:-1])
            self.params_F_gauss=[popt[0],np.sqrt(pcov[0,0]),popt, pcov]
            F_gauss=False
        except:
            F_gauss=True
            self.params_F_gauss=[np.nan]*4

        try: # Fixed Parabola
            popt, pcov = curve_fit(self.parabola, rvs, corrs, p0=p0[:-1])
            self.params_F_parab=[popt[0],np.sqrt(pcov[0,0]),popt, pcov]
            F_parab=False
        except:
            F_parab=True
            self.params_F_parab=[np.nan]*4

        try: # Parabola
            popt, pcov = curve_fit(self.parabola, rvs, corrs, p0=p0)
            self.params_parabola=[popt[0],np.sqrt(pcov[0,0]),popt, pcov]
            parab=False
        except:
            parab=True
            self.params_parabola=[np.nan]*4

        try: # Gaussian
            popt, pcov = curve_fit(self.gaussian,rvs,corrs, p0=p0)
            self.params_gauss=[popt[0],np.sqrt(pcov[0,0]),popt, pcov]
            gauss=False
        except:
            gauss=True
            self.params_gauss=[np.nan]*4
            
        
        the_fits=[F_gauss,F_parab,parab,gauss]

        if sum(the_fits)!=0:
            print('WARNING: Some of the fits didn\'t converge')
        

        RV,ERR,M_cent,std_cent,M_errs=self.get_errs_info(ALL_info=True)#,PRINT=True)

        self.RV_error        = ERR
        self.RV_median       = M_cent
        self.RV_std          = std_cent
        self.median_err_fits = M_errs
        # Return the Rv sampled, the correlation for each one of the Rv sampled, the best RV and the Err
        return RV_sample,self.Correlation,RV,ERR



    #######################################################################
    
    # Function to plot CCF and Template spectra and Obs spectra shifted
    def PLOT(self, RV        = None,             # RV to be diiplayed, if not given default result from CCF
                   window    = None,             # Range of RV displayed in RV vs Corr plot
                   window1   = None,             # Range of Wavelenght displayed in plot spectra right upper
                   window2   = None,             # Range of Wavelenght displayed in plot spectra right middle
                   window3   = None,             # Range of Wavelenght displayed in plot spectra right lower
                   zoom      = True,             # If you would like to include the lines that indicates were is the zoom 
                   PAUSE     = False,            # If we want to display quickly the figure (whe we call repeated times this function)
                   time      = 1,                # in case of PAUSE=True indicate the time that you want the figure to be displayed (seconds)
                   Save      = False,            # If you want to save the figure
                   Save_name ='Output_CCF.png',  # Name of figure in case of being Saved
            ) -> None:                           # This function does not return anything
        """
        # Function to plot CCF and Template spectra and Obs spectra shifted
    
        Arguments of Function:data type = Default

        RV       :float       = None,            # RV to be diiplayed, if not given default result from CCF
        window   :tuple[float]= None,            # Range of RV displayed in RV vs Corr plot
        window1  :tuple[float]= None,            # Range of Wavelenght displayed in plot spectra right upper
        window2  :tuple[float]= None,            # Range of Wavelenght displayed in plot spectra right middle
        window3  :tuple[float]= None,            # Range of Wavelenght displayed in plot spectra right lower
        zoom     :bool        = True,            # If you would like to include the lines that indicates were is the zoom 
        PAUSE    :bool        = False,           # If we want to display quickly the figure (whe we call repeated times this function)
        time     :float       = 1,               # in case of PAUSE=True indicate the time that you want the figure to be displayed (seconds)
        Save     :bool        = False,           # If you want to save the figure
        Save_name:str         ='Output_CCF.png', # Name of figure in case of being Saved
        """
        
        # Data stored
        template_wavelength  = self.t_wave
        template_flux        = self.t_flux
        observed_wavelength  = self.obs_wave
        observed_flux        = self.obs_flux
        RV_sample            = self.RV_sample
        Correlation          = self.Correlation
        best_RV              = self.RV
        idx                  = self.idx_peak

        N=len(observed_flux)
        # Initialize default values within the method if they are not provided
        window  = window  if window  is not None else [None]
        window1 = window1 if window1 is not None else [None]
        window2 = window2 if window2 is not None else [observed_wavelength[int(5*N/21)],observed_wavelength[int(8*N/21)]]
        window3 = window3 if window3 is not None else [observed_wavelength[int(13*N/21)],observed_wavelength[int(16*N/21)]]   

        if RV is None:
            RV   = best_RV
            TITLE= f'CCF found best RV={RV:.3f}$\pm${self.RV_error:.3f} [km/s]'
        else:
            TITLE= f'CCF found best RV={RV:.3f} [km/s]'
        
        #------------------------------------------
        # Image distribution
        Fig= plt.figure(17,figsize=(13,8))
        gs = gridspec.GridSpec(
                5,7,                    # Number of axis y,x
                height_ratios = [2,0.2,1,0.05,1],        # relatives ratios of heigh
                width_ratios  = [0.05,1,1,0.15,1,1,0.2], # relatives ratio of with
                left  = 0.05,     # Space to the edge of left   from the nearest axis
                right = 0.94,     # Space to the edge of right  from the nearest axis
                bottom= 0.1,      # Space to the edge of bottom from the nearest axis
                top   = 0.95,     # Space to the edge of top    from the nearest axis
                wspace= 0.2,      # Space horizontal between each of the axis
                hspace= 0.2)      # Space vertical   between each of the axis

                              # y,x
        ax0=Fig.add_subplot(gs[0,:3])
        AX0=Fig.add_subplot(gs[0,4:-1])
        ax1=Fig.add_subplot(gs[2,:])
        ax2=Fig.add_subplot(gs[4,:3])
        ax3=Fig.add_subplot(gs[4,4:-1],sharey=ax2)
        #------------------------------------------

        
        # Plot RV sample vs Correlation 
        ax0.axvline(RV,color='gray',ls=':',lw=2)
        ax0.plot(RV_sample,Correlation,color='#384dff',lw=1.1,zorder=-1)
        ax0.scatter(RV_sample[idx[0]:idx[-1]],Correlation[idx[0]:idx[-1]],color='blue',marker='x',s=1)
        ax0.scatter(RV_sample[idx[1]],Correlation[idx[1]],color='red',marker='x',s=1)
        ax0.set_xlabel('RV [km/s]')
        ax0.set_ylabel('Correlation')
        ax0.set_title(TITLE, fontsize=14)

        # Plot fits to the peak for error determination
        self.plot_fit_peak(AX0)
        AX0.set_title(f'Median RV = {self.RV_median:.3f}')
        
        # Plot Template Spectrum
        ax1.plot(template_wavelength,template_flux,color='#5c6b9c',lw=1.2,label='Template')
        ax2.plot(template_wavelength,template_flux,color='#5c6b9c',lw=1.2,label='Template')
        ax3.plot(template_wavelength,template_flux,color='#5c6b9c',lw=1.2,label='Template')

        # Plot Shifted Observed Spectrum
        shift_obs_spectra=observed_wavelength/((1 + RV / 299792.458))

        ax1.plot(shift_obs_spectra,observed_flux,color='red',lw=0.6)
        ax2.plot(shift_obs_spectra,observed_flux,color='red',lw=0.6)
        ax3.plot(shift_obs_spectra,observed_flux,color='red',lw=0.6)

        
        # Set tthe windows that we wan to display
        ax0.set_xlim(window [0],window [-1])
        ax1.set_xlim(window1[0],window1[-1])
        ax2.set_xlim(window2[0],window2[-1])
        ax3.set_xlim(window3[0],window3[-1])

        # Aesthetics
        ax1.set_ylabel('Flux')
        ax2.set_ylabel('Flux')
        # ax2.set_xlabel('Wavelength')
        # ax3.set_xlabel('Wavelength')
        Fig.supxlabel('Wavelength')
        plt.setp(ax3.get_yticklabels(), visible=False)
        ax3.plot([],[],color='red',lw=2,label='Observed\nSpectrum') # Just to se better the line in the Legend
        # ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 0.975), ncol=2) # Legend in the title position
        ax3.legend(loc='center left', bbox_to_anchor=(1.001, 0.75), ncol=1)

        if zoom:
            self.INSET_ZOOM(ax0,AX0,lw=3,color='lightseagreen',c1=True,c2=True,c3=False,c4=False)
            self.INSET_ZOOM(ax1,ax2,lw=3,color='lightseagreen')
            self.INSET_ZOOM(ax1,ax3,lw=3,color='lightseagreen')
        
        if Save:
            plt.savefig(Save_name,dpi=150)
        elif PAUSE:
            plt.pause(time)
            plt.clf()
        else:
            plt.show()
    
    #############################################################################
    
    # Easy and repeateable way to show the fit performed
    def plot_fit_peak(self,AX0=None,RV=None):

        RV_sample            = self.RV_sample
        Correlation          = self.Correlation
        best_RV              = self.RV
        idx                  = self.idx_peak

        # if not ax is provide we generate anothe figure
        if AX0 is None:
            fig, AX0 = plt.subplots()

        # For posible plots of shifts
        if not RV is None:
            check=np.argmin(abs(RV_sample-RV))
            if check<idx[0]:
                idx[0]=check
            elif check>idx[-1]:
                idx[-1]=check            
        else:
            RV=best_RV
        #-----------------------------
        # Data
        ColorI = 'gray'
        ColorII= '#b9722f'
        rvs,corrs=RV_sample[idx[0]:idx[-1]],Correlation[idx[0]:idx[-1]] # The selected peak

        Mean = lines.Line2D([], [], color=ColorI, marker='|', linestyle='None',
                          markersize=10, markeredgewidth=1.5,label='Best RV')
        Median = lines.Line2D([], [], color=ColorII, marker='|', linestyle='None',
                                markersize=10, markeredgewidth=1.5,label='Median RV')

       

        AX0.plot(rvs,corrs,color='blue',lw=2,zorder=10,label='Measured\n Correlation\n  Peak')
        AX0.axvline(best_RV,color=ColorI,lw=1)
        AX0.axvline(self.RV_median,color=ColorII,lw=1)
        
        AX0.scatter(RV,1,color='red',marker='X',s=12,zorder=11)
        
        # Fixed Gaussian
        try:
            center,err,popt,cov=self.params_F_gauss
            AX0.plot(rvs,self.gaussian(rvs,popt[0],popt[1]),color='green',label='Gaussian\n fixed')
        except:
            pass

        # Fixed Parabola
        try:
            center,err,popt,cov=self.params_F_parab
            AX0.plot(rvs,self.parabola(rvs,popt[0],popt[1]),color='#ff7c00',label='Parabola\n Fixed ')
        except:
            pass
        
        # Parabola
        try:
            center,err,popt,cov=self.params_parabola
            AX0.plot(rvs,self.parabola(rvs,popt[0],popt[1],popt[2]),color='yellow',label='Parabola')
        except:
            pass

        # Gaussian
        try:
            center,err,popt,cov=self.params_gauss
            AX0.plot(rvs,self.gaussian(rvs,popt[0],popt[1],popt[2]),color='red',label='Gaussian')
        except:
            pass

        leg=AX0.legend()

        for legobj in leg.legend_handles:
            legobj.set_linewidth(2.0)

        AX0.legend(handles=[Mean,Median]+leg.legend_handles,loc='center left',bbox_to_anchor=(1.001,0.5),ncol=1,labelspacing=1.5)

        AX0.set_xlabel('RV [km/s]')
        

    #############################################################################
    
    # Easy way to generate conection of zoom on plots
    def INSET_ZOOM(self,AX1,AX2,color='fuchsia',lw=5,alpha=0.5,c1=True,c2=False,c3=True,c4=False):
        # this generates the box (controls the with and color) and the conectors
        rectpatch, connects=AX1.indicate_inset_zoom(AX2,linewidth=lw, edgecolor=color,alpha=alpha)

        AX2.spines[:].set_color(color) # Color of the axes in the zoomed area
        AX2.spines[:].set_linewidth(lw)   # Linewith of axes in the zoomed plot
        AX2.spines[:].set_alpha(alpha)

        # Select Corners to conect
        connects[0].set_visible(c1)
        connects[1].set_visible(c2)
        connects[2].set_visible(c3)
        connects[3].set_visible(c4)

        for line in connects: # Linewidth of the line conecting the zoomed area
            line.set_linewidth(lw)

    #############################################################################
    # # Function to plot a gif of the CCF
    # # You need to provide or IDX_RV_sample or RV_range 

    # def GIF(self,RV_RANGE      = None,             # RV range that we want to display on the gif
    #              IDX_RV_sample = None,             # IDX of RV sample that we would like to plot
    #              path          = './',             # path to save the gif (also a temp directory will be created here)
    #              name          = 'CCF.gif',        # Name of output Gif
    #              window        = None,             # Range of RV displayed in RV vs Corr plot
    #              window1       = None,             # Range of Wavelenght displayed in plot spectra right upper
    #              window2       = None,             # Range of Wavelenght displayed in plot spectra right middle
    #              window3       = None,             # Range of Wavelenght displayed in plot spectra right lower
    #              loop          = True,
    #              duration      = 100,):
       
    #     tempdir=path+'creating_gif_do_not_touch'  # Name of temporal directory that will contain images to create the gif

    #     os.system(f'mkdir {tempdir}')
        

    #     # List of image file paths
    #     image_files = [f'{tempdir}/{i}.png' for i in range(250)]

    #     # Load the images
    #     images = [Image.open(img) for img in image_files]

    #     # Save as GIF
    #     images[0].save('/home/ian/Desktop/output.gif', save_all=True, append_images=images[1:], loop=loop, duration=duration)

    #     os.system(f'rm -rf {path}{tempdir}')

    ###############################################################

    # Function used to find peak in a 1D array (Used to find the peak of CCF)
    def find_peaks(self,
                   Y,                   # Data where we want to find a peak (should be just values)
                   npeak     = 1,       # Quantity of peaks that we want to search
                   peak      = 'max',   # Kind of peak that we want to search ('max' or 'min')
                   tolerance = 60,      # Percentage of tolerance (relative to the max and min height in the sample) that will be considered as part of the peak
                   SHOW      = False    # Variable introduced just to show how code works
                    ):
        #-------
        npeak = int(npeak)
        N     = len(Y)-1
        #-------
        if not 0<tolerance<100:
            raise Exception('The tolerance should be in percentage (between 0 and 100)')
        #-------
        if peak=='max':
            fun  = np.argmax
            nfun = np.argmin
        elif peak=='min':
            fun  = np.argmin
            nfun = np.argmax
        else:
            raise Exception('The peak variable should be \"max\" or \"min\"')
        #-------
        LIM = ((Y[fun(Y)]-Y[nfun(Y)])*tolerance/100)  # Limit where we will consider points belongin to the same peak
        sign= LIM/abs(LIM)
        LIM+= Y[nfun(Y)]
        #-------
        # if npeak<=1 :
        #     return fun(Y)
        #-----------------------
        if SHOW:
            X=np.arange(len(Y))
            Y_o=Y.copy()

        peak_idx   = []
        width_peak = []

        for i in range(npeak):
            #-------
            if SHOW:
                plt.plot(X,Y)
            #-------
            if all(Y<=LIM):
                print(f'Found every peak at the level of tolerance defined {tolerance}%, you can decrease this number in order to found more peak')
                print(f'Stoping at peak {i+1}')
                break
            #-------
            # Find the peak
            peak=fun(Y)
            peak_idx.append(peak)
            #-------
            # We move around the peak till we reach the tolerance level or the edge
            L,R=0,0
            l,r=(peak-1)>0,(peak+1)<N

            while l or r:
                if l:
                    L+=1
                    aux=(Y[peak-L]-LIM)*sign
                    l=(abs(aux)==aux) and ((peak-L)>0)
                if r:
                    R+=1
                    aux=(Y[peak+R]-LIM)*sign
                    r=(abs(aux)==aux) and ((peak+R)<N)
  

            Y[peak-L:peak+R+1]= LIM
            width_peak.append([peak-L,peak+R+1])
            #-------

        peak_idx   = np.array(peak_idx)
        width_peak = np.array(width_peak)

        if SHOW:
            plt.plot(X,Y)
            plt.scatter(X[peak_idx],Y_o[peak_idx],s=50,marker='*',facecolor='gray',edgecolor='k',zorder=npeak+1)
            plt.plot(X,np.repeat(LIM,len(Y)))
            plt.show()


        return peak_idx,width_peak