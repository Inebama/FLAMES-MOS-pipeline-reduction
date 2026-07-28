#---------------------------------------------------------------
# Packages
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d       

from astropy.io import fits

import os
import sys
# from typing import Iterable # Not gonna used due deprecation
#---------------------------------------------------------------

class cont_fit:
    """
    # Class defined to perform continuum fitting to spectral data (fits files only)\n
    
    ### Parameters\n
    #
    spectra        :Iterable[str] ,     # Iterable containing names of fits files containing spectra\n
    path           :str  ='./',         # path where the spectras are allocated\n
    RUN            :bool = True,        # if True it will automatically start with the continuum fit\n
    Continue       :bool = False,       # if given a parameters stored complete path will continue from the last spectrum fitted\n
    wave_range     :float|tuple= (-1,np.inf), # wavelenth range (in case you want to remove the extremes)\n
    sigma_range    :tuple= (0.1,5),     # allowed range to test the sigma clipping\n
    iter_max       :int  = 50,          # allowed max number of iterations for sigma clipping\n
    smooth_range   :tuple= (0,5),       # allowed range for smooth parameter UnivariateSpline\n
    iter_default   :int  = 3,           # default initial value of iterations\n
    sigma_L_default:float= 1.25,        # default initial value of sigma below the continuum\n
    sigma_U_default:float= 2.5,         # default initial value of sigma above the continuum\n
    smooth_default :float= 1,           # default initial value of smooth parameter\n
    INP_OUT_FILE   :str  = 'Stored_params.npy', # File name for saving and reading the parameters (in case of continue from previous run).\n
    READ_FUNCTION        = None,        # We can provide function that reads our specific files, if None we will use a built in\n
    wave_units     :str  ='Angstrom'    # just to add in headers\n
                                         
    """
    def __init__(self,
                 spectra        :object,             # Iterable of name files containing spectra
                 path           :str  ='./',         # path where the spectras are allocated
                 RUN            :bool = True,        # if True it will automatically start with the continuum fit
                 Continue       :bool = False,       # if given a parameters stored complete path will continue from the last spectrum fitted
                 wave_range     :tuple= (-1,np.inf), # wavelenth range (in case you want to remove the extremes)
                 sigma_range    :tuple= (0.1,5),     # allowed range to test the sigma clipping
                 iter_max       :int  = 50,          # allowed max number of iterations for sigma clipping
                 smooth_range   :tuple= (0,5),       # allowed range for smooth parameter UnivariateSpline
                 iter_default   :int  = 0,           # default initial value of iterations
                 sigma_L_default:float= 1.25,        # default initial value of sigma below the continuum
                 sigma_U_default:float= 2.5,         # default initial value of sigma above the continuum
                 smooth_default :float= 0.03,        # default initial value of smooth parameter
                 INP_OUT_FILE   :str  = 'Stored_params.npy', # File name for saving and reading the parameters (in case of continue from previous run).
                 READ_FUNCTION        = None,        # We can give a proper function that reads our specific files, if None we will use a built in
                 wave_units     :str  ='Angstrom'     # just to add in headers
                 ):
   
        
        # Save the values given
        self.spectra        = spectra
        self.path           = path
        self.sigma_range    = sigma_range
        self.iter_max       = iter_max
        self.smooth_range   = smooth_range
        self.iter_default   = iter_default
        self.sigma_L_default= sigma_L_default
        self.sigma_U_default= sigma_U_default
        self.smooth_default = smooth_default
        self.wave_units     = wave_units
        self.wave_range     = wave_range

        #-----------------------------
        # I will assume that if you gave somenthing in fuction you have tested that works properly
        if READ_FUNCTION is None:
            self.Function_read_fits=self.default_read
        else:
            self.Function_read_fits=READ_FUNCTION
        #-----------------------------
        
        # Values to get median flux/Initial guess
        self.seed    = 170499
        np.random.seed(self.seed) # Used to get replicable outputs
        self.shift   = 180
        self.Nsamples= 2500
        self.RANDOM  = np.random.normal(0,self.shift,self.Nsamples)
        self.sign    = np.sign(self.RANDOM).astype('int')+1
        self.percentile_of_cont=86
        self.percentile_range=(80,95)

        # Auxiliar variable when change to different spectra
        self.updating_spectrum=False

        try:
            if Continue:
                n=self.read(INP_OUT_FILE)
                # If the file was not found (example running this by first time) we will display a message informing that it was not found 
                if n is None:
                    if False:
                        ax = plt.figure(figsize=(12,4)).add_subplot(xticks=[], yticks=[])
                        # The first word, created with text().
                        text = ax.text(-.1, .5, 'It has not been found the stored parameters file in the current directory,\n starting a new continuum substraction in 5s', color="k",fontsize=20)
                        ax.spines[:].set_visible(False)
                        plt.pause(5)
                        plt.close()
                    else:
                        print('\n'+'='*60)
                        print(f'\n\nStarting a New continuum fitting since no file as specified on INP_OUT_FILE={INP_OUT_FILE} was found in this folder\n\n')
                        print('='*60+'\n')
                    Continue=False
                else:
                    self.RUN(n)
        
            if RUN and not Continue:
                # Initialize storage for each spectra (parameters and states)
                self.stored_params   = {name: None for name in spectra}
                self.RUN()
        
        except Exception as E:
            print(E)
            print('\nEither way we will try to store the data of the fitting that has been already made\n')

        finally: # No matters whats happens we will always execute this
            self.save(INP_OUT_FILE) # We will always store the values 

    #----------------------------------------------------------------------
    
    # Function to read the fits files that contains the spectra given the name, feel free to update in case of need
    def default_read(self,hdul):
        try:# Multidimensional data format
            wave=hdul[1].data['wave']
            flux=hdul[1].data['flux']
            # We scale the values of teh spectra around one for easier computation
        except:# 1D format standard IRAF
            flux =hdul[0].data
            CRVAL=hdul[0].header['CRVAL1']
            CDELT=hdul[0].header['CDELT1']
            wave =CRVAL+np.arange(len(flux))*CDELT
        return wave,flux
    
    # Main function that operates over the readed fits/ it will add an extension to the fits file in order to save time and get reproductable data
    def get_spectra(self,name):
        with fits.open(self.path+name) as hdul:    
            wave,flux= self.Function_read_fits(hdul)

            # Keep this lines
            self.current_spectrum= name
            
            wave_range=self.wave_range
            if isinstance(wave_range,tuple):
                mask= (wave_range[0]<wave) & (wave<wave_range[1])
            elif 0<wave_range<1:
                EDGE=(np.nanmax(wave)-np.nanmin(wave))*wave_range
                mask=(np.nanmin(wave)+EDGE<wave) & (wave<np.nanmax(wave)-EDGE)
            flux      = flux[mask]
            self.wave = wave[mask]
            self.flux = flux/np.median(flux)

            try:
                self.init_mask =hdul[-1].data['init_mask']
                self.cont_guess=hdul[-1].data['Cont_guess']
            except:
                self.CONT_GUESS()
    
                # Generate new extension to not compute this every time
                columns=[
                    fits.Column(name='wave',      format='D', unit=self.wave_units,array= self.wave      ),
                    fits.Column(name='init_mask', format='L', unit='BOOL',         array= self.init_mask ),
                    fits.Column(name='Cont_guess',format='D', unit='Relative Flux',array= self.cont_guess)
                ]
                NEWFITS = fits.BinTableHDU.from_columns(columns)
                # Set the EXTNAME keyword to label the binary table HDU
                NEWFITS.header['EXTNAME'] = 'CONTINUUM GUESS'
                # Information about how it was performed the init guess
                NEWFITS.header['Nsamp']   = (self.Nsamples,'Number of samples used to get inti guess')
                NEWFITS.header['Seed']    = (self.seed,'Seed used in numpy in order to get same values')
                NEWFITS.header['Shift']   = (self.shift,'Characteristical shift used to get init guess')
                NEWFITS.header['wMin']    = (self.wave.min(),'Min wave used')
                NEWFITS.header['wMax']    = (self.wave.max(),'Max wave used')
                NEWFITS.header['PerC']    = (self.percentile_of_cont,'Percentile used as representative of initguess')
                NEWFITS.header['PerL']   = (self.percentile_range[0],'Percentile lower range used to fit')
                NEWFITS.header['PerU']   = (self.percentile_range[1],'Percentile upper range used to fit')
                
                hdul.append(NEWFITS)
                hdul.writeto(self.path+name, overwrite=True)

    
    # Funtion that looks for an initial guess of the continuum 
    # This is based in performing sligth shift an look for
    # Representative flux (high number percentile) to guide the eye when fitting (or directly use) and find initial mask of points that could be used to get the continnum with the classical sigma clipping
    def CONT_GUESS(self):
        wave=self.wave
        flux=self.flux
        
        M1 = np.median(flux[wave<wave.min()*(1+5*self.shift/299792.458)])
        M2 = np.median(flux[wave>wave.max()*(1-5*self.shift/299792.458)])
        AUX= np.array([M1,(M1+M2)/2,M2])
        
        FLUXES=[]
        for rv,edge in zip(self.RANDOM,AUX[self.sign]):
            temp_wave=wave/(1+rv/299792.458)
            temp_flux=interp1d(temp_wave, flux, bounds_error=False, fill_value=edge)(self.wave)

            FLUXES.append(temp_flux)

        cont=np.percentile(FLUXES,self.percentile_of_cont,axis=0,method='median_unbiased') # 1 sigma

        NORM=flux/cont

        a=np.nanpercentile(NORM,self.percentile_range[1])
        b=np.nanpercentile(NORM,self.percentile_range[0])

        mask=(b<NORM)&(NORM<=a)

        self.init_mask=mask
        self.cont_guess=cont
    
    #----------------------------------------------------------------------
    
    # Functions that defines the frame of the figure the interactive actions and the first plot
    def Figure_setup(self,N=0):
        # We create the figure where we will be seeing the spectra and the continuum
        # We create figure and axes using GridSpec
        fig = plt.figure(figsize=(14, 8))
        gs  = gridspec.GridSpec(6, 8,
                            height_ratios=[1,0.3, 0.15, 0.15, 0.15, 0.15], 
                            width_ratios =[0.1,0.1,0.1,0.2, 1, 0.15,0.05,0.05],
                            left   = 0.05, 
                            right  = 0.97, 
                            bottom = 0.08, 
                            top    = 0.9, 
                            wspace = 0.0, 
                            hspace = 0.0)

        # Axes for spectrum plot, sliders, and spectra list
        ax_spectrum = fig.add_subplot(gs[0, 3:])   # Large plot for spectrum and continuum
        ax_slider1 = fig.add_subplot(gs[2, 4])      # Sigma Min slider
        ax_slider2 = fig.add_subplot(gs[3, 4])      # Sigma Max slider
        ax_slider3 = fig.add_subplot(gs[4, 4])      # Iteration slider
        ax_slider0 = fig.add_subplot(gs[5, 4])      # Smothing slider

        # Scrollable list of spectra
        ax_list = fig.add_subplot(gs[0:4, 0:2])        
        ax_prev = fig.add_subplot(gs[5, 0])
        ax_next = fig.add_subplot(gs[5, 1])

        # Buttons for incrementing and decrementing sliders
        ax_inc1 = fig.add_subplot(gs[2, 7])
        ax_dec1 = fig.add_subplot(gs[2, 6])
        ax_inc2 = fig.add_subplot(gs[3, 7])
        ax_dec2 = fig.add_subplot(gs[3, 6])
        ax_inc3 = fig.add_subplot(gs[4, 7])
        ax_dec3 = fig.add_subplot(gs[4, 6])
        ax_inc4 = fig.add_subplot(gs[5, 7])
        ax_dec4 = fig.add_subplot(gs[5, 6])


   
        # We start using the first spectrum
        self.get_spectra(self.spectra[N]) 

        # display of intial plot and labels
        self.line_flux,      = ax_spectrum.plot(self.wave, self.flux,label='Flux')
        self.line_cont_guess,= ax_spectrum.plot(self.wave, self.cont_guess,color='lime',label='Init Guess')
        self.line_continuum, = ax_spectrum.plot(self.wave, self.flux, color='red' ,label='Continuum')

        ax_spectrum.set_xlabel('Wavelength')
        ax_spectrum.set_ylabel('Flux')
        ax_spectrum.set_xlim(min(self.wave),max(self.wave))
        ax_spectrum.set_ylim(0,np.percentile(self.flux,99)*1.2)
        ax_spectrum.legend(loc='lower right',bbox_to_anchor=(1.01,0.97),framealpha=0)

        # Set up sliders
        self.sigma_min_slider = Slider(ax_slider1, 'Sigma below'     , self.sigma_range[0] ,self.sigma_range[1] , valinit=self.sigma_L_default, valstep=0.01)
        self.sigma_max_slider = Slider(ax_slider2, 'Sigma above'     , self.sigma_range[0] ,self.sigma_range[1] , valinit=self.sigma_U_default, valstep=0.01)
        self.iter_slider      = Slider(ax_slider3, 'Iterations'      , -1                  ,self.iter_max       , valinit=self.iter_default   , valstep=1   )  # Int values only
        self.smother_slider   = Slider(ax_slider0, 'Smoother factor' , self.smooth_range[0],self.smooth_range[1], valinit=self.smooth_default , valstep=0.01)

        # Set up buttons 
        # Buttons for each slider (increment and decrement)
        self.btn_inc1 = Button(ax_inc1, '+')
        self.btn_dec1 = Button(ax_dec1, '-')
        self.btn_inc2 = Button(ax_inc2, '+')
        self.btn_dec2 = Button(ax_dec2, '-')
        self.btn_inc3 = Button(ax_inc3, '+')
        self.btn_dec3 = Button(ax_dec3, '-')
        self.btn_inc4 = Button(ax_inc4, '+')
        self.btn_dec4 = Button(ax_dec4, '-')

        self.btn_inc1.on_clicked(lambda event: self.increment_slider(self.sigma_min_slider,  0.01))
        self.btn_dec1.on_clicked(lambda event: self.increment_slider(self.sigma_min_slider, -0.01))
        self.btn_inc2.on_clicked(lambda event: self.increment_slider(self.sigma_max_slider,  0.01))
        self.btn_dec2.on_clicked(lambda event: self.increment_slider(self.sigma_max_slider, -0.01))
        self.btn_inc3.on_clicked(lambda event: self.increment_slider(self.iter_slider,       1   ))
        self.btn_dec3.on_clicked(lambda event: self.increment_slider(self.iter_slider,      -1   ))
        self.btn_inc4.on_clicked(lambda event: self.increment_slider(self.smother_slider,     .01))
        self.btn_dec4.on_clicked(lambda event: self.increment_slider(self.smother_slider,    -.01))

        self.btn_prev = Button(ax_prev, 'Previous')
        self.btn_next = Button(ax_next, 'Next')

        # Spectra list display 
        fig.suptitle(self.spectra[0], fontsize=14)
        self.spectra_slider = Slider(ax_list, 'Current Spectra', 0, len(self.spectra)-1, valinit=0, valstep=1, orientation='vertical')
        ax_list.axis('off')

        # Always store the current parameters for each spectrum
        self.params = {'sigma_min': self.sigma_min_slider.val,
                       'sigma_max': self.sigma_max_slider.val, 
                       'iterations':self.iter_slider.val,
                       'smooth':    self.smother_slider.val}

        self.fig        =fig
        self.ax_spectrum=ax_spectrum


        # Connect sliders to the update function
        self.sigma_min_slider.on_changed(self.update)
        self.sigma_max_slider.on_changed(self.update)
        self.iter_slider.on_changed(self.update)
        self.smother_slider.on_changed(self.update)

        # Connect buttoms to changue spectrum
        self.btn_prev.on_clicked(self.prev_spectrum)
        self.btn_next.on_clicked(self.next_spectrum)

        self.spectra_slider.on_changed(self.change_spectrum)

    
    #----------------------------------------------------------------------
    
    # Increment/Decrement button logic (keep values within slider limits)
    def increment_slider(self,slider, amount):
        new_val = np.clip(slider.val + amount, slider.valmin, slider.valmax)
        slider.set_val(new_val)

    
    # Previous button callback
    def prev_spectrum(self,event):
        current_val = self.spectra_slider.val
        if current_val > 0:
            self.spectra_slider.set_val(current_val - 1)


    # Update the slider value when "Next" is clicked
    def next_spectrum(self,event):
        current_val = self.spectra_slider.val
        if current_val < len(self.spectra) - 1:
            self.spectra_slider.set_val(current_val + 1)



    #----------------------------------------------------------------------
    
    # Update plot and store values
    def update(self,val=None):
        if not self.updating_spectrum:
            params=self.params
            params['sigma_min']  = self.sigma_min_slider.val
            params['sigma_max']  = self.sigma_max_slider.val
            params['iterations'] = self.iter_slider.val
            params['smooth']     = self.smother_slider.val
            
            # Recompute the continuum with updated values
            continuum = self.CONTINUUM(self.wave,self.flux,
                                sigma_clip      =(params['sigma_min'],params['sigma_max']),
                                max_iter        = params['iterations'],
                                smoothing_factor= params['smooth'])
            
            self.line_continuum.set_ydata(continuum)
            
            # Store the parameters for the current spectrum
            self.stored_params[self.current_spectrum] = [params['sigma_min'],  # 0
                                                         params['sigma_max'],  # 1
                                                         params['iterations'], # 2
                                                         params['smooth']]     # 3
            self.params=params
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()


    #----------------------------------------------------------------------
    
    # Function to change spectrum
    def change_spectrum(self,val):
        spectrum_name=self.spectra[int(val)]
        # Update current spectrum
        self.get_spectra(spectrum_name) # this updates the self.current_spectra/wave/flux
        # Change title
        self.fig.suptitle(spectrum_name, fontsize=14)
        # Update plot with the new spectrum
        self.line_flux.set_xdata(self.wave)
        self.line_flux.set_ydata(self.flux)
        self.ax_spectrum.set_xlim(min(self.wave),max(self.wave))
        self.ax_spectrum.set_ylim(0,np.percentile(self.flux,99)*1.2)

        self.line_cont_guess.set_xdata(self.wave)
        self.line_cont_guess.set_ydata(self.cont_guess)
        
        # Update wavelength of continumm in plot
        self.line_continuum.set_xdata(self.wave)


        # Reset sliders based on the new spectrum's stored parameters (if available)
        self.updating_spectrum= True

        if self.stored_params[spectrum_name]:
            self.sigma_min_slider.set_val(self.stored_params[spectrum_name][0])
            self.sigma_max_slider.set_val(self.stored_params[spectrum_name][1])
            self.iter_slider.set_val(     self.stored_params[spectrum_name][2])
            self.smother_slider.set_val(  self.stored_params[spectrum_name][3])
        else:
            self.sigma_min_slider.set_val(self.sigma_L_default) # Default value if not processed
            self.sigma_max_slider.set_val(self.sigma_U_default) # Default value if not processed
            self.iter_slider.set_val(     self.iter_default   ) # Default value for iterations
            self.smother_slider.set_val(  self.smooth_default ) # Default for smooth
        
        self.updating_spectrum= False
        self.update()


    #----------------------------------------------------------------------
    
    # Function used to compute the continuum
    def CONTINUUM(self,wavelength, flux, sigma_clip=(1.25, 2.5), max_iter=15, smoothing_factor=1):
        # If you set the iterations to -1 you get the init guess
        if max_iter==-1:
            return self.cont_guess.copy()
        # If you set the iterations to 0 you get an smoothed version of the init guess
        if max_iter==0:
            spline = UnivariateSpline(wavelength, self.cont_guess, s=smoothing_factor)
            return spline(wavelength)
        # if you set the iterations to 1, you make an spline estimation with points around the init guess
        # for max_iter > 1 you perform sigma clipping
        mask      = self.init_mask.copy()
        for _ in range(max_iter):
            # we fit a line acording to the points in the mask and input parameter of smooth
            spline = UnivariateSpline(wavelength[mask], flux[mask], s=smoothing_factor)
            continuum = spline(wavelength)
            residuals = flux - continuum
            std_dev = np.std(residuals[mask])
            # We apply an asymetric mask in order to obtain just the upper points an not the contribution from spectral lines
            mask = (residuals <= sigma_clip[1] * std_dev) & (residuals >= -sigma_clip[0] * std_dev)
            del spline
        return continuum
    #----------------------------------------------------------------------
    
    # Function that stores the values used
    def save(self,Output='Stored_params.npy'):
        np.save(self.path+Output,self.stored_params)
    #----------------------------------------------------------------------
    
    # Function that read the values previously used
    def read(self,Input='Stored_params.npy'):
        try:
            self.stored_params=np.load(self.path+Input,allow_pickle=True).item()
        except:
            return None
        # Check that we are talking about the same spectras
        check=  all([key in self.spectra for key in self.stored_params]) and len(self.spectra)==len(self.stored_params)

        n=0
        if check:
            for key in self.spectra:
                if self.stored_params[key]:
                    n+=1
                else:
                    return n
            else:# If never ends we have alredy check all and we plot the last one checked
                return len(self.stored_params)-1
        else:
            print('\n\nThe spectras given does not match the spectras loaded (maybe not the same number or at least one extra/missing spectra), please check this before continue\n\n')
            sys.exit('If you want to start from zero just set Continue=False \n')
    #----------------------------------------------------------------------
    
    # General way of run this code 
    def RUN(self,n=0):
        # Create the figure
        self.Figure_setup(n)

        # We will start working with the next spectra given on the list
        self.current_spectrum= self.spectra[n]
        # We set the slider in the spectra that we will analyse
        self.spectra_slider.set_val(n) # This updates the figure
        plt.show() # To display the figure and be able to start intercat with it

    #----------------------------------------------------------------------
    
    # Storage the Spectra with the continuum substracted
    # Create a new fits
    def UPDATE_fits(self,path,name,key):
        with fits.open(self.path+name) as hdul:
            header=hdul[0].header.copy()
            wave,flux=self.Function_read_fits(hdul)
            # We apply the same range that we used to fit the continuum
            wave_range=self.wave_range
            if isinstance(wave_range,tuple):
                mask= (wave_range[0]<wave) & (wave<wave_range[1])
            elif 0<wave_range<1:
                EDGE=(np.nanmax(wave)-np.nanmin(wave))*wave_range
                mask=(np.nanmin(wave)+EDGE<wave) & (wave<np.nanmax(wave)-EDGE)

            flux = flux[mask]

            wave = wave[mask]
            flux = flux/np.median(flux)

            try:
                self.init_mask =hdul[-1].data['init_mask']
                self.cont_guess=hdul[-1].data['Cont_guess']
            except:
                # After reviewing all spectras this should not happend
                print('\nHey apparently we were not able to found the init guess, I will compute anyway but this should not be happening,\n so maybe you should check ',name,'\n')
                self.get_spectra(name)
        
        # We get the continuum obtained with the parameters fitted (stored)
        params    = self.stored_params[name]
        continuum = self.CONTINUUM(wave,flux,
                                   sigma_clip      =(params[0],params[1]),
                                   max_iter        = params[2],
                                   smoothing_factor= params[3])

        Nflux=flux/continuum


        phdu=fits.PrimaryHDU(header=header)

        try:
            CRVAL=header['CRVAL1']
            CDELT=header['CDELT1']
            phdu.data=Nflux
        except:
            pass
        
        # Create the columns of data
        columns=[
            fits.Column(name='wave', format='D', unit=self.wave_units,array= wave ),
            fits.Column(name='nflux', format='D', unit='Normalized flux',array= Nflux ),
            fits.Column(name='cont', format='D', unit='Relative flux',array= continuum ),
            fits.Column(name='flux', format='D', unit='Relative flux',array= flux ),
            fits.Column(name='Cont_guess',format='D', unit='Relative Flux',array= self.cont_guess),
            fits.Column(name='init_mask', format='L', unit='BOOL',         array= self.init_mask ),
            ]
        # Generate the secondary fits with the data
        NEWFITS = fits.BinTableHDU.from_columns(columns)

        # Set the EXTNAME keyword to label the binary table HDU
        NEWFITS.header['EXTNAME'] = 'SPECTRUM (cont sub)'

        # Final parameters used in the application
        NEWFITS.header['SIGMAL']= (params[0] ,'lower range for sigma cliping')
        NEWFITS.header['SIGMAU']= (params[1] ,'upper range of sigma cliping' )
        NEWFITS.header['ITERS'] = (params[2] ,'iterations of sigma cliping'  )
        NEWFITS.header['SMOOTH']= (params[3] ,'Smooth parameter spline3'     )
        # Information about how it was performed the init guess
        NEWFITS.header['Nsamp']  = (self.Nsamples,'Number of samples used to get inti guess')
        NEWFITS.header['Seed']   = (self.seed,'Seed used in numpy in order to get same values')
        NEWFITS.header['Shift']  = (self.shift,'Characteristical shift used to get init guess')
        NEWFITS.header['wMin']   = (self.wave.min(),'Min wave used')
        NEWFITS.header['wMax']   = (self.wave.max(),'Max wave used')
        NEWFITS.header['PerC']   = (self.percentile_of_cont,'Percentile used as representative of initguess')
        NEWFITS.header['PerL']   = (self.percentile_range[0],'Percentile lower range used to fit')
        NEWFITS.header['PerU']   = (self.percentile_range[1],'Percentile upper range used to fit')
                
        # Comments about the columns (Topcat)
        NEWFITS.header['TCOMM1']  = 'wavelength '+self.wave_units
        NEWFITS.header['TCOMM2']  = 'Flux normalized' 
        NEWFITS.header['TCOMM3']  = 'Continnum used'
        NEWFITS.header['TCOMM4']  = 'Original Flux' 
        NEWFITS.header['TCOMM5']  = 'Statistical guess of cont' 
        NEWFITS.header['TCOMM6']  = 'Mask of nearby points assumed as cont' 
        # Save the fits file
        hdul = fits.HDUList([phdu, NEWFITS])
        hdul.writeto(path+key+name, overwrite=True)
    #----------------------------------------------------------------------
    
    # Store the fits files in a subfolder
    def store_cont_sub(self,key='Continuum_sub/',path=None,Force=False):
        # This will only work if we have end to fill the stored params
        check= all(map(lambda x: x is not None,self.stored_params.values()))
        if not check and not Force:
            sys.exit('\n\nThe Stored params are not complete please fill them up before storing in the new fits files\n\n')
            return None

        if path is None:
            path=self.path
        
        print(f'\n\n\nCreating the new fits with continuum substracted, please wait ...\n')
        # Create a SubFolder where we will place the spectrums with the continuum substracted
        try:
            os.makedirs(f'{path}{key}', exist_ok=True)
        except:
            pass
        
        for name in self.spectra:
            if not self.stored_params[name] is None: # Only those with fitting performed
                self.UPDATE_fits(path,name,key)
        
        print(f'\n\n\nAll new fits file are saved to the folder {path+key}\n\n\n')
            
            
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# if __name__=='__main__':

#     path='/home/ian/Desktop/testing_again_and_again/OB-1/'
#     # path='/home/ian/Desktop/OUTPUT_DATA/TESTING-A-288-289-290/'


#     # Synthetic spectra
#     # path_syn_spec  = '/home/ian/Desktop/Test-Autokur/kidH_R37000.dat'

#     # syn_spec=np.genfromtxt(path_syn_spec,skip_header=47)

#     # Template_wave=syn_spec[:,0]
#     # Template_flux=syn_spec[:,1]


#     spectrums=np.genfromtxt(path+'Oficial_list.list',dtype='U100')


#     # CONT=cont_fit(spectrums,path,wave_range=(6391.4,6613))

#     CONT=cont_fit(spectrums,path,wave_range=(6391.4,6613),Continue=True)#.store_cont_sub()



if __name__=='__main__':

    path='/home/ian/Desktop/Test-Autokur/'



    spectrums=['kidJB_R47000.fits','UVES_3047=6.035874999999899=-72.06522222222101=2019-09-14T10:02:14.945.fits']


    # Function used to create a new fits file containing the spectra after being sky substrcated and merged (if applies)
    # def new_fits(Fname,
    #             OBJ,
    #             Flux,Wave,
    #             CRVAL= None,
    #             CDELT= None,
    #             MAG  = None,
    #             path ='./',
    #             wave_to_angstrom=1):
        
    #     # Create a FITS header
    #     header = fits.Header()

    #     # Info about the star
    #     header['OBJECT']  = (OBJ,   'Object name in parent file 1')
        


    #     # Comments Example
    #     header['COMMENT'] = "File done by I. Baeza"

    #     # Create primary hdu
    #     primary_hdu = fits.PrimaryHDU(header=header)
        
    #     # If we provide the values we also write fits file in iraf like 
    #     if not (CRVAL is None) and not (CDELT is None):
    #         header['CRVAL1']  = (CRVAL, 'wavelenght first pixel ')
    #         header['CDELT1']  = (CDELT, 'Step of wavelenght ')

    #         primary_hdu.data=Flux

    #     # Create the columns of data
    #     columns=[
    #         fits.Column(name='wavelength', format='D', unit='Angstrom',  array= Wave*wave_to_angstrom ),
    #         fits.Column(name='Flux',       format='D', unit='adu', array= Flux )
    #         ]

    #     # Generate the secondary fits with the data
    #     NEWFITS = fits.BinTableHDU.from_columns(columns)

    #     # Set the EXTNAME keyword to label the binary table HDU
    #     NEWFITS.header['EXTNAME'] = 'SPECTRUM'

    #     # Save the fits file
    #     hdul = fits.HDUList([primary_hdu, NEWFITS])
    #     hdul.writeto(path+Fname, overwrite=True)
    
    # file=path+"kidJ_R47000.dat"
    # wave,flux,line,cont=np.genfromtxt(file,skip_header=41,unpack=True)

    # new_fits(Fname=spectrums[0],OBJ='synth_kidJ',Wave=wave,Flux=line,path=path)

    CONT=cont_fit(spectrums,path,Continue=True).store_cont_sub()




# PROBLEMS

#---------------------
# OBS 3
# MED RGB 370
# MED RGB 402

#---------------------
# OBS 6
# Descartada

#---------------------
# OBS
