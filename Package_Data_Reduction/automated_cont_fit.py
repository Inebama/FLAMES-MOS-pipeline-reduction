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

from typing import Iterable,Callable

#____________________________________________________________________________________________________________________________________
class automated_continuum:
    """"
    # Class that allow to search convergence on relative continuum based on distribution of points\n

    The main idea is to shift several times the spectra in order to get for each pixel a distribution of points that are representative of the nearby values. We measure a defined high percentile of this distribution (per pixel) and in this way we get a continuum guess, at this point we create a mask around this guess based on the dispersion of the normalized spectra and we cut the point that are farthen than the limit stablished. We iterate this untill we reach a convergence (i.e. the dispersion measures start to fllaten). With our final selection we use a spline cubic function just to smooth the guess (since came form statistical inference it will show small fluctuations)  

    ### Parameters\n
        #\n
        wave: np.array,  # Iterable containing the wave of the spectra\n
        flux: np.array,  # Iterable containing a flux of the spectra (same length as wave)\n
        name: str,       # The cont will be stored in a fits file, if exists already it will overwrite it \n
        subfolder: str = 'CONT_SUB/',   # subfolder where the new file will be stored\n
        path: str      = './',          # path where the file will be stored\n
        wave_range: float|tuple[float,float] = (-1,np.inf), # wavelenth range (in case you want to remove the extremes)\n
        percentile_cont: int|tuple[int,int]  = (84,93),     # percentile used to search the continuum, also a range can be set this will be increasing alongside with the iterations\n 
        percentile_range: tuple              = (75,95),     # range search convergence in number of cont points \n
        percentile_edges: tuple              = (65,98),     # Relaxed percentile for edges of the spectrum (they will converge slower thatn the rest of the spectra)
        iter_max: int                        = 10,          # Number of max iteration to try convergence on the continuum\n
        RUN: bool                            = True,        # by default when calling it it will trigger the automatic run calculation, however if for example the smoothing went wrong you can avoid the calculation and make just use of the already computed values\n
        wave_units: str                      = 'Angstrom',  # units of wavelength to add in fits file\n
        save_plots: bool                     = False        # If you want there are some plots relative to the convergence of this continuum that you can store\n
        primary_header: fits.Header|None     = None,        # If you want to keep an original primary header you can provide it and it will be added to the output fits
    .          
    """
    def __init__(self,
                wave            : np.array,                           # Iterable containing the wave of the spectra\n
                flux            : np.array,                           # Iterable containing a flux of the spectra (same length as wave)\n
                name            : str,                                # the cont will be stored in a fits file, it will overwrite it\n
                path            : str                  = '',          # path where the file will be stored\n
                subfolder       : str                  = 'CONT_SUB/', # subfolder where the new file will be stored\n
                wave_range      : float|tuple          = (-1,np.inf), # wavelenth range (in case you want to remove the extremes), if you provide a number between 0-1 it will cut the proportional edge %\n
                percentile_cont : int|tuple[int,int]   = (84,92),     # percentile used to search the continuum, also a range can be set this will be increasing alongside with the iterations\n 
                percentile_range: tuple                = (75,93),     # range search convergence in number of cont points \n
                percentile_edges: tuple                = (65,95),     # Relaxed percentile for edges of the spectrum (they will converge slower thatn the rest of the spectra)\n
                iter_max        : int                  = 10,          # Number of max iteration to try convergence on the continuum\n
                RUN             : bool                 = True,        # by default when calling it it will trigger the automatic run calculation, however if for example the smoothing went wrong you can avoid the calculation and make just use of the already computed values\n
                wave_units      : str                  = 'Angstrom',  # units of wavelength to add in fits file\n    
                save_plots      : bool                 = False,       # If you want there are some plots relative to the convergence of this continuum that you can store\n
                primary_header  : fits.Header|None     = None,        # If you want to keep an original primary header you can provide it and it will be added to the output fits
       ): 
        # Save the values given
        if isinstance(wave_range,tuple):
            mask= (wave_range[0]<wave) & (wave<wave_range[1])
        elif 0<wave_range<1:
            EDGE=(np.nanmax(wave)-np.nanmin(wave))*wave_range
            mask=(np.nanmin(wave)+EDGE<wave) & (wave<np.nanmax(wave)-EDGE)
            print(r'Cut by % lead to cut of original points ',len(mask),'->',sum(mask))
        else:
            raise ValueError('The waverange provided does not meet the needs provided by the code')
        

        self.wave            = wave[mask]
        self.median_flux     = np.nanmedian(flux[mask])
        self.flux            = (flux/self.median_flux)[mask]
        self.path            = path
        self.name            = name
        self.wave_units      = wave_units
        self.percentile_range= percentile_range
        self.percentile_edges= percentile_edges
        self.sub_folder      = subfolder
        self.save_plots      = save_plots
        self.primary_header  = primary_header

        if isinstance(percentile_cont,tuple): 
            self.percentile_of_cont=np.linspace(min(percentile_cont),max(percentile_cont),iter_max).astype(int)  # We define the continuum percentile
        elif isinstance(percentile_cont,int):
            self.percentile_of_cont=np.repeat(percentile_cont,iter_max)

        if percentile_cont[-1] > max(percentile_range):
            raise ValueError('The percentile of the continuum can no be higher than the larger percentile used to determine the continuum')

        if iter_max<3:
            print('Iter max can not be less than 3, it will be forced to this value in this run')
            iter_max=3
        self.iter_max = iter_max

        self.cont_mask   = np.repeat(True,len(self.wave)) # We start using all datapoints
        self.cont_guess  = np.repeat(1,len(self.wave))
        self.current_flux= self.flux # we will work over this copy of teh flux
        self.npoints     = [len(self.wave)] # We will store the convergence of the number of points considered as part of continuum

        # Values to get median flux/Initial guess (migth be editted for you)
        self.seed    = 170499
        np.random.seed(self.seed) # Used to get replicable outputs
        self.shift   = 200        # Range in velocity of shift in km/s
        self.Nsamples= 2500       # Number of samples that will be created for spectra

        self.positive_shift_waves=self.wave<wave.min()*(1+3*self.shift/299792.458)
        self.negative_shift_waves=self.wave>wave.max()*(1-3*self.shift/299792.458)

        if RUN:
            self.RUN_app()

    #_______________________________________________________________________________________________________________________________
    # Funtion that looks for an initial guess of the continuum 
    # This is based in performing sligth shift an look for
    # Representative flux (high number percentile) to guide the eye when fitting (or directly use) and find initial mask of points that could be used to get the continnum with the classical sigma clipping
    def CONT_GUESS(self,iter):
        wave=self.wave

        RANDOM  = np.random.normal(0,self.shift,self.Nsamples) # in order to minimize impact on the extremes we precompute the random sampling

        wave=self.wave[self.cont_mask]
        flux=self.current_flux[self.cont_mask]
        FLUXES=[]
        # For the edges we will keep a conservative view of which initial value we will use for the first iter
        if iter==0 or True:
            sign   = np.sign(RANDOM).astype('int')+1   # we identify if the sampled RV is blue or red shifted

            # We estimate representative values of the edge of the spectrum
            M1 = np.median(self.current_flux[:int(sum(self.negative_shift_waves)/3)])
            M2 = np.median(self.current_flux[-int(sum(self.positive_shift_waves)/3):])
            
            # accordingly to the shift we compute a representative value of the current flux considered outside of the interpolation range
            edges= np.array([M1,(M1+M2)/2,M2])[sign] 

            # We sample the spectrum over shifts
            for rv,edge in zip(RANDOM,edges):
                temp_wave=wave/(1+rv/299792.458)
                temp_flux=interp1d(temp_wave, flux, bounds_error=False, fill_value=edge)(self.wave)
                FLUXES.append(temp_flux)

        # On the second iteration we will keep the first continuum guess on the edges
        cont=np.nanpercentile(FLUXES,self.percentile_of_cont[iter],axis=0,method='median_unbiased') # 1 sigma
   
        # Re-normalize and look for point in the vicinity of this continuum
        NORM=self.current_flux/cont
        a=np.nanpercentile(NORM,self.percentile_range[1])
        b=np.nanpercentile(NORM,self.percentile_range[0])

        mask=(b<NORM)&(NORM<a)

        # Again we are a bit more conservative on the edges since whe we shift they do not have the same informations the center wavelengths
        # Unless we are in the final iteration
        if iter!=self.iter_max-1:
            a=np.nanpercentile(NORM[self.negative_shift_waves],self.percentile_edges[1])
            b=np.nanpercentile(NORM[self.negative_shift_waves],self.percentile_edges[0])

            nshift= self.negative_shift_waves & (b<NORM)&(NORM<a)

            a=np.nanpercentile(NORM[self.positive_shift_waves],self.percentile_edges[1])
            b=np.nanpercentile(NORM[self.positive_shift_waves],self.percentile_edges[0])

            pshift= self.positive_shift_waves & (b<NORM)&(NORM<a)

            mask= mask | nshift | pshift

            # Store the number of points considered continuum
            self.npoints.append(sum(mask))

        if iter==0:
            
            return mask,cont
        else:
            self.current_flux=NORM
            return mask,cont

    #_______________________________________________________________________________________________________________________________ 
    def search_convergence(self,plot=False,save=False,key='',percentage_of_convergence=0.005):
        for i in range(self.iter_max):
            self.cont_mask,self.cont_guess=self.CONT_GUESS(i)
            if plot and i<self.iter_max-1:
                plt.figure(98765)
                plt.plot(self.wave,self.cont_guess,label=str(i+1),zorder=i+1)
                plt.figure(98764)
                plt.scatter(self.wave[self.cont_mask],self.flux[self.cont_mask],label=str(i+1),marker='o',s=(self.iter_max-i)*5) 
            if len(self.npoints)>3:
                if  (self.npoints[-3] - self.npoints[-2])/self.npoints[0]<percentage_of_convergence and  (self.npoints[-2] - self.npoints[-1])/self.npoints[0]<percentage_of_convergence:
                    self.cont_mask,self.cont_guess=self.CONT_GUESS(self.iter_max-1)
                    break
        if plot:
            plt.figure(98765)
            plt.plot(self.wave,self.flux,lw=0.5,color='k',zorder=0)
            plt.plot(self.wave,self.cont_guess,label='final',zorder=self.iter_max+1,color='cyan')
            plt.title('Cont Guess (just statistical)')
            plt.ylabel('relative flux')
            plt.xlabel('wavelength')
            plt.legend()
            plt.tight_layout()
            plt.figure(98764)
            plt.plot(self.wave,self.flux,lw=0.5,color='k',zorder=0)
            plt.scatter(self.wave[self.cont_mask],self.flux[self.cont_mask],label='final',color='cyan',marker='+') 
            plt.title('Cont points')
            plt.ylabel('relative flux')
            plt.xlabel('wavelength')
            plt.legend()
            plt.tight_layout()
            plt.figure(98763)
            plt.bar(np.arange(len(self.npoints)),self.npoints)
            plt.title('Convergence of continuum points')
            plt.ylabel('N points')
            plt.xlabel('N iterations')
            plt.tight_layout()
            if save:
                plt.figure(98765)
                plt.savefig(self.path+key+self.name+'_continuum.png',dpi=150)
                plt.figure(98764)
                plt.savefig(self.path+key+self.name+'_points_cont.png',dpi=150)
                plt.figure(98763)
                plt.savefig(self.path+key+self.name+'_converge.png',dpi=150)
            else:
                plt.show()


    #_______________________________________________________________________________________________________________________________ 
    # Final smoothing of the continuum converged
    def final_cont(self,smoothing_factor=0.04):
        self.cont_guess=self.flux/self.current_flux
        spline = UnivariateSpline(self.wave, self.cont_guess, s=smoothing_factor)
        return spline(self.wave)
    
    def RUN_app(self):
        self.search_convergence(key=self.sub_folder,plot=self.save_plots,save=self.save_plots)
        os.makedirs(self.path+self.sub_folder,exist_ok=True)
        self.UPDATE_fits(key=self.sub_folder)

    #_______________________________________________________________________________________________________________________________ 
    # Storage the Spectra with the continuum substracted
    # Create a new fits
    def UPDATE_fits(self,key=''):
        path=self.path
        name=self.name
        flux=self.flux
        wave=self.wave

        # If we provide a header fromold file we can keep it
        if not self.primary_header is None:
            header=self.primary_header
        else:
            # Create a FITS header
            header = fits.Header()
        
        # We get the continuum obtained with the parameters fitted (stored)
        continuum=self.final_cont()

        Nflux=flux/continuum

        phdu=fits.PrimaryHDU(header=header)
        
        # Create the columns of data
        columns=[
            fits.Column(name='wave', format='D', unit=self.wave_units,array= wave ),
            fits.Column(name='nflux', format='D', unit='Normalized flux',array= Nflux ),
            fits.Column(name='cont', format='D', unit='Relative flux',array= continuum ),
            fits.Column(name='flux', format='D', unit='Relative flux',array= flux ),
            fits.Column(name='Cont_guess',format='D', unit='Relative Flux',array= self.cont_guess),
            fits.Column(name='init_mask', format='L', unit='BOOL',         array= self.cont_mask ),
            ]
        # Generate the secondary fits with the data
        NEWFITS = fits.BinTableHDU.from_columns(columns)

        # Set the EXTNAME keyword to label the binary table HDU
        NEWFITS.header['EXTNAME'] = 'SPECTRUM (cont sub)'

        NEWFITS.header['medext']  = (self.median_flux,'Median flux divided at the beggining')
        # Information about how it was performed the init guess
        NEWFITS.header['Nsamp']  = (self.Nsamples,'Number of samples used to get inti guess')
        NEWFITS.header['Seed']   = (self.seed,'Seed used in numpy in order to get same values')
        NEWFITS.header['Shift']  = (self.shift,'Characteristical shift used to get init guess')
        NEWFITS.header['wMin']   = (self.wave.min(),'Min wave used')
        NEWFITS.header['wMax']   = (self.wave.max(),'Max wave used')
        NEWFITS.header['PerC']   = (self.percentile_of_cont[-1],'Final Percentile used as representative of initguess')
        NEWFITS.header['PerL']   = (self.percentile_range[0],'Percentile lower range used to fit')
        NEWFITS.header['PerU']   = (self.percentile_range[1],'Percentile upper range used to fit')
                

                
        # Comments about the columns (Topcat)
        NEWFITS.header['TCOMM1']  = 'wavelength '+self.wave_units
        NEWFITS.header['TCOMM2']  = 'Flux normalized' 
        NEWFITS.header['TCOMM3']  = 'Continnum used'
        NEWFITS.header['TCOMM4']  = 'Original Flux' 
        NEWFITS.header['TCOMM5']  = 'Statistical guess of cont' 
        NEWFITS.header['TCOMM6']  = 'Mask of points assumed as cont' 
        # Save the fits file
        hdul = fits.HDUList([phdu, NEWFITS])
        hdul.writeto(path+key+name, overwrite=True)
    


######################################################################################################################

# Function to read the fits files that contains the spectra given the name, feel free to update in case of need
def default_read(full_path,table=1,return_header=True):
    with fits.open(full_path) as hdul:
        header=hdul[0].header.copy()
        try:# Multidimensional data format
            wave=hdul[table].data['wave']
            flux=hdul[table].data['flux']
            # We scale the values of teh spectra around one for easier computation
        except:# 1D format standard IRAF
            flux =hdul[0].data
            CRVAL=hdul[0].header['CRVAL1']
            CDELT=hdul[0].header['CDELT1']
            wave =CRVAL+np.arange(len(flux))*CDELT
    if return_header:
        return wave,flux,header
    return wave,flux

# Function to obtain severl continuum of different spectra stored in the same subfolder just by giving the same inputs 
def multiple_spectra_cont(spectra         :Iterable,                             # Iterable of name files containing spectra
                          path            :str                     = './',        # path where the file will be stored
                          wave_range      :float|tuple[float,float]= (-1,np.inf), # wavelenth range (in case you want to remove the extremes)\n
                          percentile_cont :int|tuple[int,int]      = (84,93),     # percentile used to search the continuum, also a range can be set this will be increasing alongside with the iterations (recommended)
                          percentile_range:tuple                   = (75,95),     # range search convergence in number of cont points \n
                          percentile_edges:tuple                   = (65,98),     # Relaxed percentile for edges of the spectrum (they will converge slower thatn the rest of the spectra)
                          iter_max        :int                     = 10,          # Number of max iteration to try convergence on the continuum\n
                          wave_units      :str                     = 'Angstrom',  # units of wavelength to add in fits fileREAD_FUNCTION                        = None,        # We can give a proper function that reads our specific files, if None we will use a built in
                          READ_FUNCTION   :Callable|None           = None,        # Function to read you type of fits
                          subfolder       :str                     = 'CONT_SUB/'
                 ):

        if READ_FUNCTION is None:
            READ_FUNCTION=default_read

        for file in spectra:
            print('---\nGetting cont of file '+file)
            wave,flux,header=READ_FUNCTION(path+file)
            automated_continuum(wave            = wave,
                                flux            = flux,
                                name            = 'CS-'+file,              
                                path            = path,          
                                wave_range      = wave_range, 
                                percentile_cont = percentile_cont, 
                                percentile_range= percentile_range,    
                                percentile_edges= percentile_edges,    
                                iter_max        = iter_max,          
                                wave_units      = wave_units, 
                                subfolder       = subfolder,
                                primary_header  = header,
            )
    
    

if __name__=='__main__':
    from astropy.io import fits
    path='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_median_48h_new/OB-1/RESTFRAME/Continuum_sub/'
    target='MED_RGB_64=5.981916666666565=-71.95788888888768=2018-10-19T14:19:29.020.fits'
    path='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_median_48h_new/OB-0/RESTFRAME/Continuum_sub/'
    target='47TUC_R00198=5.403541666666575=-72.0786111111099=2011-11-26T01:12:09.821.fits'
    path='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_median_48h_new/OB-8/RESTFRAME/Continuum_sub/'
    target='MED_RGB_107=5.994791666666566=-71.90288888888769=2019-09-14T10:01:57.799.fits'
    with fits.open(path+target) as hdul:
        header=hdul[0].header.copy()
        wave=hdul[-1].data.wave
        flux=hdul[-1].data.flux
        Nflux=hdul[-1].data.Nflux
    
    if 1:
        automated_continuum(wave,flux,'Testing.fits',save_plots=True,wave_range=0.05,primary_header=header)
    else:
        path='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_median_48h_new/OB-8/RESTFRAME/'
        spectrums=np.genfromtxt(path+'Oficial_list.list',dtype='U100')
        
        multiple_spectra_cont(spectra    = spectrums, # Iterable of name files containing spectra
                              path       = path,      # path where the spectrums are allocated
                              wave_range = 0.05,      # wavelenth range (in case you want to remove the extremes)\n
                            )