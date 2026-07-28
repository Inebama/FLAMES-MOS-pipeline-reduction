# Packages 
import sys # Progress bar
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from tkinter import filedialog as fd
from prettytable import PrettyTable
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, RangeSlider
import multiprocessing as mp

# Just for typing annotation
try:
    from typing import Iterable
except:
    from collections.abc import Iterable



#___________________________________________________________________________

# I recommend to not touch this if you are unsure of what you are doing
def log_likelihood_cmd( delta_mag, 
                        delta_col, 
                        inv_cov_matrix,
                        N=2,
                        )-> Iterable:
    """
    Compute the log-likelihood for a set of observed points compared to model points
    using inverse covariance matrices. 

    Parameters:
    - delta_mag: array, delta ovserved-modeled magnitudes
    - delta_col: array, delta observed-modeled colors
    - inv_cov_matrix: array of shape (2,2)
    - N: Degree inside the exponential (higher values steep decay implies more inportante distance than density, N=2 Gaussian like)

    Returns:
    - normalized weigths of points: array of exp10(log-likelihoods) (one per point of the model)
    """
    
    delta = np.stack([delta_mag, delta_col], axis=1)  # (n_points, 2)

    if inv_cov_matrix.ndim == 2:
        # Single covariance matrix applied to all points
        logL = -0.5* N * np.sum(delta @ inv_cov_matrix * delta, axis=1)
    else:
        raise IndexError("The functions does not take a n dimensional covariance matrix, if you want to provide individual values for each point please give it ordered per each pair of point an unique matrix")

    MAX     = logL.max()   # We will use this to store the min distance
    logL   -= MAX          # We set Max Value = 0
    weights = np.exp(logL) # Exponential decay of values contribution

    # mask= (weights <= 0)

    N_w = len(weights)#sum(~mask)       # We will store the number of model points used 

    # Contributions of less than 1% of the maximum weigth are uninformative and could (will) bias our results
    # weights[mask] = 0

    # Normalization
    weights/= np.sum(weights)

    return weights,(MAX/(-0.5*N))**0.5,N_w

#___________________________________________________________________________

def weighted_mean_and_std(variable, weights):
    """
    Computes the weighted mean and weighted standard deviation.
    
    Parameters:
        variable : ist or array of the values (e.g., mass, teff,logg)
        weights  : array-like Associated weights (likelihoods)
        
    Returns:
        mean : float, Weighted mean of the distribution.
        std  : float, Weighted standard deviation of the distribution.
    """
    mask= (weights>0)

    weights= weights[mask]
    variable= variable[mask]

    # Weighted mean
    mean = np.nansum(weights * variable) / np.nansum(weights)

    # # Asymetric errors
    X   = (variable - mean)


    mask     = (X>=0)
    variance = np.sum((weights[mask]*X[mask]**2)/np.sum(weights[mask]))
    std_u    = np.sqrt(variance)

    mask     = (X<=0)
    variance = np.sum((weights[mask]*X[mask]**2)/np.sum(weights[mask]))
    std_l    = np.sqrt(variance)

    return mean, -std_l,std_u
#___________________________________________________________________________
# Function used to determine the nearest point in the line of sight of a point around it (default in 360 degrees bins)
def line_of_sight_of_point(delta_y,delta_x,nbins=360):
    r     = np.sqrt(delta_x**2 + delta_y**2)
    theta = np.arctan2(delta_y, delta_x)  # en radianes, de -pi a pi

    # we just divide in degrees
    theta_bins = np.linspace(-np.pi, np.pi, nbins + 1)
    digitized  = np.digitize(theta, theta_bins)

    selected_indices = []
    for i in range(1, nbins + 1):  # bins están indexados desde 1
        idx_in_bin = np.where(digitized == i)[0]
        if len(idx_in_bin) == 0:
            continue
        idx_closest = idx_in_bin[np.argmin(r[idx_in_bin])]
        selected_indices.append(idx_closest)

    selected_indices = np.array(selected_indices)

    return selected_indices

#___________________________________________________________________________
# Function that helps to the split all data to be calculated in chunks per core
def process_chunk(chunk_indices):
    results = {}
    for i in chunk_indices:
        results[i] = process_single_point(i)
    return results

# Function that globalizes the variables in order to meke them accesible from multiprocessing
def init_worker(_MAG, _COLOR, _N_iso, _iso_mag, _iso_col,
                _max_mag_dist, _max_color_dist, _min_mag_dist, _min_color_dist,
                _inferred_keys, _Iso_data, _Iso_mask,
                _check_inv, _inv_cov, _degree,_N_bins):

    global MAG, COLOR, N_iso, iso_mag, iso_col
    global max_mag_dist, max_color_dist, min_mag_dist, min_color_dist
    global inferred_keys, Iso_data, Iso_mask
    global check_inv, inv_cov, degree, N_bins

    MAG            = _MAG
    COLOR          = _COLOR
    N_iso          = _N_iso
    iso_mag        = _iso_mag
    iso_col        = _iso_col
    max_mag_dist   = _max_mag_dist
    max_color_dist = _max_color_dist
    min_mag_dist   = _min_mag_dist
    min_color_dist = _min_color_dist
    inferred_keys  = _inferred_keys
    Iso_data       = _Iso_data
    Iso_mask       = _Iso_mask
    check_inv      = _check_inv
    inv_cov        = _inv_cov
    degree         = _degree
    N_bins         = _N_bins


#___________________________________________________________________________
# Function that makes the main calculations for each point, requires global variables from init_worker
def process_single_point(i):
    result = {key: np.nan for key in inferred_keys}
    result.update({key+'_err_l': np.nan for key in Iso_data.dtype.names})
    result.update({key+'_err_u': np.nan for key in Iso_data.dtype.names})
    result['good']  = -1
    result['d_min'] = np.nan
    result['NP_mod']= np.nan

    # try:
    mag, color = MAG[i], COLOR[i]
    mask = ((iso_mag < mag + max_mag_dist) & (iso_mag > mag - max_mag_dist) &
            (iso_col < color + max_color_dist) & (iso_col > color - max_color_dist))

    if not np.any(mask):
        result['good'] = 0
        return i,result

    d_mag = iso_mag[mask] - mag
    d_col = iso_col[mask] - color

    if (abs(d_mag / min_mag_dist) + abs(d_col / min_color_dist)).min() > 1.3: # Similar area of elipse but less computationally expensive
        result['good'] = 0
        return i,result

    # ===  ===
    weights = np.zeros(N_iso)

    idx = line_of_sight_of_point(d_mag, d_col,N_bins)
    if check_inv:
        weights[np.where(mask)[0][idx]],MAX,N_w = log_likelihood_cmd(d_mag[idx], d_col[idx], inv_cov, degree)
    else:
        weights[np.where(mask)[0][idx]],MAX,N_w = log_likelihood_cmd(d_mag[idx], d_col[idx], inv_cov[i], degree)

    for key in Iso_data.dtype.names:
        mean, std_l, std_u = weighted_mean_and_std(Iso_data[key][Iso_mask], weights)
        result[key]          = mean
        result[key+'_err_l'] = std_l
        result[key+'_err_u'] = std_u
    
    result['d_min'] = MAX
    result['NP_mod']= N_w
    result['good']  = 1
    
    return i,result

    # except Exception:
    #     print(Exception)
    #     return i,result


class ISO_FIT:
    """
    Class thought to easyly fit parameters need on an Isochrone and retrieve atmospheric parameters
    Make sure your sample of Targets is a relatively clean sample 

    
    ### PARAMETERS : data type = default value , Meaning
    

        MAG       : Iterable ,                Magnitude of stars to be used
        COLOR     : Iterable ,                Color of stras to be used
        err_mag1  : float|Iterable = 0.05,    You could provide error for magnitude (global or individual per each source) err 5% on flux \\implies +2.5*\\log10(1.05) or -2.5\\log10(0.95) \\approx +/- 0.05 mag
        err_mag2  : float|Iterable = 0.05,    You could provide error for magnitude that compose color mag2-mag3 (global or individual)
        err_mag3  : float|Iterable = 0.05,    You could provide error for magnitude that compose color mag2-mag3 (global or individual)
        err_dist  : float|Iterable = 0.03,    Error in mag produced by error in distance => \\delta(m-M)\sim 2.17 \\times ((\\delta d)/d) \\implies 1%  0.027 mag
        err_redd  : float|Iterable = 0.05,    Error in visual absorption => \\delta A_v \\propto \\delta color \\implies \\delta E(B-V)= 0.015 \\approx delta A_v=0.05                
        iso_path  : str      | None = None,   If you want you can provide directly a path otherwise it will ask you to point one with an emergent window
        degree    : float           = 5.00,   Degree of decaiment on the relevance of point (higher degree priorizes nearest data, lower degree priorizes density of models see likelihood function)         
        OBJ       : Iterable | None = None,   [Optional] an iterable with some identifier to be add to the output file          
        RA        : Iterable | None = None,   [Optional] an iterable with Right Ascention of the stars to be add in the output file
        DEC       : Iterable | None = None,   [Optional] an iterable with the Declination of the stars to be add on the output file
        photometry: str  = 'Gaia',            gaia or ubvir are allowed by default (does not matter upper or lower case)
        Header    : int  = 13,                Number of columns of header in the isochrone (old version 13, PARSEC V2 ->14)
        Stages    : list = [                  List Of stages that you want to see and fit according to the PARSEC isochrones 'label' keyword 
                            0,   Pre Main Sequence - by default not used
                            1,   Main Sequence
                            2,   Sub Giant
                            3,   RGB
                            4,   CHEB
                            5,   Red-CHEB          - by default not used
                            6,   Blue-CHEB         - by default not used
                            7,   E-AGB
                            8,   TP-AGB
                            9,   Post AGB          - by default not used
                            ],
                 ):
    
    """

    def __init__(self,
                 MAG       : Iterable ,                # Magnitude of stars to be used
                 COLOR     : Iterable ,                # Color of stars to be used
                 err_mag1  : float|Iterable   = 0.05,  # You could provide error for magnitude (global or individual per each source)
                 err_mag2  : float|Iterable   = 0.05,  # You could provide error for magnitude that compose color mag2-mag3 (global or individual)
                 err_mag3  : float|Iterable   = 0.05,  # You could provide error for magnitude that compose color mag2-mag3 (global or individual)
                 err_dist  : float|Iterable   = 0.03,  # Error in mag produced by error in distance => \delta(m-M)\sim 2.17 \times ((\delta d)/d) \implies 1%  0.027 mag
                 err_redd  : float|Iterable   = 0.05,  # Error in visual absorption => \delta A_v \propto \delta color \implies \delta E(B-V)= 0.015 \approx delta A_v=0.05                
                 iso_path  : str      | None  = None,  # If you want you can provide directly a path otherwise it will ask you to point one with an emergent window
                 degree    : float            = 2.0,  # Degree of decaiment on the relevance of point (higher degree priorizes nearest data, lower degree priorizes density of models see likelihood function)
                 N_bins    : int              = 360,   # Cones to look around each point for the nearest model point (default search per 1 degree)
                 OBJ       : Iterable | None  = None,  # [Optional] an iterable with some identifier to be add to the output file          
                 RA        : Iterable | None  = None,  # [Optional] an iterable with Right Ascention of the stars to be add in the output file
                 DEC       : Iterable | None  = None,  # [Optional] an iterable with the Declination of the stars to be add on the output file
                 photometry: str  = 'Gaia',            # gaia or ubvir are allowed by default (does not matter upper or lower case)
                 Header    : int  = 13,                # Number of columns of header in the isochrone (old version 13, PARSEC V2 ->14)
                 Stages    : list = [                  # List Of stages that you want to see and fit according to the PARSEC isochrones 'label' keyword 
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
        mask = ()#(MAG<14)#&(MAG>11)
        nnn  = -1#200#11
        self.MAG      = MAG  [mask][:nnn]  # Target Magnitude
        self.COLOR    = COLOR[mask][:nnn]  # Target Color

        # Errors
        self.err_mag1=err_mag1      
        self.err_mag2=err_mag2 
        self.err_mag3=err_mag3        
        self.err_dist=err_dist       
        self.err_redd=err_redd      
        
        # Suggested extra info
        self.OBJ      = OBJ[mask][:nnn]      # Array of objject if given
        self.RA       = RA [mask][:nnn]      # rigth ascention to identify targets
        self.DEC      = DEC[mask][:nnn]      # declination to identify targets

        # Depends on the version of PARSEC isochrones
        self.Stages = Stages     # Evolutionary stages to be used
        self.HEADER = Header

        # Number of bins used to look for model points
        self.N_bins=N_bins

        # Degree of decaiment in contribution of the likelihood
        self.degree = degree

        # Stablishes the photometric system and their absorption coeficients
        self.photometric_system(photometry)
        
        # Run the application
        self.run_app(path=iso_path)

    #____________________________________________________________
    def run_app(self,path=None):
        try:
            plt.close()
        except:
            print('Not figure')
        # We read the isochrone and define the variables
        self.isochrone_read(path)

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
    # I recommend to not touch this if you are unsure of what you are doing
    def inverse_covariance_matrix_cmd(self,
                                    err_mag1:float|Iterable=0.05, 
                                    err_mag2:float|Iterable=0.05,
                                    err_mag3:float|Iterable=0.05,
                                    err_dist:float|Iterable=0.03,  # \delta(m-M)\sim 2.17 \times ((\delta d)/d) \implies 1%  0.027 mag
                                    err_redd:float|Iterable=0.05,  # \delta A_v \propto \delta color \implies \delta E(B-V)= 0.015 \approx delta A_v=0.05  
                                    )->Iterable:
        """
        ## Function that retrieves the inverse covariance matrix for the photometric system being used

        Parameters:
        - Default values are representative of ~10% of error in flux, ~1% error in distance and a direct variation of 0.05 in the visual Absorption 
        - All parameters can be float or iterable (must be consistent):
            err_mag1 : Error in mag1 (for CMD magnitude)
            err_mag2 : Error in mag2 (for CMD color asuming mag2-mag3)
            err_mag3 : Error in mag3 (for CMD color asuming mag2-mag3)
            err_dist : Error in distance modulus (affects magnitude only)    units:[mag]
            err_redd : Error in reddening (affects both magnitude and color) units:[Av]

        Returns:
            Iterable of inverse covariance matrices, one per source.
        """
        # A_mag1, A_mag2, A_mag3 : Extinction coefficients for the bands
        A_mag1=self.absorption_coeff[self.system['Mag1']]
        A_mag2=self.absorption_coeff[self.system['Mag2']]
        A_mag3=self.absorption_coeff[self.system['Mag3']]

        # I f we use the same magnitude for m1 and the color this will correlate the errors of both measures
        if A_mag1==A_mag2:
            factor= 1
        elif A_mag1==A_mag2:
            factor=-1
        else:
            factor= 0

        # Calculate variances
        var_mag = err_mag1**2 + err_dist**2 + (A_mag1 * err_redd)**2             # Affected by intrinsic error in measure, error in the distance estimation and by reddening
        var_col = err_mag2**2 + err_mag3**2 + ((A_mag2 - A_mag3) * err_redd)**2  # affected by intrinsic errro and reddening

        # Calculate covariance
        cov_mag_col = factor*err_mag1**2+ (A_mag1)*(A_mag2 - A_mag3)*(err_redd**2) # both proportionally to Av therefore 100% correlated (this error)

        # Now build covariance matrices
        try:
            cov_matrices = np.zeros((len(var_mag), 2, 2))
            cov_matrices[:, 0, 0] = var_mag
            cov_matrices[:, 0, 1] = cov_mag_col
            cov_matrices[:, 1, 0] = cov_mag_col
            cov_matrices[:, 1, 1] = var_col
        except:
            cov_matrices = np.zeros((2, 2))
            cov_matrices[0, 0] = var_mag
            cov_matrices[0, 1] = cov_mag_col
            cov_matrices[1, 0] = cov_mag_col
            cov_matrices[1, 1] = var_col

        # Invert each matrix
        inv_cov_matrices = np.linalg.inv(cov_matrices)

        return inv_cov_matrices

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
        ax_MH    = fig.add_subplot(gs[2, 6:8])        
        label_mh = fig.add_subplot(gs[4, 6:9])

        # Avail AGE 
        ax_AGE    = fig.add_subplot(gs[2, 9:11])        
        label_age = fig.add_subplot(gs[4, 9:11])

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
        self.iso_scatter  = ax_CMD.scatter(self.Fit_Color_iso,self.Fit_Mag_iso,c=self.Iso_mask,marker='+',s=10,vmin=0,vmax=1,cmap='bwr',rasterized=True)
        self.data_scatter = ax_CMD.scatter(self.COLOR, self.MAG,label='Targets',marker='o',color='gray',s=5,alpha=0.5,rasterized=True)
        
        ax_CMD.set_xlabel(f'{self.system["Mag2"]}-{self.system["Mag3"]}')
        ax_CMD.set_ylabel(f'{self.system["Mag1"]}')

        # just for display purposes
        ax_CMD.set_xlim(min(self.COLOR),max(self.COLOR))
        ax_CMD.set_ylim(max(self.MAG),min(self.MAG))

        ax_CMD.legend(loc='lower right',bbox_to_anchor=(1.01,0.97),framealpha=0)

        # Set up sliders
        self.AV_slider  = Slider(ax_slider1, 'Visual\nAbsorption', 0 , 3 , valinit=self.AV, valstep=0.01)
        self.UM_slider  = Slider(ax_slider2, 'Distance\nModulus' , 0 , 25, valinit=round(self.UM,2), valstep=0.01)
        self.MH_slider  = RangeSlider(ax_MH, '', valmin=0, valmax=len(self.MH) - 1, valinit=(0, len(self.MH) - 1), valstep=1,orientation='vertical')
        self.AGE_slider = RangeSlider(ax_AGE,'', valmin=0, valmax=len(self.AGE)- 1, valinit=(0, len(self.AGE)- 1), valstep=1,orientation='vertical')

        # Set up buttons 
        # Buttons for each slider (increment and decrement)
        self.btn_inc1   = Button(ax_inc1, '+')
        self.btn_dec1   = Button(ax_dec1, '-')
        self.btn_inc2   = Button(ax_inc2, '+')
        self.btn_dec2   = Button(ax_dec2, '-')    

        self.btn_inc1.on_clicked(lambda event: self.increment_slider(self.AV_slider,  0.01))
        self.btn_dec1.on_clicked(lambda event: self.increment_slider(self.AV_slider, -0.01))
        self.btn_inc2.on_clicked(lambda event: self.increment_slider(self.UM_slider,  0.02))
        self.btn_dec2.on_clicked(lambda event: self.increment_slider(self.UM_slider, -0.02))
        
        self.fig   =fig
        self.ax_CMD=ax_CMD

        # Connect sliders to the update function
        self.AV_slider.on_changed(self.update)
        self.UM_slider.on_changed(self.update)
        self.MH_slider.on_changed(self.update_mh_age)
        self.AGE_slider.on_changed(self.update_mh_age)
        
        # Set labels of slider AGE
        label_age.axis('off')
        self.label_age = label_age.text(0.5, 0.5, '', ha='center', va='center',fontsize=10)
        self.label_age.set_text(f"Age range\n{self.string_age(self.AGE[0])}\n{self.string_age(self.AGE[-1])}")

        # Set labels of slider MH
        label_mh.axis('off')
        self.label_mh = label_mh.text(0.5, 0.5, '', ha='center', va='center',fontsize=10)
        self.label_mh.set_text(f"MH range\n{self.MH[0]} dex\n{self.MH[-1]} dex")

        # Special Buttons
        self.btn_new_iso=Button(ax_new_iso,'Open other Isochrone')
        self.btn_kiel_dg=Button(ax_kiel,   'Obtain\n stellar parameters')

        self.btn_new_iso.on_clicked(lambda event: self.run_app())
        self.btn_kiel_dg.on_clicked(lambda event: self.match_nearest())

    #___________________________________________________________________
    # Function that creates a readable Age (asumes logAge=log(Age[yr]) )
    def string_age(self,logAge:float) -> str: 
        if logAge>=9:
            return f'{10**(logAge-9):6.2f} Gyr'
        else:
            return f'{10**(logAge-6):6.2f} Myr'
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

    # Function that updates the Isochrones used when changed the MH or AGE
    def update_mh_age(self,event=None):
        mh_min, mh_max = map(int, self.MH_slider.val)
        mh_min, mh_max = self.MH[mh_min],self.MH[mh_max]
        self.label_mh.set_text(f"MH range\n\t{mh_min:^6.2f} dex\n{mh_max:^6.2f} dex")

        age_min, age_max = map(int, self.AGE_slider.val)
        age_min, age_max = self.AGE[age_min],self.AGE[age_max]
        self.label_age.set_text(f"Age range\n{self.string_age(age_min)}\n{self.string_age(age_max)}")
        
        self.Iso_mask = (age_min<=self.all_ages)&(age_max>=self.all_ages)&(mh_min<=self.all_mh)&(mh_max>=self.all_mh)

        self.iso_scatter.set_array(self.Iso_mask)
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()   

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
                       #
                       ) -> None:
        """
        Function thougth to work with PARSEC ISOCHRONES

        It will use the stages defined on __init__  

        it has the option update for the firt iteration 
        
        """
        print("\n\n\nIsochrones with multiples models on the same file are large files, please wait till the entire files is loaded into the ram (readed)\n\n\n")

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

        self.all_ages=self.Iso['logAge']
        self.all_mh  =self.Iso['MH']

        # set up variables used in the plots and interpolations
        self.Mag_iso  =  self.Iso[self.system['Mag1']]
        self.Color_iso= (self.Iso[self.system['Mag2']]-self.Iso[self.system['Mag3']])
        
        self.Iso_mask = np.repeat(True,len(self.Mag_iso))

        # Init values for display purposes
        AV= (np.mean(self.COLOR)-np.mean(self.Color_iso))/self.abs_col
        if AV<0: # you can not have negative absorption
            AV=0

        self.AV= 0.01#AV 
        self.UM= 13.37#(np.mean(self.MAG  )-np.mean(self.Mag_iso)  )

        # Shift the isochrone near the data
        self.fitted_isochrone()

    #______________________________________________________________________
    # Funtion to ask for path if you desire to open another isochrone
    def ask_path(self) -> str:
        path: str =fd.askopenfilename(title='Select the file of PARSEC isochrone') 
        return path
        

    #____________________________________________________________________________________
    # Function that computes representative values from isochrone points to the data
    # It will store the data in a file and after it will plot the interpolated results that will pop up in another window
    # Min distance to the to point is defined (if the points to be fitted are not in middle of the models has no sense to infere stellar parameters)
    # Max distance to the point to be fitted is defined 
    # (for sources that are farther than these values form the nearest point of the models can not be retrived reliable estimations and computationally expensive)
    def match_nearest(self,Output_file:str='Isochrone_values',
                      min_mag_dist=0.075,min_color_dist=0.035,
                      max_mag_dist=1 ,max_color_dist=0.3,
                      ) -> None:
            
            self.MH_range =tuple(map(lambda x:   self.MH[int(x)],   self.MH_slider.val))
            self.AGE_range=tuple(map(lambda x:  self.AGE[int(x)],  self.AGE_slider.val))

            plt.close()
            # self.waiting_figure()

            # Target data
            MAG   = self.MAG
            COLOR = self.COLOR

            # Data of isochrone to be used form the CMD
            iso_mag=self.Fit_Mag_iso  [self.Iso_mask]
            iso_col=self.Fit_Color_iso[self.Iso_mask]

            inv_cov=self.inverse_covariance_matrix_cmd(
                        err_mag1=self.err_mag1,
                        err_mag2=self.err_mag2,
                        err_mag3=self.err_mag3,
                        err_dist=self.err_dist,
                        err_redd=self.err_redd,
            )

            check_inv=(inv_cov.ndim==2)
            
            # while adding the data we will use a dictionary
            inferred_params=dict()
            inferred_params['good']  = np.empty(len(MAG), dtype=np.float64)
            inferred_params['d_min'] = np.empty(len(MAG), dtype=np.float64)
            inferred_params['NP_mod']= np.empty(len(MAG), dtype=np.float64)
            for key in self.Iso.dtype.names:
                inferred_params[key         ]=np.empty(len(MAG), dtype=np.float64)
                inferred_params[key+'_err_l']=np.empty(len(MAG), dtype=np.float64)
                inferred_params[key+'_err_u']=np.empty(len(MAG), dtype=np.float64)

            #------------------
            N_data=len(MAG)
            N_iso=len(iso_mag)

            cpu_count = mp.cpu_count()
            chunk_size = (N_data + cpu_count - 1) // cpu_count
            chunks = [list(range(i, min(i + chunk_size, N_data))) for i in range(0, N_data, chunk_size)]
            print(cpu_count)
    
            with mp.Pool(processes=cpu_count,
                initializer=init_worker,
                initargs=(MAG, COLOR, N_iso, iso_mag, iso_col,
                          max_mag_dist, max_color_dist, min_mag_dist, min_color_dist,
                         list(inferred_params.keys()), self.Iso,
                         self.Iso_mask, check_inv, inv_cov, self.degree,self.N_bins) 
                         ) as pool:

                # chunk_results = pool.map(process_chunk, chunks)
            
            # for chunk in chunk_results:
            #     for i, res in chunk.items():
            #         for key in res:
            #             inferred_params[key][i] = res[key]
                for i, res in pool.imap_unordered(process_single_point, range(N_data)):
                    for key in res:
                        inferred_params[key][i] = res[key]
                    sys.stdout.write(f'\r{100*i/N_data:.3f}%')
                    sys.stdout.flush()

            
            inferred_params['Skepticism_score']= ((inferred_params['d_min']/(min_mag_dist**2+min_color_dist**2))**0.5)*((1-(inferred_params['NP_mod']/self.N_bins))**0.5)
            inferred_params['scale_err']=1+(inferred_params['Skepticism_score']**0.5)

            # for i in range(N_data):
            #     # We create an array that stores the weigths of the isochrones for each data point
            #     weigths   = np.zeros(N_iso)
            #     # Data point
            #     mag,color = MAG[i],COLOR[i]
                
            #     # We trim the isochrone in the nerby zone (easier to compute)
            #     mask  = ((iso_mag<mag+max_mag_dist)&
            #              (iso_mag>mag-max_mag_dist)&
            #              (iso_col<color+max_color_dist)&
            #              (iso_col>color-max_color_dist))

            #     # Distance in data and isochrone
            #     d_mag = (iso_mag[mask]-mag)
            #     d_col = (iso_col[mask]-color)

            #     # If data not embedded inside the parameter space it can be inferred stellar parameters
            #     if sum(mask)==0 or (abs(d_mag/min_mag_dist)+abs(d_col/min_color_dist)).min()>2:
            #         for key in inferred_params:
            #             inferred_params[key][i]=np.nan
            #         inferred_params['good'][i]=0 # We will inform in our output file that we do not compute for this source
            #         continue # We skip this point
                
            #     try:# If anything goes wrong we will skip the point
            #         # We will work with just the neares point around the datapoint in around 360 degrees
            #         idx=self.line_of_sight_of_point(delta_y=d_mag,delta_x=d_col)

            #         # We obtain the weigths for each model point
            #         if check_inv:
            #             weigths[np.where(mask)[0][idx]]= self.log_likelihood_cmd(
            #                                                         delta_mag      = d_mag[idx],
            #                                                         delta_col      = d_col[idx],
            #                                                         inv_cov_matrix = inv_cov,
            #                                                         N              = self.degree)
            #         else: # If we have N dimensional covariance matrix
            #             weigths[np.where(mask)[0][idx]]= self.log_likelihood_cmd(
            #                                                         delta_mag      = d_mag[idx],
            #                                                         delta_col      = d_col[idx],
            #                                                         inv_cov_matrix = inv_cov[i],
            #                                                         N              = self.degree)
            #         # We get representative values an errors
            #         for key in self.Iso.dtype.names:
            #             mean,std_l,std_u=self.weighted_mean_and_std(self.Iso[key][self.Iso_mask],weigths)
            #             # We store the data
            #             inferred_params[key][i]          = mean
            #             inferred_params[key+'_err_l'][i] = std_l
            #             inferred_params[key+'_err_u'][i] = std_u
                    
            #         # We store that we calculate this column
            #         inferred_params['good'][i]=1
                
            #     except:
            #         for key in inferred_params:
            #             inferred_params[key][i]=np.nan
            #         inferred_params['good'][i]=-1 # We will inform in our output file that we do not compute for this source
            #         continue # We skip this point

            #     # print('\n#______\n',i)
            #     # print(sum(mask),len(idx))
            #     # # if (i==5) or (i==11):
            #     # if (i==11):
            #     #     print(i,i/N_data,flush=True)
            #     #     # plt.figure()
            #     #     # plt.scatter(iso_col[mask],weigths[mask],color='cyan')
            #     #     # plt.figure()
            #     #     # plt.scatter(iso_mag[mask],weigths[mask],color='cyan')
            #     #     plt.figure()
            #     #     plt.scatter(iso_col[mask],iso_mag[mask],cmap='nipy_spectral_r',c=weigths[mask])
            #     #     # plt.scatter(iso_col[mask][idx],iso_mag[mask][idx],color='green')
        
            #     #     plt.scatter(color,mag,color='r',marker='X')
            #     #     plt.gca().invert_yaxis()
            #     #     plt.figure()
            #     #     plt.scatter(self.Iso['logg'][self.Iso_mask],weigths,color='blue')
            #     #     mean,std_l,std_u=self.weighted_mean_and_std(self.Iso['logg'][self.Iso_mask],weigths)
            #     #     plt.axvline(mean,color='cyan')
            #     #     plt.axvline(mean+std_u,color='fuchsia')
            #     #     plt.axvline(mean+std_l,color='fuchsia')
            
            #     #     plt.figure()
            #     #     plt.scatter(10**self.Iso['logTe'][self.Iso_mask],weigths,color='red')
            #     #     mean,std_l,std_u=self.weighted_mean_and_std(10**self.Iso['logTe'][self.Iso_mask],weigths)
            #     #     plt.axvline(mean,color='cyan')
            #     #     plt.axvline(mean+std_u,color='fuchsia')
            #     #     plt.axvline(mean+std_l,color='fuchsia')
            #     #     plt.pause(300)
                
            #     # if i>=11:
            #     #     import sys
            #     #     plt.close()
            #     #     sys.exit()

            #     # Progress Bar
            #     symbol = spinner[i % 4]
            #     # print(symbol,i/N_data,flush=True)
            #     sys.stdout.write(f'\r{symbol}  {100*i/N_data:.3f}%')
            #     sys.stdout.flush()

            #---------For each data point 
            #__________________________________________
            
            self.save_fits(params=inferred_params,Fname=Output_file)
            # self.save_plain_txt(params=inferred_params,Fname=Output_file)
            #-------------------
            # we end our calculations with displaying a kiel diagram with our inferred params
            # self.kiel_diagram(10**inferred_params['logTe'],inferred_params['logg'])




    #----------------------------------
    # Function that stores the inferred params in a fits file
    def save_fits(self,
                  params:dict,
                  Fname:str,
                  path ='./'):
        
        header = fits.Header()

        # Info about the isochrones
        header['Iso']    = (self.iso_path    , 'Isochrone used')
        header['MH_l']   = (self.MH_range[0] , 'MH min of range used')
        header['MH_u']   = (self.MH_range[1] , 'MH max of range used')
        header['logA_l'] = (self.AGE_range[0], 'logAge min of range used')
        header['logA_u'] = (self.AGE_range[1], 'logAge min of range used')
        header['AV']     = (self.AV          , 'Visual Absorption used' )
        header['UM']     = (self.UM          , 'Distance Modulus used' )
        header['AbsMag'] = (self.abs_mag     , 'Absortion coeff Ay/Av in mag')
        header['AbsCol'] = (self.abs_col     , 'Absortion coeff Ay/Av in color')
        
        # Comments Example
        header['COMMENT'] = "File done by I. Baeza"

        # Create primary hdu
        primary_hdu = fits.PrimaryHDU(header=header)
            
        columns=[]
        if not self.OBJ is None:
            columns.append(fits.Column(name='OBJ',format='A99', unit='id'    ,  array= self.OBJ ))
        if not self.RA is None:
            columns.append(fits.Column(name='RA', format='D'  , unit='degree',  array= self.RA ))
        if not self.DEC is None:
            columns.append(fits.Column(name='DEC',format='D'  , unit='degree',  array= self.DEC ))
         
        columns.append(fits.Column(name='MAG'  ,format='D'  , unit='mag',  array= self.MAG ))
        columns.append(fits.Column(name='COLOR',format='D'  , unit='mag',  array= self.COLOR ))

        for key in params:
            columns.append(fits.Column(name=key, format='D', array= params[key] ))
            
        # Generate the secondary fits with the data
        NEWFITS = fits.BinTableHDU.from_columns(columns)

        # Set the EXTNAME keyword to label the binary table HDU
        NEWFITS.header['EXTNAME'] = 'Data'

        # Save the fits file
        hdul = fits.HDUList([primary_hdu, NEWFITS])
        hdul.writeto(path+Fname+'.fits', overwrite=True)

    #----------------------------------
    # Function that stores the inferred params in a plain txt
    def save_plain_txt(self,
                       params:dict,
                       Fname:str,
                       path ='./'):
        
        table = PrettyTable()

        if not self.OBJ is None:
            table.add_column('Obj',self.OBJ)
        if not self.RA is None:
            table.add_column('RA',self.RA)
        if not self.DEC is None:
            table.add_column('DEC',self.DEC)
        
        table.add_column('MAG'  ,self.MAG)
        table.add_column('COLOR',self.COLOR)

        for key in params:
            table.add_column(key,params[key])

        # Write the file with the info
        with open(path+Fname+'.txt','w') as file:
            file.write(f'# DATA USED TO INFER THESE PARAMETERS\n#\n')
            file.write(f'# Isochrone Used ={self.iso_path} \n')
            file.write(f'# MH             ={self.MH_range} \n')
            file.write(f'# logAGE         ={self.AGE_range}\n')
            file.write(f'# AV             ={self.AV}\n')
            file.write(f'# UM             ={self.UM}\n')
            file.write(f'# abs_coeff_Mag  ={self.abs_mag}\n')
            file.write(f'# abs_coeff_Color={self.abs_col}\n')
            file.write(f'#\n#----------------------------------------------------\n')
        
            # Write the data in a file
            file.write('#'+table.get_string(border=False)[1:])
            
    #_____________________________________________________________________
    # warning plot to explicitly say that the computer is doing somenthing
    def waiting_figure(self,id=1704) -> None:
        plt.close()
        FIG= plt.figure(id,figsize=(12,4))
        ax =FIG.add_subplot(xticks=[], yticks=[])
        # The first word, created with text().
        text = ax.text(-.1, .5, 'The stellar parameters are being estimated,\n please wait till another window pop up,\n A warning might appear saying "matplotlib is not responding" but just wait', color="k",fontsize=20)
        ax.spines[:].set_visible(False)
        plt.show(block=True)
        plt.pause(0.1)

        print("\n\nRunning the stellar parameters ...\n\n")

    # Figure to create the Kiel diagramm 
    def kiel_diagram(self,teff_inferred,logg_inferred,id=1704):
        plt.close()
        plt.figure(id)
        plt.scatter(teff_inferred,logg_inferred)
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.xlabel('Teff')
        plt.ylabel('Log g')
        plt.title('Kiel Diagram')
        plt.show(block=False)


if __name__=='__main__':
    from astropy.io import fits

    path_data='/home/ian/Desktop/Data_base/'
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

    
    iso3='PARSEC_V2_10-13.5Gyr__-1_0.51MH_n0.48_stetson.txt'
    iso3='Parsec_V2_10_13Gyr__-0.78_-0.41MH__n_0.35.dat'
    iso3='Parsec_V2_10_13Gyr__-0.78_-0.41MH__n_0.4.txt'

    ISO_FIT(MAG=I,COLOR=V-I,OBJ=id,RA=ra,DEC=dec,iso_path='/home/ian/Downloads/'+iso3,photometry="ubvir",Header=14)

# MH =-0.59
# AGE=10.16879
# AV =0.01
# UM =13.219999999999997
# abs coeff Mag  =0.59893
# abs coeff Color=0.4020300000000001

# exp10(logTe+logTe_err_U*pow(scale_err,1))
# 
# 
# logg+logg_err_L*pow(scale_err,1)
# 
# MAG- I_fSBmag+I_fSBmag_err_L*pow(scale_err,1)*3
#  