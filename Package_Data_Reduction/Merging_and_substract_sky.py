import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Angle
from astropy import units as u

from scipy.interpolate import interp1d

import os
import sys



###########################################################################3

# This is though to merge spectra of stars taken with the same telescope and with no significant variations in between epochs of their radial velocities. Barycentric correction is applied to the wavelengths, also this assumes that similar data is given due we will take the corrected wavelenghts of the first epoch/spectrum as reference
def MERGING_MULTIPLE_SPECTRA(MJD,
                             path,
                             Wave,
                             Flux,
                             OFiles,
                             RAS,
                             DECS,
                             GLAT,
                             GLON,
                             HTEL,
                             function_used_to_merge=np.median,
                             NOT_SKY   = None,
                             OBJS      = None,
                             MAGS      = None,
                             SAVE      = True,
                             SHOW      = True,
                             SHOW_SKIES= True,
                             window    = [None],
                             AXIS_TO_MERGE=0,
                             wave_to_angstrom=1,
                             coord_sys='fk5', # usually used by ESO
                             ):
    #-----------------------------------------------------------
    # Check options of the function
    if NOT_SKY is None:
        NOT_SKY=np.repeat(True,len(MJD))
    
    no_obj=False
    if OBJS is None:
        OBJS   = np.arange(len(MJD))
        no_obj = True
    
    no_mags = False
    if MAGS is None:
        MAGS   = [np.median(i) for i in Flux[0]]
        no_mags= True

    #-----------------------------------------------------------
    # We will Check that the objects match in shape and Observing block

    try:
        Flux[0][NOT_SKY[0]]
    except:
        for i in range(len(Flux)):
            Flux[i]=Flux[i].T

    SORT_ACCORDINGLY_TO_FIRST_EPOCH=[np.arange(len(Flux[0][NOT_SKY[0]]))]
    # to check that every file matches the whole sample of science targets
    AUX=[] # we will use this auxiliar variable

    first_measure= SkyCoord(ra=RAS[0]  [NOT_SKY[0]]  *u.degree, dec=DECS[0]  [NOT_SKY[0]]  *u.degree)
    
    for i in range(len(MJD)-1):
        next_measure = SkyCoord(ra=RAS[i+1][NOT_SKY[i+1]]*u.degree, dec=DECS[i+1][NOT_SKY[i+1]]*u.degree)

        indexes,d2d,d3d=first_measure.match_to_catalog_sky(next_measure)
        # Threshold of Separation
        sep= d2d>0.5*u.arcsec
        # If the all the targets (not sky) matches the position and the observations has the same elements we confirm that they belong to the same Obs Block
        if (not sum(sep)) and (len(first_measure)==len(next_measure)):
            SORT_ACCORDINGLY_TO_FIRST_EPOCH.append(indexes)
            AUX.append(False)
        else:
            AUX.append(True)
            
    # If any of the science targets does not match or we have at least one traget different we will stop
    if sum(AUX)!=0: # In case that does not match the RADEC independently of fiber position or number of sky fibers
        AUX=np.array(AUX)
        print('These files: ',MJD[~AUX])
        print('Does not match these files:',MJD[AUX])
        raise Exception('The Values of RA DEC does not match between files (different number of science targets or not the same targets)')

    # We will store the science target and the sky separately
    DATAS=[]
    SKIES=[]

    AUX=[] #This will help us to keep trace of the nearest sky fiber in ccd (assuming that were given in this order, applies for girrafe)
    # Separate The data from the sky
    for i in range(len(MJD)):
        or_idx=np.arange(len(Flux[i]))
        data,sky=Flux[i][NOT_SKY[i]],Flux[i][~NOT_SKY[i]]

        sci_idx = or_idx[NOT_SKY[i]]          # (Nscience,)
        sky_idx = or_idx[~NOT_SKY[i]]         # (Nsky,)

        # dist = np.abs(sci_idx[:, None] - sky_idx[None, :]).astype(float)   # (Nscience, Nsky)
        pos_stars=SkyCoord(ra=RAS[i][NOT_SKY[i]]*u.degree, dec=DECS[i][NOT_SKY[i]]*u.degree)
        dist=[]
        
        for ra,dec in zip(RAS[i][~NOT_SKY[i]],DECS[i][~NOT_SKY[i]]):
            d2d=pos_stars.separation(SkyCoord(ra=ra  *u.degree, dec=dec *u.degree))
            dist.append(d2d)
            # print(d2d)
        
        dist=np.array(dist).T

        print(dist.shape)

        dist = (np.abs(np.max(dist, axis=1)[:, None] - dist))**40  # invert distance
        dist /= np.sum(dist, axis=1)[:, None]                      # normalize
                
        # print(dist)
        AUX.append(dist.copy())

        DATAS.append(data[SORT_ACCORDINGLY_TO_FIRST_EPOCH[i]])
        print(data.shape,end='\t')
        if not NOT_SKY is None and sum(~NOT_SKY[i])>0:
            print(sky.shape)
            SKIES.append(sky)
        else:
            print(f'\n\n# WARNING: This file:\"{OFiles[i]}\" was encountered to not have sky...\n\t I will proceed anyway but you were warned about\n')
            SKIES.append(np.zeros_like(data))

    try:
        array   = np.array(DATAS) # Must be the same shape
        DATAS   = array.copy()
    except:
        DATAS   = MISMATCH_SECOND_DIMENSION(DATAS,Wave) # Thought just in the case of UVES fibre
        SKIES   = MISMATCH_SECOND_DIMENSION(SKIES,Wave)
    # Not necessary to have the same number of fibers
    array   = np.empty(len(SKIES), dtype=object)
    array[:]= SKIES
    SKIES   = array 
    #--------------------------------------------------------------
    #--------------------------------------------------------------
    # Here we merge the files
    Data = []
    Sky = []
    counts = []
    RVCORR = []
    
    sc = SkyCoord(ra=RAS[0][NOT_SKY[0]]*u.degree, dec=DECS[0][NOT_SKY[0]]*u.degree, frame=coord_sys)
    TEL = EarthLocation(lon=GLON[0], lat=GLAT[0], height=HTEL[0])
    
    WMIN, WMAX = min(list(map(np.nanmin, Wave))), min(list(map(np.nanmax, Wave)))
    
    # SHAPE FIX 1: Use [-1] to get the number of pixels, not the number of epochs
    N_pixels = len(Wave[0])

    
    samp_wave = extend_linear_range(WMIN*0.9999999, WMAX*1.000000001, 3 * N_pixels)
    
    # SHAPE FIX 2: Rename this to avoid overwriting it inside the loop later
    # master_temp_wave = extend_linear_range(WMIN, WMAX, 3*N_pixels) 
    master_temp_wave = samp_wave


    for i in range(len(MJD)):
        # 1. Get Barycentric Correction 
        # Force it to be at least a 1D array so we can safely check its length/shape
        rv_val = sc.radial_velocity_correction(kind='barycentric', obstime=Time(MJD[i], format='mjd'), location=TEL).to(u.km/u.s).value
        rvcorr = np.atleast_1d(rv_val) 
        RVCORR.append(rvcorr)

        # SHAPE FIX 3: Reshape rvcorr to (N_spectra, 1) so it broadcasts properly across the N_pixels axis
        # doppler_factor shape is now: (N_spectra, 1)
        doppler_factor = (1 + rvcorr / 299792.458).reshape(-1, 1) 

        try: 
            # If Wave[i] is (N_spectra, N_pixels) or (N_pixels,)
            # Multiplying by (N_spectra, 1) guarantees an output of (N_spectra, N_pixels)
            temp_waves = Wave[i] * doppler_factor
        except: 
            # If Wave[i] is just a flat 1D array of N_pixels, tile it to match N_spectra
            # Use temp_waves (plural) so the zip loop below works!
            temp_waves = np.tile(Wave[i], (len(rvcorr), 1)) * doppler_factor

        # Get sky and subtract it. current_data shape: (N_spectra, N_pixels)
        sky = function_used_to_merge(SKIES[i], axis=AXIS_TO_MERGE).astype(float)
        current_data = DATAS[i] - sky
        
        ########################
        # EXPERIMENTAL (scale every spectra to have a median equal to one)
        # var shape: (N_spectra, 1)
        var = np.nanmedian(current_data, axis=1).reshape(-1, 1)
        counts.append(var)
        
        # Avoid division by zero
        var[var == 0] = 1
        
        # Broadcasting (N_spectra, N_pixels) / (N_spectra, 1) works perfectly here
        current_data = current_data / var

        # Barycentric correction interpolation
        CURRENT_DATA = []
        # Zip safely unpacks the first dimension: wave and spectrum are both 1D (N_pixels,)
        for wave, spectrum in zip(temp_waves, current_data):
            # Interpolate onto the expanded samp_wave grid
            interp_func = interp1d(wave, spectrum, bounds_error=False, fill_value=np.nan)
            CURRENT_DATA.append(interp_func(samp_wave))
            
        ########################
        Data.append(np.array(CURRENT_DATA))

        # We already calculated sky above, but if you need to append it again:
        Sky.append(sky)

    # Flatten the counts
    COUNTS = function_used_to_merge(counts, axis=AXIS_TO_MERGE).ravel()

    # Merge the final Data arrays
    Data = function_used_to_merge(Data, axis=AXIS_TO_MERGE)
    
    Data2 = []
    for i in Data:
        # i is 1D array (len(samp_wave),)
        # We interpolate back onto master_temp_wave! (Using the protected variable name)
        y = interp1d(samp_wave, i, bounds_error=True)(master_temp_wave)
        y[np.isnan(y)] = 0
        Data2.append(y)
    
    Data = np.array(Data2)
    
    #--------------------------------------------------------------

    RVCORR=np.array(RVCORR)
    ###########################################################################3

    if SAVE:
        for i in range(len(SKIES)):
            RA =RAS [i][~NOT_SKY[i]]
            DEC=DECS[i][~NOT_SKY[i]]
            
            Fname=f'SKY_{MJD[i]}.fits'
      
            save_sky(skies_obs=SKIES[i],
                    Fname  = Fname,
                    OBJ    = "SKY",
                    RA     = RA,
                    DEC    = DEC,
                    Ofile  = OFiles[i],
                    MJD    = MJD[i],
                    flux   = SKIES[i],
                    wave   = Wave[i],
                    GLAT   = GLAT[i],
                    GLON   = GLON[i],
                    HTEL   = HTEL[i],
                    path   = path+"skies/",
                    wave_to_angstrom=wave_to_angstrom)

        # Create an oficial list of files
        LIST_FILES=open(path+'/Oficial_list.list','w')
        INFO_FILES=open(path+'/INFO_FILES.txt','w')
        INFO_FILES.write('#File_name\tRA\tDEC\tMAGS\tN_counts\n')
        mjd =np.mean(MJD)

        date=Time(mjd,format='mjd')
        date=date.isot
        
        # We just keep the data of the first epoch considered in the observation block to fill the header 
        GLAT=GLAT[0]   
        GLON=GLON[0]   
        HTEL=HTEL[0]  

        RAS =RAS [0][NOT_SKY[0]]
        DECS=DECS[0][NOT_SKY[0]]
        copy_wave=Wave[0].copy()
        Wave =[master_temp_wave]
        if not no_obj:
            OBJS=OBJS[0][NOT_SKY[0]]
        if not no_mags:
            MAGS=MAGS[0][NOT_SKY[0]]
        
        for i in range(len(Data)):
            OBJ=OBJS[i]
            RA =RAS [i]
            DEC=DECS[i]
            
            MAG=MAGS[i]
            COUNT=COUNTS[i]

            Fname=str(OBJ)+'='+str(RA)+'='+str(DEC)+'='+date+'.fits'

            new_fits(Fname  = Fname,
                     OBJ    = OBJ,
                     RA     = RA,
                     DEC    = DEC,
                     MJD    = MJD,
                     MAG    = MAG, 
                     Ofile  = OFiles,
                     MJDFile= MJD,
                     Flux   = Data[i],
                     Wave   = Wave[0],
                     GLAT   = GLAT,
                     GLON   = GLON,
                     HTEL   = HTEL,
                     path   = path,
                     count  =COUNT,
                     wave_to_angstrom=wave_to_angstrom,
                     RVCORR = RVCORR[:,i])
    
            LIST_FILES.write(Fname+'\n')
            INFO_FILES.write(f'{Fname}\t{RA}\t{DEC}\t{MAG}\t{COUNT}\n')
        
        LIST_FILES.close() # Close the oficial list of files
        INFO_FILES.close()

    ###########################################################################3

    if SHOW or SAVE:
        Data=Data*COUNTS.reshape(-1,1)
        Fig= plt.figure(1999,figsize=(13,6))
        gs = gridspec.GridSpec(1,2,             # Number of axis y,x
                        height_ratios = [1],    # relatives ratios of heigh
                        width_ratios  = [1,0.05],  # relatives ratio of with
                        left  = 0.1,      # Space to the edge of left   from the nearest axis
                        right = 0.93,     # Space to the edge of right  from the nearest axis
                        bottom= 0.08,     # Space to the edge of bottom from the nearest axis
                        top   = 0.95,     # Space to the edge of top    from the nearest axis
                        wspace= 0.2,      # Space horizontal between each of the axis
                        hspace= 0.2)      # Space vertical   between each of the axis

        ax1=Fig.add_subplot(gs[0,0])
        cax=Fig.add_subplot(gs[0,1])
        cax.invert_yaxis()
        # We will color the spectras accordigly to the their magnitudes
        try:
            MAGS  = np.array(list(map(float,MAGS[0][NOT_SKY[0]])))
        except:
            MAGS  = np.array(list(map(float,MAGS)))

        colors= (MAGS-MAGS.min())/(MAGS.max()-MAGS.min())
        cmap  = plt.cm.gist_rainbow(colors)
        
        for i in range(len(Data)):
            ax1.plot(Wave[0],Data[i],color=cmap[i])
        
        try:
            ax1.set_xlim(window[0],window[-1])
            ax1.set_ylim(np.nanmin(COUNTS)*0.95,np.nanmax(COUNTS)*1.05)
        except ValueError as VL:
            print('\n\n\n')
            print(VL)
            print('\n\n\n')
            print(Fname,Data.shape,sky.shape)
            print(sum(np.isnan(Data)),sum(np.isnan(sky)))
            print(DATAS,SKIES)
            plt.show()
            sys.exit()

        norm= mpl.colors.Normalize(vmin=min(MAGS), vmax=max(MAGS), clip=False)
        SM  = mpl.cm.ScalarMappable(norm=norm, cmap='gist_rainbow')
        SM.set_array([])
        
        LABEL= 'median Counts' if no_mags else 'Mags' 
        
        plt.colorbar(mappable=SM,label=LABEL,extend='both',cax=cax)

        if SAVE:
            print('Saving Figure')
            plt.savefig(path+'/Substract_SKY.png',dpi=150)
            plt.close(1999)
        if SHOW:
            plt.show()
        else:
            plt.close(1999)

    ###########################################################################3

    if SHOW_SKIES or SAVE and not NOT_SKY is None:
        Wave=[copy_wave]
        #-----------------
        Fig= plt.figure(1999,figsize=(17,9))
        gs = gridspec.GridSpec(len(MJD),2,                 # Number of axis y,x
                        height_ratios = [1]*len(MJD), # relatives ratios of heigh
                        width_ratios  = [1,1],          # relatives ratio of with
                        left  = 0.1,      # Space to the edge of left   from the nearest axis
                        right = 0.95,     # Space to the edge of right  from the nearest axis
                        bottom= 0.08,     # Space to the edge of bottom from the nearest axis
                        top   = 0.95,     # Space to the edge of top    from the nearest axis
                        wspace= 0.2,      # Space horizontal between each of the axis
                        hspace= 0.2)      # Space vertical   between each of the axis
        
        for q in range(len(MJD)):
            ax=Fig.add_subplot(gs[q,0])
            for i in range(sum(~NOT_SKY[q])):
                ax.plot(Wave[0],SKIES[q][i],label=str(q))
            ax=Fig.add_subplot(gs[q,1],sharey=ax,sharex=ax)
            ax.plot(Wave[0],Sky[q],color='black',label=OFiles[q])
            ax.legend()
            ax.set_ylim(Sky[q].min()*0.9,Sky[q].max()*1.5)
            # ax.legend()
        #-----------------
        #---------------------
        if SAVE:
            print('Saving Figure')
            plt.savefig(path+'SKY.png',dpi=150)
            plt.close(1999)
        if SHOW:
            plt.show()
        else:
            plt.close(1999)
    



###########################################################################3

# Function used in case that data had a diferrent shape in one dimension (e.g. different night of calibration on UVES)
def MISMATCH_SECOND_DIMENSION(fluxes:list,waves:list,tolerance:float=0.02) -> list:
    print('WARNING:\n\t The files that you have inputed to merge does not match in their dimension, I will assume that you know what you are doing and I will interpolate the data to match the first file (do not trust the edges)')
    print('\n\tThis has been only though at the begining just for UVES FIBRE that under different calibrations ends up with slightly different number of pixels but either way for each exposure all files will share the same CRVAL and CDELT\n')
    if len(fluxes)!=len(waves):
        raise Exception('You have not inputed the required wavelengths to perform this workaround')
    
    min_w=np.array([np.nanmin(i) for i in waves])
    max_w=np.array([np.nanmax(i) for i in waves])
    diff_min=max(min_w)-min(min_w)
    diff_max=max(max_w)-min(max_w)
    m_width=np.nanmedian(max_w-min_w)
    if tolerance<diff_min/m_width or tolerance<diff_max/m_width:
        raise ValueError(f"The ends and/or the begining of the spectra that you wanted to force to merge differ more than the tolerance accepted of {tolerance*100}%\n\n Please check this and if you want to continue relax the tolerance\n\nThe new required value for this to work must be for the min waves {100*diff_min/m_width}%, and for the max waves {100*diff_max/m_width}% ")

    new_fluxes=[fluxes[0]]
    for i in range(1,len(fluxes)):
        edge=np.median(fluxes[i])
        temp_flux= interp1d(waves[i],fluxes[i],bounds_error=False,fill_value=edge)(waves[0])
        
        new_fluxes.append(temp_flux)

    return new_fluxes



###########################################################################3

# Function used to create a new fits file containing the spectra after being sky substrcated and merged (if applies)
def new_fits(Fname,
             OBJ,
             RA,DEC,
             MJD,
             Ofile,
             MJDFile,
             Flux,Wave,
             GLAT,GLON,HTEL,

             CRVAL= None,
             CDELT= None,
             MAG  = None,
             path ='./',
             count=1,
             wave_to_angstrom=1,
             RVCORR=[0]):
    # Create a FITS header
    header = fits.Header()

    # Info about the star
    header['OBJECT']  = (OBJ,   'Object name in parent file 1')
    header['RA']      = (RA,    'Right Ascension /'+str(Angle(float(RA)*u.degree).to('hourangle')))
    header['DEC']     = (DEC,   'Declination     /'+str(Angle(float(DEC)*u.degree).to('hourangle')))
    date=Time(np.mean(MJD),format='mjd')
    date=date.isot
    header['MJD']     = (np.mean(MJD),   'MJD (mean of files)| ' + str(date))
            
    if not MAG is None:
        header['MAG']     = (MAG,   'MAG Associated in original file')
    
    #Info about the parent files
    for i,j in enumerate(Ofile):
        header['OFIL'+str(i)]   = (j,f'Original parent File N:{i+1} of {len(Ofile)}')
    for i,j in enumerate(MJDFile):
        header['DFIL'+str(i)]  = (j, f'Original parent File MJD N:{i+1} of {len(Ofile)}')

    # Info about telescope
    header['GEOLON']  = (GLAT,  'Geo Latitude of Telescope')
    header['GEOLAT']  = (GLON,  'Geo Longitude of Telescope')
    header['GEOLEV']  = (HTEL,  'Height over sea of Telescope')
    
    # original number of proportionally counts that you have to multiply your file to get the real (not scaled to 1 spectra)
    header['MCOUNT'] = (count,'Multiply flux by this value to get counts rather than relative flux')

    header['RVBARY'] = (np.mean(RVCORR),f'RV barycentric Correction (mean)')
 
    for j,rvcorr in enumerate(RVCORR): 
        header[f'BARY{j}']  = (RVCORR[j],f'RV barycentric Correction {j}')

        # We keep our track of the orignial times
        date=Time(MJD[j],format='mjd')
        date=date.isot
        header[f'MJD{j}']     = (MJD[j],   f'MJD epoch {j} | ' + str(date))

    # Comments Example
    header['COMMENT'] = "File done by I. Baeza"

    # Create primary hdu
    primary_hdu = fits.PrimaryHDU(header=header)
    
    # If we provide the values we also write fits file in iraf like 
    if not (CRVAL is None) and not (CDELT is None):
        header['CRVAL1']  = (CRVAL, 'wavelenght first pixel ')
        header['CDELT1']  = (CDELT, 'Step of wavelenght ')

        primary_hdu.data=Flux

    # Create the columns of data (corrected by barycentric velocity)
    columns=[
        fits.Column(name='wave', format='D', unit='Angstrom',      array= (Wave*wave_to_angstrom)),
        fits.Column(name='flux', format='D', unit='relative flux', array= Flux )
        ]

    # Generate the secondary fits with the data
    NEWFITS = fits.BinTableHDU.from_columns(columns)

    # Set the EXTNAME keyword to label the binary table HDU
    NEWFITS.header['EXTNAME'] = 'SPECTRUM (sky sub)'
    # Comments about the columns (topcat)
    NEWFITS.header['TCOMM1']  = 'wavelenght in angstrom corrected by barycentric vel'
    NEWFITS.header['TCOMM2']  = 'relative Flux wiht sky substracted' 
    NEWFITS.header['COMMENT'] = f'Object: {OBJ}' 

    # Save the fits file
    hdul = fits.HDUList([primary_hdu, NEWFITS])
    hdul.writeto(path+Fname, overwrite=True)


###########################################################################3

# Function used to get a good distribution of the plots
def malla(x,PRINT=0):
    aux   = (x**0.5)+0.5
    alpha = int(aux)
    beta  = alpha + round(aux-int(aux))
    if PRINT:
        print('\n',x,round(x**0.5,3))
        print(f'({alpha},{beta})')
        for i in range(alpha):
            print('|=======|'*beta)
            print('|---O---|'*beta)
            print('|_______|'*beta)
        print('\n')
    return alpha,beta

###########################################################################3

def save_sky(skies_obs,
             Fname,
             OBJ,
             RA,DEC,
             MJD,
             Ofile,
             flux,wave,
             GLAT,GLON,HTEL,
             path ='./',
             function_used_to_merge=np.median,
             AXIS_TO_MERGE=0,
             wave_to_angstrom=1):

    os.makedirs(f'{path}', exist_ok=True)

    # Create a FITS header
    header = fits.Header()
    header['OBJECT']  = (OBJ,   'Type of file')

    # Info about the star
    
    date=Time(MJD,format='mjd')
    date=date.isot
    header['MJD']     = (MJD,   'MJD | ' + str(date))
  
    #Info about the parent files
    header['OFIL']   = (Ofile,f'Original parent File')
 

    # Info about telescope
    header['GEOLON']  = (GLAT,  'Geo Latitude of Telescope')
    header['GEOLAT']  = (GLON,  'Geo Longitude of Telescope')
    header['GEOLEV']  = (HTEL,  'Height over sea of Telescope')

    # Comments Example
    header['COMMENT'] = "File done by I. Baeza"

    # Create primary hdu
    primary_hdu = fits.PrimaryHDU(header=header)

    FIBERS=[primary_hdu]

    for nnn,fiber in enumerate(skies_obs):
        temp_header = fits.Header()
        
        wave_samp=np.linspace(min(wave),max(wave),len(fiber))
        fiber =interp1d(wave,fiber)(wave_samp)
        wave_samp*=wave_to_angstrom

        temp_header['RA']      = (RA[nnn],    'Right Ascension /'+str(Angle(float(RA[nnn])*u.degree).to('hourangle')))
        temp_header['DEC']     = (DEC[nnn],   'Declination     /'+str(Angle(float(DEC[nnn])*u.degree).to('hourangle')))
        temp_header['EXTNAME'] = f'SKY-{nnn}'
        temp_header['CRVAL1']  = (wave_samp[0], 'wavelenght first pixel ')
        temp_header['CDELT1']  = (wave_samp[1]-wave_samp[0], 'Step of wavelenght ')
        FIBERS.append(fits.ImageHDU(header=temp_header,data=fiber))


    # Create the columns of data
    columns=[
        fits.Column(name='wave', format='D', unit='Angstrom',array= wave_samp ),
        fits.Column(name='counts',format='D', unit='Relative',array= function_used_to_merge(skies_obs,axis=AXIS_TO_MERGE) )
        ]
    # Generate the secondary fits with the data
    NEWFITS = fits.BinTableHDU.from_columns(columns)

    # Set the EXTNAME keyword to label the binary table HDU
    NEWFITS.header['EXTNAME'] = 'Merged Sky'
    # Comments about the columns (Topcat)
    NEWFITS.header['TCOMM1']  = 'wavelenght'
    NEWFITS.header['TCOMM2']  = 'Flux' 

    FIBERS.append(NEWFITS)
    # Save the fits file
    hdul = fits.HDUList(FIBERS)
    hdul.writeto(path+Fname, overwrite=True)

#######################################################################
def extend_linear_range(min_val, max_val, num_points, rv_kms=50.0):
    """
    Extends a range by a given rv value while maintaining 
    the original point density (proportional spacing).
    """
    # 1. Calculate the original step size (delta)
    # Spacing = (max - min) / (n - 1)
    original_range = max_val - min_val
    delta = original_range / (num_points - 1)
    
    # 2. Define the new boundaries
    new_min = min_val*(1-rv_kms/299792.458)
    new_max = max_val*(1+rv_kms/299792.458)
    
    # 3. Calculate the new number of points to keep density constant
    # New N = (New Range / delta) + 1
    new_range = new_max - new_min
    new_num_points = int(round(new_range / delta)) + 1
    
    # 4. Generate the new array
    return np.linspace(new_min, new_max, new_num_points)