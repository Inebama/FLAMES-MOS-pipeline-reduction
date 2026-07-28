# For continuum sharing same parameters across observations (copying the Stored params on each folder wiht the respective names)
if 0: 
    import numpy as np 

    path='/scratch/home/ibaeza/Desktop/MY_OUTPUT/GIRAFFE_median_1h_new/'

    data_or=np.load(path+'OB-15/RESTFRAME/Stored_params.npy',allow_pickle=True).item()

    Folders=['OB-1','OB-2','OB-3','OB-4','OB-5','OB-6','OB-7','OB-8','OB-9','OB-10','OB-11','OB-12','OB-13','OB-14','OB-16']

    names=dict()
    for i in data_or.keys():
        names[i.split('=')[0]]=i


    # path='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_mean/'


    for i in Folders:
        spectrums=np.genfromtxt(path+i+'/Oficial_list.list',dtype='U100')
        params = dict()
        for j in spectrums:
            params[j]=data_or[names[j.split('=')[0]]]
        # BEWARE of overwriting this
        np.save(path+i+'/RESTFRAME/Stored_params.npy',params)

# All files per target
if 0:
    import numpy as np
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as U
    # PATH='/home/ian/Desktop/MY_OUTPUT/GIRAFFE/'
    # PATH='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_mean_48h/'
    # PATH='/home/ian/Desktop/MY_OUTPUT/GIRAFFE_median_48h_new/'
    PATH='/scratch/home/ibaeza/Desktop/MY_OUTPUT/GIRAFFE_median_15min_new/'
    # targets='/home/ian/Downloads/TARGETS_all_photometry_rvs.fits'
    targets='/scratch/home/ibaeza/Desktop/DATA/TARGETS_all_photometry_rvs.fits'
    with fits.open(targets) as F:
        I    =F[1].data.I
        V    =F[1].data.V
        sort = np.argsort(I+V)[::-1]
        V,I  = V[sort],I[sort]
        OBJ  =F[1].data.Object[sort]
        Prob =F[1].data.Prob[sort]
        ra   =F[1].data.RA[sort]
        dec  =F[1].data.DEC[sort]
        OBS_mul=F[1].data.Multiple_Obs[sort]

        SC1=SkyCoord(ra*U.degree,dec*U.degree)
    
    # LIST=['OB-0','OB-1','OB-2','OB-3','OB-4','OB-6','OB-7','OB-8','OB-9','OB-10']#
    DIST=1.75*U.arcsec

    # # Separated by 1.04 hours
    # SUPER_LIST={'OB-1':['OB-1','OB-2','OB-3'],
    #             'OB-2':['OB-4','OB-5','OB-6'],
    #             'OB-3':['OB-7','OB-8','OB-9'],
    #             'OB-4':['OB-10','OB-11','OB-12'],
    #             'OB-5':['OB-13','OB-14','OB-15'],
    #             'OB-7':['OB-17','OB-18','OB-19','OB-20','OB-21'],
    #             'OB-8':['OB-22','OB-23','OB-24','OB-25','OB-26'],
    #             'OB-9':['OB-27','OB-28','OB-29','OB-30','OB-31'],
    #             'OB-10':['OB-32','OB-33','OB-34','OB-35','OB-36'],
    #             'OB-11':['OB-37','OB-38','OB-39','OB-40','OB-41'],
    #             'OB-12':['OB-42','OB-43','OB-44','OB-45','OB-46'],
    #             'OB-13':['OB-47','OB-48','OB-49','OB-50','OB-51'],
    #             'OB-14':['OB-52','OB-53','OB-54','OB-55'],
    #             'OB-15':['OB-56','OB-57','OB-58','OB-59'],
    #             'OB-16':['OB-60','OB-61','OB-62','OB-63']}
    
    # # Separated by 48 hours
    # SUPER_LIST={'OB-1':['OB-1','OB-2','OB-3','OB-4','OB-5','OB-6'],
    #             'OB-2':['OB-7','OB-8','OB-9'],
    #             'OB-3':['OB-10','OB-11','OB-12'],
    #             'OB-4':['OB-13','OB-14','OB-15'],
    #             'OB-6':['OB-17','OB-18','OB-19','OB-20','OB-21'],
    #             'OB-7':['OB-22','OB-23','OB-24','OB-25','OB-26','OB-27','OB-28','OB-29','OB-30','OB-31'],
    #             'OB-8':['OB-32','OB-33','OB-34','OB-35','OB-36','OB-37','OB-38','OB-39','OB-40','OB-41','OB-42','OB-43','OB-44','OB-45','OB-46','OB-47','OB-48','OB-49','OB-50','OB-51'],
    #             'OB-9':['OB-52','OB-53','OB-54','OB-55'],
    #             'OB-10':['OB-56','OB-57','OB-58','OB-59','OB-60','OB-61','OB-62','OB-63']}
    
    # All Merged
    SUPER_LIST={'OB-1':['OB-1','OB-2','OB-3','OB-4','OB-5','OB-6','OB-7','OB-8','OB-9','OB-10','OB-11','OB-12','OB-13','OB-14','OB-15','OB-17','OB-18','OB-19','OB-20','OB-21','OB-22','OB-23','OB-24','OB-25','OB-26','OB-27','OB-28','OB-29','OB-30','OB-31','OB-32','OB-33','OB-34','OB-35','OB-36','OB-37','OB-38','OB-39','OB-40','OB-41','OB-42','OB-43','OB-44','OB-45','OB-46','OB-47','OB-48','OB-49','OB-50','OB-51','OB-52','OB-53','OB-54','OB-55','OB-56','OB-57','OB-58','OB-59','OB-60','OB-61','OB-62','OB-63']}

    for k in SUPER_LIST.keys():
        LIST = SUPER_LIST[k]
        AUX  =[]
        for i in range(len(OBJ)):
            AUX.append({i:None for i in LIST})

        for epoch in LIST:
            files=np.genfromtxt(PATH+epoch+'/INFO_FILES.txt',usecols=(0),dtype=str)
            ra_ep,dec_ep,mag,counts=np.genfromtxt(PATH+epoch+'/INFO_FILES.txt',usecols=(1,2,3,4),unpack=True)
            SC = SkyCoord(ra_ep*U.degree,dec_ep*U.degree)

            idx,d2d,d3d= SC.match_to_catalog_sky(SC1)
            idx = idx[d2d<DIST]
            for i,j in enumerate(idx):
                AUX[j][epoch]=files[i]

        Final=dict()
        for i, star in enumerate(OBJ):
            Final[star]=AUX[i]
            
        for key in Final.keys():
            print(key,Final[key])

        np.save(PATH+'MERGED_all/'+k+'Targets_files.npy',Final)

# Test broadening
if 0:
    import numpy as np
    from scipy.signal import fftconvolve
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt

    def rotational_profile(v, vsini, epsilon=0.6):
        x = v / vsini
        profile = np.zeros_like(x)
        inside = np.abs(x) <= 1
        x_inside = x[inside]
        profile[inside] = ((2 * (1 - epsilon) * np.sqrt(1 - x_inside**2) +
                            np.pi * epsilon * (1 - x_inside**2) / 2) /
                        (np.pi * vsini * (1 - epsilon / 3)))
        return profile / np.sum(profile)

    def apply_rotational_broadening_loglambda(wavelength, flux, vsini, epsilon=0.6, oversample=4):
        c = 299792.458  # km/s

        # Create log-wavelength grid (constant velocity steps)
        loglam = np.log(wavelength)
        dloglam = np.mean(np.diff(loglam))
        dv = dloglam * c
        loglam_fine = np.arange(loglam[0], loglam[-1], dloglam / oversample)
        
        # Interpolate flux onto fine log-lambda grid
        flux_interp = interp1d(loglam, flux, kind='linear', bounds_error=False, fill_value="extrapolate")
        flux_log = flux_interp(loglam_fine)
        
        # Create velocity grid for kernel
        v_extent = 1.2 * vsini  # buffer
        v_grid = np.arange(-v_extent, v_extent, dv / oversample)
        kernel = rotational_profile(v_grid, vsini, epsilon)
        
        # Convolve in velocity space
        flux_conv = fftconvolve(flux_log, kernel, mode='same')

        # Interpolate back to original wavelength grid
        flux_final_interp = interp1d(loglam_fine, flux_conv, kind='linear', bounds_error=False, fill_value="extrapolate")
        return flux_final_interp(loglam)
    

    path_syn_spec  = '/home/ian/Desktop/Synth_spectra/MED_RGB_2430/MED_RGB_2430_R36900.dat'
    path_syn_spec_2= '/home/ian/Desktop/Synth_spectra/MED_RGB_2430/MED_RGB_2430_rot_15_R36900_ROT15.dat'
    
    # path_syn_spec  = '/home/ian/Desktop/Test-Autokur/Co_lines_onlyCo_R25000.dat'
    syn_spec=np.genfromtxt(path_syn_spec,skip_header=47)

    Template_wave = syn_spec[:,0]
    Template_flux = syn_spec[:,1]

    plt.plot(Template_wave,Template_flux,label='0 km/s')
    
    plt.plot(Template_wave, apply_rotational_broadening_loglambda(Template_wave,Template_flux, vsini=15.0),label='15 km/s profile')

    syn_spec=np.genfromtxt(path_syn_spec_2,skip_header=47)

    Template_wave = syn_spec[:,0]
    Template_flux = syn_spec[:,1]

    plt.plot(Template_wave,Template_flux,label='15 km/s AUTOKUR')
    
    # plt.plot(Template_wave, apply_rotational_broadening_loglambda(Template_wave,Template_flux, vsini=15.0)


    # plt.plot(Template_wave, apply_rotational_broadening_loglambda(Template_wave,Template_flux, vsini=100.0),label='100 km/s')

    # plt.plot(Template_wave, apply_rotational_broadening_loglambda(Template_wave,Template_flux, vsini=500.0),label='500 km/s')

    plt.legend()
    plt.show()

# Making median of every target spectra (we can explude the OB-0)
if 0:
    import os
    import numpy as np
    from astropy.io import fits
    from scipy.interpolate import interp1d
    OBS='MERGED_48h/OB-10'
    path_info_targets=f"/scratch/home/ibaeza/Desktop/MY_OUTPUT/GIRAFFE_median_15min_new/{OBS}Targets_files.npy"

    path      ="/scratch/home/ibaeza/Desktop/MY_OUTPUT/GIRAFFE_median_15min_new/"
    # subfolder ='/RESTFRAME/Continuum_sub/'
    subfolder ='/RESTFRAME/'
    new_dir   =f"{OBS}/"
    Files:dict= np.load(path_info_targets,allow_pickle=True).item()


    os.makedirs(path+new_dir, exist_ok=True)
    
    for star in Files.keys():
        STAR=Files[star]
        # Create primary hdu
        primary_hdu = fits.PrimaryHDU()
        hdulist = [primary_hdu]

        # Variables
        wave_min,wave_max=-1,np.inf
        waves,fluxes =[],[]

        # if Hours==48:
        #     if 'OB-2'==j and ((obj=='MED_RGB_370') or (obj=='MED_RGB_402')):
        #         continue
        # if Hours==1.1:
        #     if ('OB-3'==j) and ((obj=='MED_RGB_370') or (obj=='MED_RGB_402')):
        #         continue
        # if Hours==0.1:
        #     if ('OB-9'==j or 'OB-10'==j or 'OB-11'==j) and ((obj=='MED_RGB_370') or (obj=='MED_RGB_402')):
        #         continue

        for OB in STAR.keys():
            if not STAR[OB] is None and OB!='OB-0':
                with fits.open(path+OB+subfolder+STAR[OB]) as hdul:
                    try:
                        wave   = hdul[1].data.wavelength
                        flux   = hdul[1].data.Flux
                        key='Flux'
                    except:
                        try:
                            key='nflux'
                            wave   = hdul[1].data.wave
                            flux   = hdul[1].data.nflux
                        except:
                            key='flux'
                            wave   = hdul[1].data.wave
                            flux   = hdul[1].data.flux
                    wave_min=max([wave_min,min(wave)])
                    wave_max=min([wave_max,max(wave)])

                    waves.append(wave)
                    fluxes.append(flux)

                    wave_samp=np.linspace(min(wave),max(wave),len(flux)*3)
                    flux =interp1d(wave,flux)(wave_samp)

                    hdul[0].header['EXTNAME'] = OB
                    hdul[0].header['CRVAL1']  = (wave_samp[0], 'wavelenght first pixel ')
                    hdul[0].header['CDELT1']  = (wave_samp[1]-wave_samp[0], 'Step of wavelenght ')
                    hdulist.append(fits.ImageHDU(header=hdul[0].header,data=flux))
        
        if not waves:
            continue

        wave_samp=np.linspace(wave_min,wave_max,int(len(flux)*(4/3)))

        FLUXES=[]
        for wave, flux in zip(waves,fluxes):
            flux=interp1d(wave,flux)(wave_samp)
            FLUXES.append(flux)



        # Create the columns of data
        columns=[
            fits.Column(name='wave', format='D', unit='Angstrom',array= wave_samp ),
            fits.Column(name=key,format='D', unit='Relative',array= np.median(FLUXES,axis=0) )
            ]
        # Generate the secondary fits with the data
        NEWFITS = fits.BinTableHDU.from_columns(columns)

        # Set the EXTNAME keyword to label the binary table HDU
        NEWFITS.header['EXTNAME'] = 'SPECTRUM merged'
        # Comments about the columns (Topcat)
        NEWFITS.header['TCOMM1']  = 'wavelenght'
        NEWFITS.header['TCOMM2']  = 'Flux' 

        hdulist.append(NEWFITS)
        # Save the fits file
        hdul = fits.HDUList(hdulist)
        hdul.writeto(path+new_dir+star, overwrite=True)


