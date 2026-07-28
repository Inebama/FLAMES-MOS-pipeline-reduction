# Packages 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.io import fits
from scipy.interpolate import interp1d  
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from sklearn.neighbors import KernelDensity
import astropy.units as U
import glob
from scipy.optimize import curve_fit    # To get errors


import sys 

sys.path.append("/home/ian/Desktop/Codigos/Dani_Codes/")

from gls import Gls

np.random.seed(170499)

###########################################################
def gaussian(x, x0, sigma, a=1):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
    
###########################################################
def Folding(t,          # Data
            width,      # box size of the folding
            ref=None    # Reference time
            ):
    #------------------------
    if type(ref)==type(None):  # Reference time to start the folding
        ref=np.nanmin(t)       # if not set we use the minimum available

    t=(((t-ref)%width)/width)
    t[t<0]+=1
    return t

#  Returns tuple of handles, labels for axis ax, after reordering them to conform to the label order `order`, and if unique is True, after removing entries with duplicate labels.
def reorderLegend(ax=None,order=None,unique=False):
    if ax is None: ax=plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0])) # sort both labels and handles by labels
    if order is not None: # Sort according to a given list (not necessarily complete)
        keys=dict(zip(order,range(len(order))))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
    if unique:  labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels)) # Keep only the first of each handle
    ax.legend(handles, labels)
    return(handles, labels)
#________________________________

PATH='/home/ian/Desktop/MY_OUTPUT/GIRAFFE/'
LIST=['OB-1','OB-2','OB-3','OB-4','OB-6','OB-7','OB-8','OB-9','OB-10']#
order_folder='/RESTFRAME/Continuum_sub/'


# PATH='/home/ian/Desktop/MY_OUTPUT/'
# LIST=['OB-1','OB-2','OB-3','OB-4','OB-6','OB-7','OB-8','OB-9','OB-10']#
# order_folder='/RESTFRAME/'



TARGET='MED_RGB_1743*'
TARGET='MED_VAR_581*'
TARGET='MED_VAR_3405*'
TARGET='MED_RGB_953*' # pcygny profile




RV,RVE=[],[]
MJD=[]
for epoch in LIST:

    file=glob.glob(PATH+epoch+order_folder+TARGET)
    with fits.open(file[0]) as hdul:
        MJD.append(hdul[0].header['MJD'])


    ra,dec,rv,rve=np.genfromtxt(PATH+epoch+'/OUTPUT_RV.txt',unpack=True,usecols=(1,2,-2,4))
    # check=rve>np.median(rve)+2*np.std(rve)

    # print(epoch, ra[check],dec[check])
    if epoch==LIST[0]:
        OBJ=np.genfromtxt(PATH+epoch+'/OUTPUT_RV.txt',unpack=True,usecols=(0),dtype=('str'))
        SK1=SkyCoord(ra=ra*U.degree,dec=dec*U.degree)
        RV.append(rv)
        RVE.append(rve)
    else:
        SK2=SkyCoord(ra=ra*U.degree,dec=dec*U.degree)
        indexes,d2d,d3d=SK1.match_to_catalog_sky(SK2)
        
        if epoch=='OB-2':
            rv =rv[indexes]
            rve=rve[indexes]

            mask= ~((OBJ!='MED_RGB_370')&(OBJ!='MED_RGB_402'))

            rv[mask] =np.nan
            rve[mask]=np.nan
            RV.append(rv)
            RVE.append(rve)
            continue

        RV.append(rv[indexes])
        RVE.append(rve[indexes])


err=np.nanmedian(RVE) # Median error RV in our sample

Full_RV =np.array(RV)
Full_err=np.array(RVE)

RVE=(np.nanmean(RVE,axis=0)**2+np.nanstd(RV,axis=0)**2)**0.5
RV = np.nanmean(RV,axis=0)

print('\n\n')
print('Stats RV')
print('Mean\tMedian\tSTD\tMedian(ERR)\tMedian(STD)')
print(np.mean(RV[1:]),np.median(RV[1:]),np.std(RV[1:]),np.median(RVE[1:]),np.median(np.nanstd(Full_RV,axis=0)/9))
print('\n\n')





check=RVE>=1.3#err+5*np.median(np.nanstd(Full_RV,axis=0)/9)
inDEX=np.arange(len(RV))[check]



for II,i,j,k,l in zip(inDEX,OBJ[check],SK1.ra[check],SK1.dec[check],RVE[check]):
    # print(II,i,Angle(j).to_string(unit='hourangle', sep=':'),Angle(k).to_string(unit='deg', sep=':'),l)
    print(II,i,j.value,k.value,l)

# sys.exit()



fig = plt.figure(figsize=(15, 6))
gs  = gridspec.GridSpec(len(inDEX), 3,
                    height_ratios=[1]*len(inDEX), 
                    width_ratios =[1,0.7,0.7],
                    left   = 0.12, 
                    right  = 0.93, 
                    bottom = 0.1, 
                    top    = 0.9, 
                    wspace = 0.05, 
                    hspace = 0.02)

# Axes for spectrum plot, sliders, and spectra list
ax1 = fig.add_subplot(gs[:,0])  
ax1.set_title('RV range of Variables')

for j,i in enumerate(inDEX):
    if j==0:
        ax = fig.add_subplot(gs[j,1])  
        ax2= fig.add_subplot(gs[j,2])  

        ax.set_title('Time Series')
        ax2.set_title('Phase Fold')
    else:
        ax = fig.add_subplot(gs[j,1],sharex=ax)   
        ax2= fig.add_subplot(gs[j,2],sharex=ax2)   

    temp_rv =Full_RV[:,i]
    temp_rve=Full_err[:,i]
    center= np.median(temp_rv)

    if OBJ[i]== 'MED_RGB_663':
        GGG=Gls((MJD,temp_rv,temp_rve), Pbeg=800, Pend=1000, ofac=10000,)
    else:
        GGG=Gls((MJD,temp_rv,temp_rve), Pbeg=30, Pend=1000, ofac=10000,)


    print("\n\n\n")
    print(OBJ[i])
    GGG.info()

    Period= 1./GGG.hpstat["f"]
    T0    = GGG.hpstat["T0"]#-Period/2
    AMP   = GGG.hpstat["amp"]

     


    ax1.errorbar(temp_rv,np.zeros(len(temp_rv))+len(inDEX)-1-j,xerr=temp_rve,elinewidth=2,linewidth=0,marker='o',label=OBJ[i],markersize=6,ecolor='gray')

    X  = np.linspace(min(MJD),max(MJD),1000)
    SIN= np.sin(2*np.pi*(X-T0)/Period)*AMP
    ax.plot(X,SIN)

    ax.errorbar(MJD,temp_rv-center,yerr=temp_rve,elinewidth=1,linewidth=0,marker='o',label=OBJ[i],color='k',markersize=5)
    # ax.legend(loc='lower right',bbox_to_anchor=(1,1))
    # ax.set_ylim(-AMP*1.1,AMP*1.1)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)

    X  = np.linspace(0,1,100)
    SIN= np.sin(2*np.pi*X)*AMP
    ax2.plot(X,SIN)

    ax2.legend(title=f'Period\n{Period:.2f} days',loc='upper left',bbox_to_anchor=(1,1))

    fold_time=Folding(MJD,Period,T0)
    
    ax2.errorbar(fold_time,temp_rv-center,yerr=temp_rve,elinewidth=1,linewidth=0,marker='o',label=OBJ[i],color='k',markersize=5)
    

    ax2.yaxis.set_visible(False)
    ax2.xaxis.set_visible(False)


    
ax.xaxis.set_visible(True)
ax2.xaxis.set_visible(True)  

ax1.set_yticks(np.arange(len(inDEX))[::-1],labels=OBJ[inDEX],fontsize=13)
# plt.legend(loc='upper left',bbox_to_anchor=(1,1))
ax1.set_xlabel('RV [km/s]')
ax.set_xlabel('MJD')
ax2.set_xlabel('Phase')
ax1.grid(ls=':',color='gray')

plt.tight_layout()


# sys.exit()

fig = plt.figure(figsize=(10, 6))
gs  = gridspec.GridSpec(2, 1,
                    height_ratios=[1,0.2], 
                    width_ratios =[1],
                    left   = 0.09, 
                    right  = 0.97, 
                    bottom = 0.1, 
                    top    = 0.9, 
                    wspace = 0.0, 
                    hspace = 0.0)

# Axes for spectrum plot, sliders, and spectra list
ax1 = fig.add_subplot(gs[0])   # Large plot for spectrum and continuum
ax2 = fig.add_subplot(gs[1],sharex=ax1)   # Large plot for spectrum and continuum

ax2.yaxis.set_visible(False)

from scipy.stats import norm

def variable_bandwidth_kde(X, bandwidths, query_points, kernel="gaussian"):
    """
    Perform KDE with variable bandwidths for each data point and optional kernel type.

    Parameters:
        X (array-like): 1D array of data points (n_samples,).
        bandwidths (array-like): 1D array of individual bandwidths for each data point (n_samples,).
        query_points (array-like): 1D array of points where to evaluate the KDE (m_samples,).
        kernel (str): The type of kernel to use ('gaussian' or 'tophat').

    Returns:
        array-like: 1D array of density estimates at query points (m_samples,).
    """
    X = np.asarray(X).ravel()
    bandwidths = np.asarray(bandwidths).ravel()
    query_points = np.asarray(query_points).ravel()
    densities = np.zeros(query_points.shape)

    for i, qp in enumerate(query_points):
        if kernel == "gaussian":
            # Gaussian kernel
            kernels = norm.pdf((qp - X) / bandwidths) / bandwidths
        elif kernel == "tophat":
            # Top-hat kernel: Uniform density within [X - h, X + h]
            kernels = np.where(np.abs(qp - X) <= bandwidths, 1 / (2 * bandwidths), 0)
        else:
            raise ValueError("Unsupported kernel type. Use 'gaussian' or 'tophat'.")
        
        densities[i] = np.sum(kernels) / X.size

    return densities




HIST=ax1.hist(RV,bins=40,color='gray',edgecolor='None',density=True,)
WIDTH=np.mean(HIST[1][1:]-HIST[1][:-1])

HIST[-1].set_label(f'Histogram RV (bin size ={WIDTH:.0f})')

new_x=(HIST[1][1:]+HIST[1][:-1])/2

popt, pcov    = curve_fit(gaussian,new_x, HIST[0], p0=[np.median(RV),np.std(RV),max(HIST[0])])

print('\n\n')
print('Gaussian fit')
print(popt)
print('\n\n')

X   = np.linspace(-51,75,100)
Y   = variable_bandwidth_kde(RV, bandwidths=RVE,query_points=X,kernel='gaussian')
ax1.plot(X,Y,label='Gaussian Kernel Density (Bandwidth = errors measures)')

kde = KernelDensity(kernel='epanechnikov', bandwidth=5*np.mean(RVE)).fit(RV.reshape(-1,1))
Y   = kde.score_samples(X.reshape(-1,1))


ax1.plot(X,gaussian(X,*popt),color='k',label='Gaussian Fit to Histogram Values')

ax1.plot(X,np.exp(Y.ravel()),label='Epanechnikov Kernel Density (Bandwidth = 5 x Median Error)')
ax1.set_ylabel('Density',fontsize=12)
ax2.set_xlabel('RV [km/s]',fontsize=12)

ax1.axvline(np.median(RV),color='g',label='Median RV')
ax1.axvline(np.mean(RV),color='r',label='Mean RV')

ax1.errorbar([],[],xerr=[],elinewidth=2,linewidth=0,marker='|',color='k',markersize=15,label='Data Measured')
ax1.legend(loc='upper right')


reorderLegend(ax1,[f'Histogram RV (bin size ={WIDTH:.0f})','Gaussian Fit to Histogram Values','Gaussian Kernel Density (Bandwidth = errors measures)','Epanechnikov Kernel Density (Bandwidth = 5 x Median Error)', 'Median RV','Mean RV', 'Data Measured'])


ax2.errorbar(RV,np.random.normal(0,0.3,len(RV)),xerr=RVE,elinewidth=2,linewidth=0,marker='|',color='k',markersize=15)
ax2.set_ylim(-1,1)# plt.show()

plt.show()

'''
 Lebzelter 2005 


0 17  MED_VAR_3094 0:23:51.26 -72:03:49   -> LW4
1 20  MED_RGB_85   0:24:03.39 -71:55:47.4 -> ?
2 22  MED_RGB_2690 0:24:30.16 -72:04:31.4 -> ?
3 47  MED_RGB_3396 0:25:23.29 -72:02:37.5 -> ?
4 63  MED_VAR_1142 0:25:03.68 -72:09:31.7 -> V5
5 74  MED_VAR_1163 0:24:29.57 -72:09:07.6 -> V23
6 82  MED_RGB_2430 0:24:21.74 -72:04:48.2 -> A19
7 95  MED_VAR_581  0:22:58.44 -72:06:56.2 -> V13
8 105 MED_RGB_663  0:22:12.51 -72:05:45   -> ?
9 106 MED_RGB_629  0:22:17.9  -72:06:12.3 -> ?


'''


# with open('RV_SUMMARY.txt','w') as f:
#     f.write("#OBJ\tRA\tDEC\tRV\tRV_ERR\n")
#     for i,j,k,l,m in zip(OBJ,SK1.ra.value,SK1.dec.value,RV,RVE):
#         f.write(f"{i}\t{j}\t{k}\t{l}\t{m}\n")
