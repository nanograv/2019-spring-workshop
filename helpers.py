import numpy as np
import mpmath as mp
from pypulse.singlepulse import SinglePulse
from pypulse.archive import Archive
import pypulse.utils as u
from matplotlib.pyplot import *
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
from IPython.display import display, Math, Latex
from astropy.coordinates import SkyCoord
#from PAL2 import bayesutils as bu


# From startup.py
'''
import sys
sys.path.append("/home/jupyter/shared/")
from nanograv_data.api.nanograv_apy import NANOGravAPI
api = NANOGravAPI()
from nanograv_utils import NANOGravUtils
utils = NANOGravUtils()
'''

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'],'size':14})

mp.dps = 300

shiftit = u.shiftit


# For dispersion plots
def getTemplate(PSR,rcvr):
    if rcvr == "Rcvr1_2" or rcvr == "Rcvr_800":
        dirname = "guppi2"
        be = "GUPPI"
    elif rcvr == "L-wide" or rcvr == "S0wide" or rcvr == "430" or rcvr == "327":
        dirname = "puppi"
        be = "PUPPI"
    tempfilename = "/nanograv/releases/11y/%s/%s/%s.%s.%s.11y.x.sum.sm"%(dirname,PSR,PSR,rcvr,be)
    ar = Archive(tempfilename,verbose=False)
    sptemp = SinglePulse(u.normalize(ar.getData(),simple=True),windowsize=256,period=ar.getPeriod())
    return sptemp



def dispersionPlot(num,figsize=(10,8)):
    sptemp = getTemplate("J1713+0747","Rcvr1_2")
    sptemp.shiftit(-512,save=True)
    bins = np.array(sptemp.bins,dtype=np.float)
    data = sptemp.data
    P0 = sptemp.getPeriod()
    dt = 1e3 * P0/len(bins)
    bins *= dt
    bins -= (1e3*P0)/4
    fig = figure(figsize=figsize)
    if num == 1:
        shiftdata = sptemp.shiftit(0.4788/dt) #time in ms for difference between 1.40 and 1.41 GHz
        ax = fig.add_subplot(211)
        ax.plot(bins,data,'k')
        ax.set_xlim(min(bins),max(bins))
        ax.set_ylim(-0.1,1.1)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        #ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params('both',length=8,which='major')
        ax.tick_params('both',length=4,which='minor')
        ax.text(2.25,0.75,r"$\mathrm{J1713+0747}$"+"\n"+r"$\mathrm{\nu=1.40\;GHz}$",size=16)

        ax = fig.add_subplot(212)
        ax.plot(bins,shiftdata,'k')
        ax.set_xlim(min(bins),max(bins))
        ax.set_ylim(-0.1,1.1)
        ax.set_yticks([])
        ax.text(2.25,0.75,r"$\mathrm{J1713+0747}$"+"\n"+r"$\mathrm{\nu=1.39\;GHz}$",size=16)
        ax.set_xlabel(r'$\mathrm{Time\;(ms)}$')
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params('both',length=8,which='major')
        ax.tick_params('both',length=4,which='minor')
    elif num == 2 or num == 3:
        ax = fig.add_subplot(211)
        ax.plot(bins,data,'k')
        ax.set_xlim(min(bins),max(bins))
        ax.set_ylim(-0.10,1.1)
        ax.set_yticks([])
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params('both',length=8,which='major')
        ax.tick_params('both',length=4,which='minor')
        #ax.set_xticks([])
        ax.text(2.25,0.75,r"$\mathrm{J1713+0747}$"+"\n"+r"$\mathrm{\nu=1.40\;GHz}$",size=16)

        ax = fig.add_subplot(212)
        if num == 2:
            T = 2917.2336 #DM = 16.00
        elif num == 3:
            T = 2915.410 #DM = 15.99
        shiftdata = sptemp.shiftit(T/dt) #time in ms for difference between 1.40 and 0.150 GHz
        bins2 = bins + T
        ax.plot(bins2,shiftdata,'k')
        ax.set_xlim(min(bins2),max(bins2))
        ax.set_ylim(-0.10,1.1)
        ax.set_yticks([])
        ax.set_xlabel(r'$\mathrm{Time\;(ms)}$')
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params('both',length=8,which='major')
        ax.tick_params('both',length=4,which='minor')
        ax.text(2.25+T,0.75,r"$\mathrm{J1713+0747}$"+"\n"+r"$\mathrm{\nu=0.15\;GHz}$",size=16)

                       
            
            
# For Profile Evolution
def load_profile(filename):
    ar = Archive(filename,verbose=False)
    return SinglePulse(ar.getData(),windowsize=256)


            
            

def pulse(x,sigma,area=True):
    if area:
        amplitude = 1.0/np.sqrt(2*np.pi*sigma**2)
    else:
        amplitude = 1.0
    return amplitude*np.exp(-0.5*(x/sigma)**2)

def scattered_pulse(x,sigma,taud,area=True,log=True):
    if area:
        amplitude = (1.0/taud)*(1.0/np.sqrt(2*np.pi*sigma**2)) #convolution of gaussian of unit area and exponential of unit area        
    else:
        amplitude = 1.0#convolution of gaussian of unit amplitude and exponential of unit peak

    xbar = sigma/(taud*np.sqrt(2)) - x/(sigma*np.sqrt(2))


    if log: #compute in the natural log
        retval = np.zeros_like(xbar)
        for i in range(len(x)):
            retval[i] = mp.exp(mp.log(amplitude * (sigma/2)*np.sqrt(np.pi/2)) + 0.5*sigma**2/taud**2 -x[i]/taud + float(mp.log(mp.erfc(xbar[i]))))
        return retval*2#factor of 2?
    else:
        return amplitude * (sigma/2)*np.sqrt(np.pi/2) * np.exp(0.5*sigma**2/taud**2) * np.exp(-x/taud) * (1-special.erf(xbar))
    
def add_noise(prof,SN):
    return prof + np.random.normal(0,np.ptp(prof)/SN,len(prof))
    

def fit_plot(temp,prof,verbose=True,preroll=True):
    if isinstance(temp,SinglePulse):
        sptemp = temp
    else:
        sptemp = SinglePulse(temp,windowsize=256)
    rollval = sptemp.getNbin()/2 - np.argmax(sptemp.data)     
    if isinstance(prof,SinglePulse):
        spprof = prof
    else:
        spprof = SinglePulse(prof,mpw=sptemp.mpw)
    if preroll:
        sptemp.shiftit(rollval,save=True)
        spprof.shiftit(rollval,save=True)
    tauccf, tauhat, bhat, sigma_tau, sigma_b, snr, rho = spprof.fitPulse(sptemp.data)
    if verbose:
        display(Math(r'\hat{\tau} = %0.1f \pm %0.1f'%(tauhat,sigma_tau)))
        display(Math(r'\hat{b} = %0.1f \pm %0.1f'%(bhat,sigma_b)))
    if preroll:
        sptemp.shiftit(-1*rollval,save=True)
        spprof.shiftit(-1*rollval,save=True)
    rotatedtemp = sptemp.shiftit(tauhat)
    figure(figsize=(10,8))
    subplot(211)
    plot(prof.data,'k')
    plot(bhat*rotatedtemp,'r',lw=2)
    xlim(0,len(prof.data))
    ylabel("Profile")
    subplot(212)
    plot(prof.data-bhat*rotatedtemp,'k')
    xlim(0,len(prof.data))
    xlabel("Phase Bins")
    ylabel("Difference Profile")
    return tauhat,bhat



def bu_triplot(chain, color='k', weights=None, interpolate=False, smooth=True, \
           labels=None, figsize=(11,8.5), title=None, inj=None, tex=True, \
            incMaxPost=True, cmap='YlOrBr', lw=1.5, ranges=False, axarr=None):

    """
    Make Triangle plot
    """

    # rcParams settings
    if chain.shape[1] < 10:
        ticksize = 10
        #plt.rcParams['ytick.labelsize'] = 10.0
        #plt.rcParams['xtick.labelsize'] = 10.0
    else:
        ticksize = 8
        #plt.rcParams['ytick.labelsize'] = 8.0
        #plt.rcParams['xtick.labelsize'] = 8.0
    if tex:
        plt.rcParams['text.usetex'] = True


    # get number of parameters
    ndim = chain.shape[1]
    parameters = np.arange(ndim, dtype=np.int)
    
    if axarr is not None:
        f = gcf()
        #fig, axarr = plt.subplots(nrows=len(parameters), ncols=len(parameters),figsize=figsize)
    else:
        f, axarr = subplots(nrows=len(parameters), ncols=len(parameters),figsize=figsize)

    for i in range(len(parameters)):
        # for j in len(parameters[np.where(i <= parameters)]:
        for j in range(len(parameters)):
            ii = i
            jj = len(parameters) - j - 1

            # get ranges
            if ranges:
                xmin, xmax = confinterval(chain[:, parameters[ii]], sigma=0.95, 
                                          type='equalProb')
                x_range = [xmin, xmax]
                xmin, xmax = confinterval(chain[:, parameters[jj]], sigma=0.95, 
                                          type='equalProb')
                y_range = [xmin, xmax]

            else:
                x_range = [chain[:, parameters[ii]].min(), chain[:, parameters[ii]].max()]
                y_range = [chain[:, parameters[jj]].min(), chain[:, parameters[jj]].max()]


            axarr[ii, jj].tick_params(axis='both', which='major', labelsize=10)

            xmajorLocator = matplotlib.ticker.MaxNLocator(nbins=4,prune='both')
            ymajorLocator = matplotlib.ticker.MaxNLocator(nbins=4,prune='both')

            if j <= len(parameters)-i-1:
                axarr[jj][ii].xaxis.set_minor_locator(NullLocator())
                axarr[jj][ii].yaxis.set_minor_locator(NullLocator())
                axarr[jj][ii].xaxis.set_major_locator(NullLocator())
                axarr[jj][ii].yaxis.set_major_locator(NullLocator())

                axarr[jj][ii].xaxis.set_minor_formatter(NullFormatter())
                axarr[jj][ii].yaxis.set_minor_formatter(NullFormatter())
                axarr[jj][ii].xaxis.set_major_formatter(NullFormatter())
                axarr[jj][ii].yaxis.set_major_formatter(NullFormatter())
                xmajorFormatter = FormatStrFormatter('%g')
                ymajorFormatter = FormatStrFormatter('%g')

                if ii == jj:
                    # Make a 1D plot
                    bu.makesubplot1d(axarr[ii][ii], chain[:,parameters[ii]], \
                                  weights=weights, interpolate=interpolate, \
                                  smooth=smooth, color=color, lw=lw, range=x_range)
                    axarr[ii][jj].set_ylim(ymin=0)
                    if incMaxPost:
                        mx = bu.getMax(chain[:,parameters[ii]], weights=weights)
                        axarr[ii][jj].set_title('%5.4g'%(mx), fontsize=10)

                    if inj is not None:
                        axarr[ii][ii].axvline(inj[ii], lw=2, color='k')
                else:
                    # Make a 2D plot
                    bu.makesubplot2d(axarr[jj][ii], chain[:,parameters[ii]],
                                  chain[:,parameters[jj]], cmap=cmap, 
                                  color=color, weights=weights,
                                  smooth=smooth, lw=lw, x_range=x_range,
                                  y_range=y_range)

                    if inj is not None:
                        axarr[jj][ii].plot(inj[ii], inj[jj], 'x', color='k', markersize=12, \
                                           mew=2, mec='k')

                axarr[jj][ii].xaxis.set_major_locator(xmajorLocator)
                axarr[jj][ii].yaxis.set_major_locator(ymajorLocator)
            else:
                axarr[jj][ii].set_visible(False)
                #axarr[jj][ii].axis('off')

            if jj == len(parameters)-1:
                axarr[jj][ii].xaxis.set_major_formatter(xmajorFormatter)
                if labels:
                    axarr[jj][ii].set_xlabel(labels[ii])

            if ii == 0:
                if jj == 0:
                    axarr[jj][ii].yaxis.set_major_locator(NullLocator())
                    #axarr[jj][ii].set_ylabel('Post.')
                else:
                    axarr[jj][ii].yaxis.set_major_formatter(ymajorFormatter)
                    if labels:
                        axarr[jj][ii].set_ylabel(labels[jj])

    # overall plot title
    if title:
        f.suptitle(title, fontsize=14, y=0.90)
     
    # make plots closer together 
    f.subplots_adjust(hspace=0.1)
    f.subplots_adjust(wspace=0.1)

    return axarr



# CW Related functions

def make_fake_pulsar(DIR):
    '''
    Makes a fake pulsar par file
    '''
    output = "MODE 1\n"
    
    # Sphere Point Picking
    u = np.random.uniform()
    v = np.random.uniform()
    phi = 2*np.pi*u #using standard physics notation
    theta = np.arccos(2*v-1) - np.pi/2

    c = SkyCoord(phi,theta,frame='icrs',unit='rad')
    cstr = c.to_string('hmsdms')
    #print cstr
    RAJ = cstr.split(" ")[0].replace("h",":").replace("m",":")[:-1]
    DECJ = cstr.split(" ")[1].replace("d",":").replace("m",":")[:-1]
    cstr = cstr.replace(" ","")
    name = "J"+RAJ[0:2]+RAJ[3:5]+DECJ[0]+DECJ[1:3]+DECJ[4:6]

    output += "PSR      %s\n"%name

    
    output += "PEPOCH   50000.0\n"    
    output += "POSEPOCH   50000.0\n"

    period = 0.001*np.random.uniform(1,10) #seconds
    output += "F0       %0.10f 1\n"%(1.0/period)

    output += "RAJ      %s 1\n"%RAJ
    output += "DECJ     %s 1\n"%DECJ

    dist = np.random.uniform(0.1,5) #kpc
    output += "PX       %0.5f 1\n"%(1.0/dist)

    filename = "%s/%s.par"%(DIR,name)
    with open(filename,'w') as FILE:
        FILE.write(output)

    return filename.encode('ascii','ignore')

    #pass
    


#if __name__ == '__main__':
#    print make_fake_pulsar("MTL")
    
