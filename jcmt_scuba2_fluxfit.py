#######################################
# JCMT SCUBA-2/POL-2 Analysis Script
#######################################
'''Fits a 2D Gaussian to JCMT SCUBA-2/POL-2 fits images/cubes using an MCMC algorithm, to derive flux densities for compact sources.
Additionally, performs a basic polarization analysis on POL-2 data.

INPUT:  my_dir - raw data directory
        data_dir - output data products directory
        date - date of observation (yyymmdd)
        target - target name
        cal_type - 'user' or 'auto' indicating user calculated FCFs or use of standard values
        obs_type - 'scuba2' or 'pol2'
        waveband - '8' or '4' for 850um or 450um
        integ - lightcurve timebins size in seconds (0 if no lightcurve analysis done)

OUTPUT: Best fit flux/rms for each scan/mosaic map
        Lightcurve files and plot for each timing cube (if integ !=0)
        Diagnostic plots (MCMC traces and hsitograms, best-fits overplotted on all maps/timeslices)

NOTES:  - The beamsize (bmaj, bmin, bpa) is fixed in the fitting, so only flux and position vary. 
        - LSQ algorithm is used to determine initial positions of MCMC walkers in fitting process.
        - For stokes Q and U, the position (x,y) is also fixed to the stokes I best-fit, so only flux is a free parameter.

Written by: Alex J. Tetarenko and Coni Echiburu Trujillo

Last Updated: March 20, 2020
'''

#packages to import
import astropy
from astropy.io import fits,ascii 
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.coordinates import Angle,SkyCoord
from astropy import wcs,coordinates
from astropy import units as u
from astropy.modeling import models,fitting
from astropy.time import Time
from astropy.utils.console import ProgressBar
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
import scipy.stats as ss
import emcee
from multiprocessing import Pool
import os
import glob
from data_reduction import get_scan_numbers
import pickle


def cut_map(filename,pos,types):
    '''Transform an absolute sky position into pixel coordinates
    INPUT:  filename - fits image name you want to transform to (str)
            pos - list of floats of [RA,Dec] (both in units of degrees)
    OUTPUT: list of pixel positions [X,Y]
    '''
    wmap=wcs.WCS(fits.open(filename)[0].header)
    coord=SkyCoord(ra=pos[0]*u.deg,dec=pos[1]*u.deg,frame='icrs')
    if types=='toti':
        x=float(wmap.wcs_world2pix(coord.ra,coord.dec,0,1)[0])
        y=float(wmap.wcs_world2pix(coord.ra,coord.dec,0,1)[1])
    elif types=='pol':
        x=float(wmap.wcs_world2pix(coord.ra,coord.dec,0,0,1)[0])
        y=float(wmap.wcs_world2pix(coord.ra,coord.dec,0,0,1)[1])
    return [x,y]
def confidenceInterval(y):
    '''Get median and 1 sigma confidence limits for a posterior distribution
    INPUT:  y - posterior distribution (array)
    OUTPUT: array of [median, lower bound, upper bound]
    '''
    median=np.median(y)
    pct15=np.percentile(y,15)
    pct85=np.percentile(y,85)
    list1=np.array([median,median-pct15,pct85-median])
    return list1
def lp(p,data,error,fixp,guess,pixsize,types):
    '''Log probability for 2D Gaussian models
    INPUT:  p - parameter array of floats consisting of [amp, x, y, bmaj(fwhm), bmin(fwhm), bpa] (units mJy, pix, pix, arcsec, arcsec, deg)
            data - 2D fits image (array)
            error - 2D fits image variance (array)
            fixp - boolean array indicating which parameters are fixed in fit; same dimensions as p
            guess - array of floats of initial guesses for parameters; same dimensions as p
            pixsize - pixel size in the image in arcsec (float)
            types - 'toti' or 'pol' (str)
    OUTPUT: log probability
    '''
    amp,xx,yy,bmaj,bmin,bpa=p[0],p[1],p[2],p[3],p[4],p[5]
    mod0=models.Gaussian2D(amp,xx,yy,bmaj/(2.*pixsize),bmin/(2.*pixsize),bpa)
    xval=np.arange(0,len(data[0,:]))
    yval=np.arange(0,len(data[:,0]))
    Xval, Yval = np.meshgrid(xval, yval)
    mod1=mod0(Xval,Yval)
    re=-0.5*np.nansum(np.log(2*np.pi*error**2))-np.nansum((mod1-data)**2/(2*error**2))
    prior=prior_func(p,fixp,guess,types)
    return(re+prior)
def prior_func(pval,fixp,guess,types):
    '''Prior probability function for 2D Gaussian model
    INPUT:  pval - parameter array of floats consisting of [amp, x, y, bmaj(fwhm), bmin(fwhm), bpa] (units mJy, pix, pix, arcsec, arcsec, deg)
            fixp - boolean array indicating which parameters are fixed in fit; same dimensions as pval
            guess - array of floats of initial guesses for parameters; same dimensions as pval
            types - 'toti' or 'pol' (str)
    OUTPUT: prior probability
    '''
    pv=[]
    for i in range(0,len(fixp)):
        if fixp[i]==True:
            pv.append(guess[i])
        elif fixp[i]==False:
            pv.append(pval[i])
        else:
            raise ValueError('The fixed param array values can only be True or False')
    p=np.array(pv)
    amp,xx,yy,bmaj,bmin,bpa=p[0],p[1],p[2],p[3],p[4],p[5]
    prior=0.0
    if types=='toti':
        prior += ss.norm.logpdf(amp,loc=guess[0],scale=100.)+ss.uniform.logpdf(amp,loc=0,scale=10e3)
    elif types=='pol':
        prior += ss.norm.logpdf(amp,loc=guess[0],scale=10.)+ss.uniform.logpdf(amp,loc=-1000,scale=2000)
    prior += ss.norm.logpdf(xx,loc=guess[1],scale=2.)+ss.uniform.logpdf(xx,loc=guess[1]-5,scale=guess[1]+10)
    prior += ss.norm.logpdf(yy,loc=guess[2],scale=2.)+ss.uniform.logpdf(yy,loc=guess[2]-5,scale=guess[2]+10)
    prior += ss.norm.logpdf(bmaj,loc=guess[3],scale=5.)+ss.uniform.logpdf(bmaj,loc=0,scale=30)
    prior += ss.norm.logpdf(bmin,loc=guess[4],scale=5.)+ss.uniform.logpdf(bmin,loc=0,scale=30)
    prior += ss.norm.logpdf(bpa,loc=guess[5],scale=10.)
    if np.isnan(prior):
        return(-np.inf)
    return(prior)
def mcmc_fit(data,error,guess,fixp,nBurn,nSample,flag,pixsize,outdir,timebin,types):
    '''MCMC fitter
    INPUT:  data - 2D fits image (array)
            error - 2D fits image variance (array)
            guess - array of floats of initial guesses for parameters; same dimensions as parameter array
            fixp - boolean array indicating which parameters are fixed in fit; same dimensions as parameter array
            nBurn - number of steps for burn in sampling (int)
            nSample - number of steps for sampling (int)
            flag - 'y' or 'n' to make diagnostic MCMC chains and histograms plots (str)
            pixsize - pixel size in the image in arcsec (float)
            outdir - output directory for diagnostic plots (str)
            timebin - size of lightcurve timebins in seconds (float)
            types - 'toti' or 'pol' (str)
    OUTPUT: bestfit parameter arrays [median, lower bound, upper bound]
            amp - best fit amplitude (mJy)
            xx - best fit X pixel position (pix)
            yy - best fit Y pixel position (pix)
            bmaj - best fit major axis (arcsec)
            bmin - best fit minor axis (arcsec)
            bpa - best fit postion angle (deg)
    '''
    ndim=6
    nwalkers=ndim*10
    p0=np.zeros((nwalkers,ndim))
    for i in np.arange(ndim):
        if fixp[i]==True:
            p0[:,i]=guess[i]
        elif fixp[i]==False:
            p0[:,i]=(((np.random.randn(nwalkers))*0.01)+guess[i])
    print('logprob:',lp(guess,data,error,fixp,guess,pixsize,types))
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lp,args=[data,error,fixp,guess,pixsize,types],pool=pool)
        with ProgressBar(nBurn) as bar:
            for i,(pos,prob,state) in enumerate(sampler.sample(p0,iterations=nBurn,skip_initial_state_check=True)):
                bar.update()
        sampler.reset()
        with ProgressBar(nSample) as bar:
            for i,(pos,prob,state) in enumerate(sampler.sample(pos,iterations=nSample,skip_initial_state_check=True)):
                bar.update()
            
    amp=confidenceInterval(sampler.flatchain[:,0])
    xx=confidenceInterval(sampler.flatchain[:,1])
    yy=confidenceInterval(sampler.flatchain[:,2])
    bmaj=confidenceInterval(sampler.flatchain[:,3])
    bmin=confidenceInterval(sampler.flatchain[:,4])
    bpa=confidenceInterval(sampler.flatchain[:,5])
    if flag=='y':
        plt.rcdefaults()
        fig=plt.figure(figsize=(10,3))
        for i in range(0, ndim):
            plt.subplot(1,6,i+1)
            patches=plt.hist(sampler.flatchain[:,i],bins=100)
        fig.subplots_adjust(hspace=0.5)
        plt.savefig(outdir+'hist_'+timebin+'.pdf')
        plt.close()
        #plt.show()
        fig=plt.figure(figsize=(10,3))
        for i in range(0, ndim):
            plt.subplot(1,6,i+1)
            plt.plot(sampler.chain[:,:,i].T)
        fig.subplots_adjust(hspace=0.5)
        plt.savefig(outdir+'trace_'+timebin+'.pdf')
        plt.close()
        #plt.show()
    return(amp,xx,yy,bmaj,bmin,bpa)
def fitting(outf,fitsim,num,cal_scan,flux_guess,pos_guess,w,rms,nburn,nsample,waveband,numscans,outdir,diag,types,plane,data_dir):
    '''Perform fitting process for a single map; LSQ used as initial guess to input in MCMC
    INPUT:  outf - output directory for fit results (str)
            fitsim - fits image of data to fit (str)
            num - target scan number (str)
            cal_scan - cal scan or 'default' if using auto calibrate option (str)
            flux_guess - intial guess for flux value (float)
            pos_guess - list of floats of initial guesses for [X,Y] pixel positions
            w - list of floats consisting of [bmaj, bmin, bpa]
            rms - error of the flux estimate for the [num] scan (float)
            nburn - number of steps for burn in sampling (int)
            nsample - number of steps for sampling (int)
            waveband - '8' or '4' indicating 850um or 450um (str)
            numscans - number of scans the fits image is made from (int; 1 is single scan, >1 is mosaic map)
            outdir - output directory for diagnostic plots (str)
            diag - 'y' or 'n' to turn on creation of mcmc diagnostic plots (str)
            types - 'toti', or 'pol' (str)
            plane - for types=='pol', I,Q,or U plane (0,1,2; int)
            data_dir - main output directory for observation date (str)
    OUPUT:  Nothing returned byt function will create,
            - txt file (fit_results.txt) with:
                column 1: MJD midpoint
                column 2: MJD error
                column 3: scan number
                column 4: Stokes
                column 5: flux in mJy/beam
                column 6: flux error in mJy/beam
                column 7: X position (RA in deg)
                column 8: Y position (Dec in deg)
            - Plot of the map with bestfit overlayed
    '''
    final_outputdir=data_dir+'data_products/'
    if not os.path.isdir(outdir+'scan_'+num+'/'):
        os.makedirs(outdir+'scan_'+num+'/')
    pixsize=Angle(abs(getheader(fitsim)['CDELT1'])*u.degree).arcsec
    header_target=getheader(fitsim)
    wmaptar=wcs.WCS(header_target)
    if os.path.isfile(outf):       
        outfile = open(outf,'a')
    elif not os.path.isfile(outf): 
        outfile = open(outf,'w')
        outfile.write('# MJD MJDerr scan stokes flux error X Y\n')
        outfile.write('# (dys) (dys) () () (mJy) (mJy) (deg) (deg)\n') 
    stokesvals={0:'I',1:'Q',2:'U'}
    #fit target by fixing beam size
    param_guess=[flux_guess,pos_guess[0],pos_guess[1],w[0],w[1],w[2]]
    if types == 'toti':
        data_target=getdata(fitsim,0)[0,:,:]
        vartarget=getdata(fitsim,1)
        tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
        bounds={'amplitude':(1,10000),'x_mean':(param_guess[1]-10,param_guess[1]+10),'y_mean':(param_guess[1]-10,param_guess[1]+10)}
        mcmc_fix=[False,False,False,True,True,True]
        fit_types='toti'
    elif types == 'pol':
        data_target=getdata(fitsim,0)[plane,0,:,:]
        vartarget=getdata(fitsim,1)[plane,0,:,:]
        if plane == 0:
            tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
            bounds={'amplitude':(1,10000),'x_mean':(param_guess[1]-10,param_guess[1]+10),'y_mean':(param_guess[1]-10,param_guess[1]+10)}
            mcmc_fix=[False,False,False,True,True,True]
            fit_types='toti'
        elif plane == 1 or plane == 2:
            tar_fixed_params={'amplitude':False,'x_mean':True,'y_mean':True,'x_stddev':True,'y_stddev':True,'theta':True}
            bounds={'amplitude':(-1*abs(flux_guess)*0.7,abs(flux_guess)*0.7)}
            mcmc_fix=[False,True,True,True,True,True]
            fit_types='pol'
    x = np.arange(0,len(data_target[0,:]))
    y = np.arange(0,len(data_target[:,0]))
    X,Y = np.meshgrid(x,y)
    Z=np.nan_to_num(data_target*1e3)
    M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3]/(2.*pixsize),\
        param_guess[4]/(2.*pixsize),param_guess[5],fixed=tar_fixed_params,bounds=bounds)
    lmff=astropy.modeling.fitting.LevMarLSQFitter()
    res=lmff(M,X,Y,Z)
    ampt=[res.parameters[0]]#,np.sqrt(np.diag(lmff.fit_info['param_cov']))[0]]
    amp_vall=ampt[0]
    x_vall=res.parameters[1]
    y_vall=res.parameters[2]
    param_guess=[amp_vall,x_vall,y_vall,w[0],w[1],w[2]]
    print('lsq:',param_guess)
    ampt,xxt0,yyt0,bmajt,bmint,bpat=mcmc_fit((data_target)*1e3,(vartarget)*1e3,param_guess,mcmc_fix,nburn,nsample,diag,\
                                                                        pixsize,outdir+'scan_'+num+'/',str(0),fit_types)
    mjderr=numscans*((Time(header_target['DATE-END'],format='isot',scale='utc').mjd-Time(header_target['DATE-OBS'],format='isot',scale='utc').mjd)/2.)
    mjdmid=Time(header_target['DATE-OBS'],format='isot',scale='utc').mjd+mjderr
    #print results
    if types=='toti':
        ra=float(wmaptar.wcs_pix2world(xxt0[0],yyt0[0],0,1)[0]) 
        dec=float(wmaptar.wcs_pix2world(xxt0[0],yyt0[0],0,1)[1])
    elif types=='pol':
        ra=float(wmaptar.wcs_pix2world(xxt0[0],yyt0[0],0,0,1)[0])
        dec=float(wmaptar.wcs_pix2world(xxt0[0],yyt0[0],0,0,1)[1])
    outfile.write('{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(mjdmid,mjderr,num,stokesvals[plane],ampt[0],rms,ra,dec))
    outfile.close()

    #show results
    mpl.rcParams['xtick.direction']='in'
    mpl.rcParams['ytick.direction']='in'
    x = np.arange(0,len(data_target[0,:]))
    y = np.arange(0,len(data_target[:,0]))
    X,Y = np.meshgrid(x,y)
    fig=plt.figure()
    ax=plt.subplot(111,projection=wmaptar.celestial)
    ax.imshow(Z,origin='lower')
    res2=models.Gaussian2D(abs(ampt[0]),xxt0[0],yyt0[0],param_guess[3]/(2.*pixsize),param_guess[4]/(2.*pixsize),param_guess[5])
    ax.contour(res2(X,Y),colors='w',levels=abs(ampt[0])*np.array([0.2,0.5,0.8]))
    ax.coords['ra'].set_axislabel('Right Ascension')
    ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
    ax.coords['ra'].set_major_formatter('hh:mm:ss')
    ax.set_title(Time(mjdmid,format='mjd').iso+': '+str('{0:.2f}'.format(ampt[0]))+'+/-'+str('{0:.2f}'.format(rms))+'mJy')
    plt.savefig(outdir+'scan_'+num+'/'+'target_fit'+waveband+'_'+num+'_'+str(integ)+'_'+stokesvals[plane]+'.pdf')
    plt.close()
    #plt.show()
    if cal_scan=='default':
        os.system('cp -r '+outdir+'scan_'+num+'/'+'target_fit'+waveband+'_'+num+'_'+str(integ)+'_'+stokesvals[plane]+'.pdf '+final_outputdir)


def fitting_timing(outf,fitsim,num,flux_guess,pos_guess,w,nburn,nsample,integ,waveband,outdir,diag):
    '''Perform fitting process for a cube; LSQ used as initial guess to input in MCMC
    INPUT:  outf - output directory for fit results (str)
            fitsim - fits image of data to fit (str)
            num - target scan number (str)
            flux_guess - intial guess for flux value (float)
            pos_guess - list of floats of initial guesses for [X,Y] pixel positions
            w - list of floats consiting of [bmaj, bmin, bpa]
            nburn - number of steps for burn in sampling (int)
            nsample - number of steps for sampling (int)
            integ - size of light curve timebins in seconds (float)
            waveband - '8' or '4' indicating 850um or 450um (str)
            outdir - output directory for diagnostic plots (str)
            diag - 'y' or 'n' to turn on creation of mcmc diagnostic plots (str)
    OUPUT:  Nothing returned but function will create,
            - txt file of lightcurve:
                column 1: MJD midpoint
                column 2: flux in mJy/beam
                column 3: flux error in mJy/beam
            - Plots of the map with bestfit overlayed for each time slice
    '''
    if not os.path.isdir(outdir+'scan_'+num+'/'):
        os.makedirs(outdir+'scan_'+num+'/')
    pixsize=Angle(abs(getheader(fitsim)['CDELT1'])*u.degree).arcsec
    header_target=getheader(fitsim)      
    data_target=getdata(fitsim,0)
    wmaptar=wcs.WCS(header_target)
    vartarget=getdata(fitsim,1)
    

    size=30
    naxis1 = fits.open(fitsim)[0].header['NAXIS1']
    naxis2 = fits.open(fitsim)[0].header['NAXIS2']
    center = np.array((int(naxis1/2.0),int(naxis2/2.0))) # gives x,y
    ranges = np.array((center[0]-size,center[0]+size,center[1]-size,center[1]+size)) # gives x1,x2,y1,y2

    mjd_start=Time(header_target['DATE-OBS'],format='isot',scale='utc').mjd
    MJDs=mjd_start+(np.arange(0,header_target['NAXIS3'])*integ)/(3600.*24.)

    param_guess=[flux_guess,pos_guess[0]*(naxis1/(400/pixsize))*(size*2/naxis1),pos_guess[1]*(naxis2/(400/pixsize))*(size*2/naxis2),w[0],w[1],w[2]]
    tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
    bounds={'amplitude':(1,10000),'x_mean':(param_guess[1]-10,param_guess[1]+10),'y_mean':(param_guess[1]-10,param_guess[1]+10)}
    mcmc_fix=[False,False,False,True,True,True]
    flux=[]
    err=[]
    for i in range(0,header_target['NAXIS3']):
        x = np.arange(0,len(data_target[i,0,ranges[2]:ranges[3]]))
        y = np.arange(0,len(data_target[i,ranges[0]:ranges[1],0]))
        X,Y = np.meshgrid(x,y)
        Z=np.nan_to_num(data_target[i,ranges[0]:ranges[1],ranges[2]:ranges[3]])*1e3
        M=models.Gaussian2D(param_guess[0],param_guess[1],param_guess[2],param_guess[3]/(2.*pixsize),\
            param_guess[4]/(2.*pixsize),param_guess[5],fixed=tar_fixed_params,bounds=bounds)
        lmff=astropy.modeling.fitting.LevMarLSQFitter()
        res=lmff(M,X,Y,Z)
        ampt=[res.parameters[0]]#,np.sqrt(np.diag(lmff.fit_info['param_cov']))[0]]
        amp_vall=ampt[0]
        x_vall=res.parameters[1]
        y_vall=res.parameters[2]
        param_guess=[amp_vall,x_vall,y_vall,w[0],w[1],w[2]]
        ampt,xxt0,yyt0,bmajt,bmint,bpat=mcmc_fit((data_target[i,ranges[0]:ranges[1],ranges[2]:ranges[3]])*1e3,\
            (vartarget[i,ranges[0]:ranges[1],ranges[2]:ranges[3]])*1e3,param_guess,mcmc_fix,nburn,nsample,diag,\
            pixsize,outdir+'scan_'+num+'/',str(i),'toti')
        amp_val=ampt[0]
        x_val=xxt0[0]
        y_val=yyt0[0]
        #show results
        mpl.rcParams['xtick.direction']='in'
        mpl.rcParams['ytick.direction']='in'
        fig=plt.figure()
        ax=plt.subplot(111,projection=wmaptar.celestial)
        ax.imshow(Z,origin='lower')
        res2=models.Gaussian2D(amp_val,x_val,y_val,param_guess[3]/(2.*pixsize),param_guess[4]/(2.*pixsize),param_guess[5])
        ax.contour(res2(X,Y),colors='w',levels=amp_val*np.array([0.2,0.5,0.8]))
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
        ax.coords['ra'].set_major_formatter('hh:mm:ss')
        rms=np.std(np.nan_to_num(data_target[i,ranges[0]:ranges[1],ranges[2]:ranges[3]]))*1e3
        ax.set_title(Time(MJDs[i],format='mjd').iso+': '+str('{0:.2f}'.format(amp_val))+'+/-'+str('{0:.2f}'.format(rms))+'mJy')
        plt.savefig(outdir+'scan_'+num+'/'+'target_fit'+waveband+'_'+str(i)+'_'+str(integ)+'.pdf')
        plt.close()
        #plt.show()
        flux.append(amp_val)
        err.append(rms)
    
   
    #print results
    outfile = open(outf,'w')
    outfile.write('#MJD flux error\n')
    outfile.write('#(dys) (mJy) (mJy)\n')
    for i in range(0,len(MJDs)):
        outfile.write('{0} {1} {2}\n'.format(MJDs[i],flux[i],err[i]))
    outfile.close()

def scuba2_analysis(data_dir,target_scan_numbers,cal_scan_numbers,waveband,integ,diag):
    ''' SCUBA-2 analysis procedure
    INPUT:  data_dir - output directory for data products (str)
            target_scan_numbers - list of strings of taget scan numbers
            cal_scan_numbers - list of strings of cal scan numbers of ['default'] if auto calibration done
            waveband - '8' or '4' indicating 850um or 450um (str) 
            integ - size of light curve timebins in seconds (float)
            diag - 'y' or 'n' to turn on creation of mcmc diagnostic plots (str)
    OUTPUT: Nothing returned by function creates,
            - file of full scan/mosaic fit results
            - light curves (data file and plot; if a timing cube is present)
            - plots of each map/time slice with best fits overlayed
            - mcmc diagnostic plots if turned on (diag)
    '''
    final_outputdir=data_dir+'/data_products/'
    for item in cal_scan_numbers:
        noise_file = data_dir+'calibrate_'+item+'/noise_log.txt'
        guesses = ascii.read(noise_file,names=['scan','type','noise','max_i','x','y','bmaj','bmin','bpa'])
        w = [float(guesses['bmaj'][0]),float(guesses['bmin'][0]),float(guesses['bpa'][0])]
        outf=data_dir+'calibrate_'+item+'/fit_results.txt'
        for targ in target_scan_numbers:
            targ_name = data_dir+'/calibrate_'+item+'/scan_'+targ+'_'+waveband+'_cal_crop.fits'
            rms = float(guesses['noise'][guesses['scan']==targ])
            flux_guess=float(guesses['max_i'][guesses['scan']==targ])
            pos_guess = [float(guesses['x'][guesses['scan']==targ]),float(guesses['y'][guesses['scan']==targ])]
            fitting(outf,targ_name,targ,item,flux_guess,pos_guess,w,rms,200,500,waveband,1,data_dir+'calibrate_'+item+'/diag_plots/',diag,'toti',0,data_dir)
            if integ!=0:
                outtiming=data_dir+'calibrate_'+item+'/target_lc'+waveband+'_'+targ+'_'+str(integ)+'_stokesI.txt'
                fits_cube=data_dir+'calibrate_'+item+'/scan_'+targ+'_'+waveband+'_shortmp_cube_cal.fits'
                full_fit=ascii.read(outf,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
                full_fit['scan']=full_fit['scan'].astype('str')
                pos_guess=cut_map(targ_name,[float(full_fit['X'][full_fit['scan']==targ]),float(full_fit['Y'][full_fit['scan']==targ])],'toti')
                flux_guess=float(full_fit['flux'][full_fit['scan']==targ])
                fitting_timing(outtiming,fits_cube,targ,flux_guess,pos_guess,w,200,500,integ,waveband,\
                    data_dir+'calibrate_'+item+'/diag_plots/',diag)
        if integ!=0:        
            make_lightcurves(item,waveband,integ,data_dir)
        mosaic_name = data_dir+'/calibrate_'+item+'/mosaic_map_scan_'+waveband+'.fits'
        rms = float(guesses['noise'][guesses['type']=='mosaic'])
        flux_guess=float(guesses['max_i'][guesses['type']=='mosaic'])
        pos_guess = [float(guesses['x'][guesses['type']=='mosaic']),float(guesses['y'][guesses['type']=='mosaic'])]
        fitting(outf,mosaic_name,'mosaic',item,flux_guess,pos_guess,w,rms,200,500,waveband,len(target_scan_numbers),data_dir+'calibrate_'+item+'/diag_plots/',diag,'toti',0,data_dir)

def make_lightcurves(item,waveband,integ,data_dir):
    '''Create lightcurve file of all the scans combined and plots full lightcurve
    INPUT:  item - cal scan or 'default' if using auto calibrate option (str)
            waveband - '8' or '4' indicating 850um or 450um (str) 
            integ - size of light curve timebins in seconds (float)
            data_dir - output directory for lightcurves (str)
    OUTPUT: Nothing returned but function will create,
            - combined scan light curve file with:
                column 1: MJD midpoint
                column 2: flux in mJy/beam
                column 3: flux error in mJy/beam
            - lightcurve plot for all scans combined 
    '''
    final_outputdir=data_dir+'data_products/'
    all_scan_lcs=glob.glob(data_dir+'calibrate_'+item+'/target_lc'+waveband+'_*_'+str(integ)+'_stokesI.txt')
    all_scan_lcs.sort()
    
    with open(data_dir+'calibrate_'+item+'/target_lc'+waveband+'_all_'+str(integ)+'_stokesI.txt','w') as writer:
        readers=[open(filename) for filename in all_scan_lcs]
        print('#MJD flux error',file=writer)
        print('#(dys) (mJy) (mJy)',file=writer)
        for red in readers:
            for lines in red:
                if not '#' in lines:
                    print((lines.strip()),file=writer)
    lc=ascii.read(data_dir+'calibrate_'+item+'/target_lc'+waveband+'_all_'+str(integ)+'_stokesI.txt',names=('MJD','Flux','Error'))
    fig=plt.figure()
    font={'family':'serif','size':14}
    rc('font',**font)
    mpl.rcParams['xtick.direction']='in'
    mpl.rcParams['ytick.direction']='in'
    ax=plt.subplot(111)
    ax.errorbar(lc['MJD'],lc['Flux'],yerr=lc['Error'],marker='o',ls='',color='m')
    ax.set_ylabel('Flux Density (mJy)',fontsize=15)
    ax.set_xlabel('Time on '+Time(lc['MJD'][0],format='mjd').iso.split(' ')[0]+' (HH:MM)',fontsize=15)
    locator=mdates.MinuteLocator(interval=30)
    locator2=mdates.MinuteLocator(interval=5)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_minor_locator(locator2)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    plt.setp(ax.get_xticklabels(),rotation=45,horizontalalignment='right')
    #ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis='x',which='major',labelsize=15,length=7,width=1.5,top='on',bottom='on',pad=7)
    ax.tick_params(axis='x',which='minor',labelsize=15,length=5,width=1.,top='on',bottom='on',pad=7)
    ax.tick_params(axis='y',which='major',labelsize=15,length=7,width=1.5,right='on',left='on',pad=7)
    ax.tick_params(axis='y',which='minor',labelsize=15,length=5,width=1,right='on',left='on',pad=7)
    plt.savefig(data_dir+'calibrate_'+item+'/'+'target_fit'+waveband+'_all_'+str(integ)+'_stokesI.pdf')
    plt.close()
    if item=='default':
        os.system('cp -r '+data_dir+'calibrate_'+item+'/target_lc'+waveband+'_all_'+str(integ)+'_stokesI.txt '+final_outputdir)
        os.system('cp -r '+data_dir+'calibrate_'+item+'/'+'target_fit'+waveband+'_all_'+str(integ)+'_stokesI.pdf '+final_outputdir)

def pol2_analysis(data_dir,target_scan_numbers,waveband,integ,diag):
    ''' POL-2 analysis procedure
    INPUT:  data_dir - output directory for data products (str)
            target_scan_numbers - list of strings of taget scan numbers
            waveband - '8' or '4' indicating 850um or 450um (str) 
            integ - size of light curve timebins in seconds (float)
            diag - 'y' or 'n' to turn on creation of mcmc diagnostic plots (str)
    OUTPUT: Nothing returned by function creates,
            - file of full IQU scan/mosaic flux fit results
            - file of pol properties (PA,LP%)
            - light curves for I (data file and plot; if a timing cube is present)
            - plots of each IQU map/I time slice with best fits overlayed
            - mcmc diagnostic plots if turned on (diag)
    '''
    noise_file = data_dir+'calibrate_default/noise_log.txt'
    guesses = ascii.read(noise_file,names=['scan','type','noise','max_i','x','y','bmaj','bmin','bpa'])
    w = [float(guesses['bmaj'][0]),float(guesses['bmin'][0]),float(guesses['bpa'][0])]
    outf=data_dir+'calibrate_default/fit_results.txt'
    for targ in target_scan_numbers:
        targ_name = data_dir+'calibrate_default/'+'stokes_i_'+waveband+'/'+date+'_000'+targ+'_stokes_cube_full.fits'
        rms = float(guesses['noise'][np.logical_and(guesses['scan']==targ,guesses['type']=='I')])
        flux_guess=float(guesses['max_i'][np.logical_and(guesses['scan']==targ,guesses['type']=='I')])
        pos_guess = [float(guesses['x'][np.logical_and(guesses['scan']==targ,guesses['type']=='I')]),\
        float(guesses['y'][np.logical_and(guesses['scan']==targ,guesses['type']=='I')])]
        fitting(outf,targ_name,targ,'default',flux_guess,pos_guess,w,rms,200,500,waveband,1,data_dir+'calibrate_default/diag_plots/',diag,'pol',0,data_dir)
        guesses2 = ascii.read(outf,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
        guesses2['scan']=guesses2 ['scan'].astype('str')
        flux_guess=float(guesses2['flux'][np.logical_and(guesses2['scan']==targ,guesses2['stokes']=='I')])
        pos_guess=cut_map(targ_name,[float(guesses2['X'][np.logical_and(guesses2['scan']==targ,guesses2['stokes']=='I')]),\
            float(guesses2['Y'][np.logical_and(guesses2['scan']==targ,guesses2['stokes']=='I')])],'pol')
        rms = float(guesses['noise'][np.logical_and(guesses['scan']==targ,guesses['type']=='Q')])
        fitting(outf,targ_name,targ,'default',flux_guess*0.1,pos_guess,w,rms,200,500,waveband,1,data_dir+'calibrate_default/diag_plots/',diag,'pol',1,data_dir)
        rms = float(guesses['noise'][np.logical_and(guesses['scan']==targ,guesses['type']=='U')])
        fitting(outf,targ_name,targ,'default',flux_guess*0.1,pos_guess,w,rms,200,500,waveband,1,data_dir+'calibrate_default/diag_plots/',diag,'pol',2,data_dir)
        pol2_calcu(outf,targ)
        if integ!=0:
            outtiming=data_dir+'calibrate_default/target_lc'+waveband+'_'+targ+'_'+str(integ)+'_stokesI.txt'
            fits_cube=glob.glob(data_dir+'calibrate_default/stokes_icube_'+waveband+'/'+'scan_'+targ+'_Imap_cube_cal.fits')[0]
            fitting_timing(outtiming,fits_cube,targ,flux_guess,pos_guess,w,200,500,integ,waveband,\
                data_dir+'calibrate_default/diag_plots/',diag)
    if integ!=0:
        make_lightcurves('default',waveband,integ,data_dir)
    mosaic_name = data_dir+'/calibrate_default/stokes_i_'+waveband+'/mosaic_map_stokes_cube_full.fits'
    rms = float(guesses['noise'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='I')])
    flux_guess=float(guesses['max_i'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='I')])
    pos_guess = [float(guesses['x'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='I')]),\
    float(guesses['y'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='I')])]
    fitting(outf,mosaic_name,'mosaic','default',flux_guess,pos_guess,w,rms,200,500,waveband,len(target_scan_numbers),data_dir+'calibrate_default/diag_plots/',diag,'pol',0,data_dir)
    guesses2 = ascii.read(outf,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
    flux_guess=float(guesses2['flux'][np.logical_and(guesses2['scan']=='mosaic',guesses2['stokes']=='I')])
    pos_guess=cut_map(targ_name,[float(guesses2['X'][np.logical_and(guesses2['scan']==targ,guesses2['stokes']=='I')]),\
            float(guesses2['Y'][np.logical_and(guesses2['scan']==targ,guesses2['stokes']=='I')])],'pol')
    rms = float(guesses['noise'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='Q')])
    fitting(outf,mosaic_name,'mosaic','default',flux_guess*0.1,pos_guess,w,rms,200,500,waveband,len(target_scan_numbers),data_dir+'calibrate_default/diag_plots/',diag,'pol',1,data_dir)
    rms = float(guesses['noise'][np.logical_and(guesses['scan']=='mosaic',guesses['type']=='U')])
    fitting(outf,mosaic_name,'mosaic','default',flux_guess*0.1,pos_guess,w,rms,200,500,waveband,len(target_scan_numbers),data_dir+'calibrate_default/diag_plots/',diag,'pol',2,data_dir)
    pol2_calcu(outf,'mosaic')
def pol2_calcu(outf,scant):
    '''Calcualtes polarization PA and LP % from stokes fitting results
    INPUT:  outf - flux fitting log file (str)
            scant - target scan number or 'mosaic' (str)
    OUTPUT: Function returns nothing, but creates:
            - log file of polarization properties
            column1: scan
            column2: I (mJy)
            column3: Q (mJy)
            column4: U (mJy)
            column5: Ierr (mJy)
            column6: Qerr (mJy)
            column7: Uerr (mJy)
            column8: pol detection (True or False)
            column9: LP %
            column10: LP% err
            column11: PA (deg)
            column12: PAerr (deg)
    '''
    full_fit=ascii.read(outf,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
    full_fit['scan']=full_fit['scan'].astype('str')
    I=np.array(full_fit['flux'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='I')])[0]
    Q=np.array(full_fit['flux'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='Q')])[0]
    U=np.array(full_fit['flux'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='U')])[0]
    Ierr=np.array(full_fit['error'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='I')])[0]
    Qerr=np.array(full_fit['error'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='Q')])[0]
    Uerr=np.array(full_fit['error'][np.logical_and(full_fit['scan']==scant,full_fit['stokes']=='U')])[0]
    if Q>2*Qerr and U>2*Uerr:
        det='True'
    else:
        det='False'
    LP=100*np.sqrt(Q**2+U**2)/I
    LPerr=(100/I)*np.sqrt((Q*Qerr/(np.sqrt(Q**2+U**2)))**2+(U*Uerr/(np.sqrt(U**2+Q**2)))**2)
    PA=0.5*np.arctan(U/Q)*(180/np.pi)
    PAerr=0.5*(180/np.pi)*np.sqrt((Uerr*Q/(1+U**2))**2+(Qerr*U/(1+Q**2))**2)
    if os.path.isfile(outf.strip('.txt')+'_pol.txt'):       
        outfile = open(outf.strip('.txt')+'_pol.txt','a')
    elif not os.path.isfile(outf.strip('.txt')+'_pol.txt'): 
        outfile = open(outf.strip('.txt')+'_pol.txt','w')
        outfile.write('#scan I Q U Ierr Qerr Uerr detection LP LPerr PA PAerr\n')
        outfile.write('#() (mJy) (mJy) (mJy) (mJy) (mJy) (mJy) () (%) (%) (deg) (deg)\n') 
    outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(scant,I,Q,U,Ierr,Qerr,Uerr,det,LP,LPerr,PA,PAerr))
    outfile.close()

def create_total_log(direc,wout,obs_type):
    '''Write a summary log, that will be updated with each days observation.
    INPUT:  direc - path to data products directory for specific date (str)
            wout - location to store log in (str)
            obs_type - 'scuba2' or 'pol2' (str)
    OUTPUT: Function returns nothing, but creates,
            - A file with:
            column1: date of observation (iso format)
            column2: MJD
            column3: I flux (mJy)
            column4: Q flux (mJy)
            column5: U flux (mJy)
            column6: error in I flux (mJy)
            column7: error in Q flux (mJy)
            column8: error in U flux (mJy)
            column9: flag indicating if polarized signal detected (True or False)
            column10: linear polarization (%)
            column11: error in linear polarization (%)
            column12: polarization PA (deg)
            column13: error in polarization PA (deg)
    '''
    logfile=direc+'data_products/fit_results.txt'
    logfile2=direc+'data_products/fit_results_pol.txt'
    if not os.path.isfile(wout+'observation_log.txt'):
        fileo=open(wout+'observation_log.txt','w')
        fileo.write('#Date MJD I Q U Ierr Qerr Uerr poldet LP Lperr PA PAerr\n')
    else:
        fileo=open(direc+'observation_log.txt','a')
    if obs_type=='scuba2':
        data=ascii.read(logfile,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
        I=np.array(data['flux'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        Ierr=np.array(data['error'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        MJD=np.array(data['MJD'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        date=Time(float(MJD),format='mjd').iso
        fileo.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}\n'.format(date,MJD,I,0,0,Ierr,0,0,'False',0,0,0,0))
    elif obs_type=='pol2':
        data=ascii.read(logfile,names=('MJD','MJDerr','scan','stokes','flux','error','X','Y'))
        data2=ascii.read(logfile2,names=('scan','I','Q','U','Ierr','Qerr','Uerr','detection','LP','LPerr','PA','PAerr'))
        I=np.array(data['flux'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        Ierr=np.array(data['error'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        Q=np.array(data['flux'][np.logical_and(data['scan']=='mosaic',data['stokes']=='Q')])[0]
        Qerr=np.array(data['error'][np.logical_and(data['scan']=='mosaic',data['stokes']=='Q')])[0]
        U=np.array(data['flux'][np.logical_and(data['scan']=='mosaic',data['stokes']=='U')])[0]
        Uerr=np.array(data['error'][np.logical_and(data['scan']=='mosaic',data['stokes']=='U')])[0]
        MJD=np.array(data['MJD'][np.logical_and(data['scan']=='mosaic',data['stokes']=='I')])[0]
        date=Time(float(MJD),format='mjd').iso
        det=np.array(data2['detection'][data2['scan']=='mosaic'])[0]
        LP=np.array(data2['LP'][data2['scan']=='mosaic'])[0]
        LPerr=np.array(data2['LPerr'][data2['scan']=='mosaic'])[0]
        PA=np.array(data2['PA'][data2['scan']=='mosaic'])[0]
        PAerr=np.array(data2['PAerr'][data2['scan']=='mosaic'])[0]
        fileo.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}\n'.format(date,MJD,I,Q,U,Ierr,Qerr,Uerr,det,LP,LPerr,PA,PAerr))

####################################
#User input
####################################
date = '20150622'
diag='y'#print diagnostic plots from fitting (y or n)
wout='/export/data2/atetarenko/PB_test/v404/'
####################################

#define input/output directories
data_dir=wout+'raw/'+date+'/'
output_dir=wout+'results/'+date+'/'

#read in shared parameters used in data reduction script
fp = open('/export/data2/atetarenko/PB_test/v404/results/'+date+'/shared.pkl','rb')
shared = pickle.load(fp)
target = shared['target']
cal_type= shared['cal_type']
obs_type=shared['obs_type']
waveband=shared['waveband']
integ=shared['integ']

#create lists of target and calibrator scan numbers based on raw data directory
target_scan_numbers,cal_scan_numbers = get_scan_numbers(data_dir,date,target,cal_type,waveband)

if obs_type=='scuba2':
    scuba2_analysis(output_dir,target_scan_numbers,cal_scan_numbers,waveband,integ,diag)
    #move the flux log (from auto calibration) into the final processing directory
    os.system('cp -r '+output_dir+'calibrate_default/fit_results.txt '+output_dir+'/data_products')
    #add to total log of all days
    create_total_log(output_dir,wout+'results/','scuba2')
elif obs_type=='pol2':
    pol2_analysis(output_dir,target_scan_numbers,waveband,integ,diag)
    #move the flux/pol logs from auto calibration into the final processing directory
    os.system('cp -r '+output_dir+'calibrate_default/fit_results.txt '+output_dir+'/data_products')
    os.system('cp -r '+output_dir+'calibrate_default/fit_results_pol.txt '+output_dir+'/data_products')
    #add to total log of all days
    create_total_log(output_dir,wout+'results/','pol2')





