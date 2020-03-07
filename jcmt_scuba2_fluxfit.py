#######################################
# JCMT SCUBA-2 Flux Fitting Code
#######################################
'''Fits a 2D Gaussian to JCMT fits images using MCMC or a simple least squares algorithm.
INPUT: my_dir: Output directory
       fitsim: Target FITS image
       w: Beam Parameters from STARLINK; [BMAJ (arcsec), BMIN (arcsec), BPA (deg)]
       cal_im: (optional) Calibrator FITS image for directly fitting beam size
       ranges/ranges_cal: pixel ranges to search for source/calibrator in FITS images
OUTPUT: Results file (my_dir/fit_results.txt)
	Resulting images/models (my_dir/calibrator_fit.png,my_dir/target_fit.png)
NOTES: It is recommended to use image RMS as an uncertainty measurment on flux.
The beam is fixed in target fitting, so only flux and position vary.

Written by: Alex J. Tetarenko
Last Updated: Apr 22, 2019'''

#needed packages
from astropy.io import fits,ascii 
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.coordinates import Angle
from astropy import wcs
from astropy import units as u
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
import matplotlib as mpl
from astropy.time import Time
import emcee
from astropy.utils.console import ProgressBar
import scipy.stats as ss
from data_reduction import get_scan_numbers
def cut_map(num,cal_scan,output_dir,size):
    filename = output_dir+'calibrate_'+cal_scan+'/'+'scan_'+num+'_cal_crop.fits'
    naxis1 = fits.open(filename)[0].header['NAXIS1']
    naxis2 = fits.open(filename)[0].header['NAXIS2']
    center = np.array((int(naxis1/2.0),int(naxis2/2.0))) # gives x,y
    ranges = np.array((center[0]-size,center[0]+size,center[1]-size,center[1]+size)) # gives x1,x2,y1,y2
    return ranges

def confidenceInterval(y,sig):
	median=np.median(y)
	pct15=np.percentile(y,15)
	pct85=np.percentile(y,85)
	list1=np.array([median,median-pct15,pct85-median])
	return list1
def lp(p,data,error,fixp,guess,pixsize):
	#amp in mJy, xx and yy in pixels, bmaj and bmin in arcsec, bpa in deg
	amp,xx,yy,bmaj,bmin,bpa=p[0],p[1],p[2],p[3],p[4],p[5]
	mod0=models.Gaussian2D(amp,xx,yy,bmaj/(2.*pixsize),bmin/(2.*pixsize),bpa)
	xval=np.arange(0,len(data[0,:]))
	yval=np.arange(0,len(data[:,0]))
	Xval, Yval = np.meshgrid(xval, yval)
	mod1=mod0(Xval,Yval)
	re=-0.5*np.nansum(np.log(2*np.pi*error**2))-np.nansum((mod1-data)**2/(2*error**2))
	prior=prior_func(p,fixp,guess)
	return(re+prior)
def prior_func(pval,fixp,guess):
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
	prior += ss.norm.logpdf(amp,loc=guess[0],scale=10.)+ss.uniform.logpdf(amp,loc=0,scale=10e3)
	prior += ss.norm.logpdf(xx,loc=guess[1],scale=5.)+ss.uniform.logpdf(xx,loc=0,scale=2.*guess[1])
	prior += ss.norm.logpdf(yy,loc=guess[2],scale=5.)+ss.uniform.logpdf(yy,loc=0,scale=2.*guess[2])
	prior += ss.norm.logpdf(bmaj,loc=guess[3],scale=5.)+ss.uniform.logpdf(bmaj,loc=0,scale=30)
	prior += ss.norm.logpdf(bmin,loc=guess[4],scale=5.)+ss.uniform.logpdf(bmin,loc=0,scale=30)
	prior += ss.norm.logpdf(bpa,loc=guess[5],scale=10.)
	if np.isnan(prior):
		return(-np.inf)
	return(prior)
def mcmc_fit(data,error,guess,fixp,nBurn,nSample,flag,pixsize):
    ndim=6
    nwalkers=ndim*2
    p0=np.zeros((nwalkers,ndim))
    for i in np.arange(ndim):
        if fixp[i]==True:
            p0[:,i]=guess[i]
        elif fixp[i]==False:
            p0[:,i]=(((np.random.randn(nwalkers))*0.01)+guess[i])
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lp,args=[data,error,fixp,guess,pixsize])
    print('Burn in')
    with ProgressBar(nBurn) as bar:
        for i,(pos,prob,state) in enumerate(sampler.sample(p0,iterations=nBurn,skip_initial_state_check=True)):
            bar.update()
    sampler.reset()
    print('Sample')
    with ProgressBar(nSample) as bar:
        for i,(pos,prob,state) in enumerate(sampler.sample(pos,iterations=nSample,skip_initial_state_check=True)):
            bar.update()
    amp=confidenceInterval(sampler.flatchain[:,0],1)
    xx=confidenceInterval(sampler.flatchain[:,1],1)
    yy=confidenceInterval(sampler.flatchain[:,2],1)
    bmaj=confidenceInterval(sampler.flatchain[:,3],1)
    bmin=confidenceInterval(sampler.flatchain[:,4],1)
    bpa=confidenceInterval(sampler.flatchain[:,5],1)
    if flag=='y':
        plt.rcdefaults()
        fig=plt.figure(figsize=(10,3))
        for i in range(0, ndim):
            plt.subplot(1,6,i+1)
            patches=plt.hist(sampler.flatchain[:,i],bins=100)
        fig.subplots_adjust(hspace=0.5)
        plt.show()
        input('Press enter to continue')
        fig=plt.figure(figsize=(10,3))
        for i in range(0, ndim):
            plt.subplot(1,6,i+1)
            plt.plot(sampler.chain[:,:,i].T)
        fig.subplots_adjust(hspace=0.5)
        plt.show()
    return(amp,xx,yy,bmaj,bmin,bpa)
def fitting(outf,cal_im,fitsim,num,flux_guess,pos_guess,w,rms):
    '''
    -----
    Parameters: - outf: path to save the .txt with the fit results
                - cal_im: path to calibration scan  file
                - fitsim: path to target scan file
                - num: target scan number
                - flux_guess: list with flux values for [cal,target/mosaic], 
                  obtained from noise_log.txt
                - pos_guess: list with position [x,y] where the flux_guess is,
                  obtained from noise_log.txt
                - w: list with [bmaj, bmin, bpa], obtained from noise_log.txt
                - rms: integer with the error of the flux estimate
    Output: - A .txt file (fit_results.txt) with
                column 1: MJD midpoint
                column 2: MJD error
                column 3: scan number
                column 4: flux in mJy/beam
                column 5: flux error in mJy/beam
            - Plots of the maps with fit on top, for each scan
    '''
    print('Starting fitting process...''')
    print('Reading in fits headers...')
    pixsize=Angle(abs(getheader(cal_im)['CDELT1'])*u.degree).arcsec
    header_cal=getheader(cal_im)
    header_target=getheader(fitsim)
    wmapcal=wcs.WCS(header_cal)
    wmaptar=wcs.WCS(header_target)
    if os.path.isfile(outf):       
        outfile = open(outf,'a')
    elif not os.path.isfile(outf): 
        outfile = open(outf,'w') 

    #fit target by fixing beam size
    data_target=getdata(fitsim,0)[0,:,:]
    vartarget=getdata(fitsim,1)#[0,ranges[2]:ranges[3],ranges[0]:ranges[1]]
    print("Fitting Target "+num+"...")
    param_guess=[flux_guess[1],pos_guess[0],pos_guess[1],w[0],w[1],w[2]]
    tar_fixed_params={'amplitude':False,'x_mean':False,'y_mean':False,'x_stddev':True,'y_stddev':True,'theta':True}
    ampt,xxt0,yyt0,bmajt,bmint,bpat=mcmc_fit(data_target*1e3,vartarget*1e3,[param_guess[0],param_guess[1],param_guess[2],w[0],\
                                                                        w[1],w[2]],[False,False,False,True,True,True],500,1500,'n',pixsize)
    mjderr=((Time(header_target['DATE-END'],format='isot',scale='utc').mjd-Time(header_target['DATE-OBS'],format='isot',scale='utc').mjd)/2.)
    mjdmid=Time(header_target['DATE-OBS'],format='isot',scale='utc').mjd+mjderr
   
    #print results
    outfile.write('{0} {1} {2} {3} {4}\n'.format(mjdmid,mjderr,num,ampt[0],rms))
    outfile.close()

    #show results
    print('Plotting result...')
    x = np.arange(0,len(data_target[0,:]))
    y = np.arange(0,len(data_target[:,0]))
    X,Y = np.meshgrid(x,y)
    fig=plt.figure()
    ax=plt.subplot(111,projection=wmaptar.celestial)
    ax.imshow(data_target,origin='lower')
    res2=models.Gaussian2D(ampt[0],xxt0[0],yyt0[0],param_guess[3]/(2.*pixsize),param_guess[4]/(2.*pixsize),param_guess[5])
    ax.contour(res2(X,Y),colors='w',levels=ampt[0]*np.array([0.2,0.5,0.8]))
    ax.coords['ra'].set_axislabel('Right Ascension')
    ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
    ax.coords['ra'].set_major_formatter('hh:mm:ss')
    ax.set_title('Target Fit '+num+' (beam size/angle fixed)')
    plt.savefig(my_dir+'target_fit'+num+'.png')
    plt.show()	


print('Reading in parameters...')
####################################
#User input
####################################
#input/output directory
my_dir='/Users/constanza/Desktop/Research/Data_products/'
data_dir ='/Users/constanza/Desktop/Research/Data/'
date = '20150622'
target = 'BHXRB V404 Cyg'
cal_name_list = ['Arp220','CRL2688']
#fits images
fitsim,cal_im = get_scan_numbers(data_dir,date,target,cal_name_list)
cal_im = ['08']
#flux guesses for cal and target in mJy
noise_file = my_dir+'calibrate_08/noise_log.txt'
guesses = ascii.read(noise_file,names=['scan','type','noise','max_i','x','y','bmaj','bmin','bpa'])



#size of gaussian beam in arcsec,pa in deg
w = [float(guesses['bmaj'][0]),float(guesses['bmin'][0]),float(guesses['bpa'][0])]

outf=my_dir+'fit_results.txt'


#########################################

for item in cal_im:
    cal_name = my_dir+'/calibrate_'+item+'/scan_'+item+'_cal_crop.fits'
    mosaic_name = my_dir+'/calibrate_'+item+'/mosaic_map_crop.fits'
    rms = float(guesses['noise'][guesses['type']=='mosaic'])
    flux_guess=[float(guesses['max_i'][guesses['type']=='cal']),float(guesses['max_i'][guesses['type']=='mosaic'])]
    pos_guess = [float(guesses['x'][guesses['type']=='mosaic']),float(guesses['y'][guesses['type']=='mosaic'])]
    fitting(outf,cal_name,mosaic_name,'mosaic',flux_guess,pos_guess,w,rms)
    for targ in fitsim:
        targ_name = my_dir+'/calibrate_'+item+'/scan_'+targ+'_cal_crop.fits'
        # save noise in the map to a parameter (represents uncertainty in flux)
        rms = float(guesses['noise'][guesses['scan']==int(targ)])
        flux_guess=[float(guesses['max_i'][guesses['type']=='cal']),float(guesses['max_i'][guesses['scan']==int(targ)])]
        # x guess, y guess
        pos_guess = [float(guesses['x'][guesses['scan']==int(targ)]),float(guesses['y'][guesses['scan']==int(targ)])]
        fitting(outf,cal_name,targ_name,targ,flux_guess,pos_guess,w,rms)

#########################################


