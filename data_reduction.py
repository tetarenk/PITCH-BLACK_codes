##################################################
#JCMT SCUBA-2/POL-2 Data Reduction Script
##################################################
'''Python script to reduce JCMT SCUBA-2/POL-2 data, using the Starlink software package.

INPUT:  data_dir - raw data directory
        output_dir - output directory
        date - date of observations (yyymmdd)
        target - target name
        cal_type - 'user' or 'auto', indicating user calculated FCFs from calibrator scans (SCUBA-2 only) or use of standard FCF values
        obs_type - 'scuba2' or 'pol2'
        waveband - '8' or '4' for 850um or 450um
        integ - lightcurve timebins size, in units of seconds (0 if no lightcurve analysis to be done)
        crop_params - radius of circular mask for cropping maps, in units of arcsec

OUTPUT: SCUBA-2 Calibrated, cropped maps for each scan and mosaic of all scans, with and without a matched filter applied (sdf and fits files)
                Calibrated timing cubes for each scan and mosaic of all scans (sdf and fits files)
        POL-2   Calibrated, cropped stokes cubes for each scan and mosaic of scans, with and without a matched filter applied (sdf and fits files)
                Calibrated timing cubes (I,Q, and U) for each scan (sdf and fits files)
        rms noise estimation log (only includes no matched filter results)

NOTES:  - Requires an installed version of the Starlink Software (http://starlink.eao.hawaii.edu/starlink) before running.
        - Script makes use of the starlink-pywrapper python package (https://github.com/Starlink/starlink-pywrapper).
        - User must set the path to their Starlink installation before running this script (line 39).

Originally written by Coni Echiburu Trujillo
Updates by Alex J. Tetarenko

Last Updated: April 17, 2020
'''

#packages to import
from starlink import wrapper
from starlink import kappa
from starlink import smurf
from starlink import picard
from starlink import convert
from starlink import polpack
wrapper.change_starpath('/star')
import glob
import numpy as np
import os
from astropy.io import ascii,fits
from astropy import wcs,coordinates
from astropy import units as u
from astropy.coordinates import Angle
from astropy.io.fits import getheader,getdata
import pickle

def create_crop_config(output_dir,radius_arcsec):
    '''Create a parameter file for use by PICARD cropping recipe.
    INPUT:  output_dir - output directory to store file in (str)
            radius_arcsec - radius of circular region in arcsec to use in cropping procedure (float)
    OUTPUT: name of cropping parameter file (str)
    '''
    filecrop=open(output_dir+'crop_parameters.lis','w')
    filecrop.write('[CROP_SCUBA2_IMAGES]\n')
    filecrop.write('CROP_METHOD=CIRCLE\n')
    filecrop.write('MAP_RADIUS='+str(radius_arcsec)+'\n')
    return(output_dir+'crop_parameters.lis')

def create_configs(starpath,output_dir,obs_type,bins):
    '''Create config files for use by the mapmaker or pol2map tasks.
    INPUT:  starpath - path to starlink installation
            output_dir - output directory to store files in (str)
            obs_type - type of observation; 'scuba2' or 'pol2' (str)
            bins - lightcurve bin size in seconds (float; 0 means no timing config file made)
    OUTPUT: file_config - name of standard config file (str)
            file_config_timing - name of timing config file (if bins!=0; str)
    '''
    if obs_type=='scuba2':
        cfile=starpath+'/share/smurf/dimmconfig_bright_compact.lis'
        file_config=output_dir+'dimmconfig_bright_compact.lis'
        file_config_timing=output_dir+'dimmconfig_bright_compact_shortmaps.lis'
    elif obs_type=='pol2':
        cfile=starpath+'/share/smurf/dimmconfig_pol2_compact.lis'
        file_config=output_dir+'dimmconfig_pol2_compact.lis'
        file_config_timing=output_dir+'dimmconfig_pol2_compact_shortmaps.lis'
    os.system('rm -rf '+file_config)
    os.system('cp -r '+cfile+' '+file_config)
    if bins!=0:
        os.system('rm -rf '+file_config_timing)
        os.system('cp -r '+cfile+' '+file_config_timing)
        fileconfig=open(file_config_timing,'a')
        fileconfig.write('#adding selected timebin size to make maps\n')
        fileconfig.write('shortmap=-'+str(bins))#note the negative value for shortmap means value is in seconds
        fileconfig.close()
    else:
        file_config_timing=''
    return(file_config,file_config_timing)

def get_scan_numbers(data_dir,date,target,cal_type,waveband):
    '''Grabs the scan numbers using the headers of all raw data files in a directory by
    looking for OBJECT and OBS_TYPE keywords.
    INPUT:  data_dir - path to directory where the raw data is stored (str)
            date - date of observation (YYYYMMDD; str)
            target - name of the target (str)
            cal_type - 'user' (outputs list of all cal scans present for manual FCF determination + 'default') or
            'auto' (outputs ['default'] for application of standard FCFs; str)
            waveband - '8' or '4' for 850um or 450um (str)
    OUTPUT: list of strings with target scan numbers
            list of strings with calibration scan numbers or ['default'] if auto cal_type selected
    '''
    names = glob.glob(data_dir+'s'+waveband+'*.sdf')
    scan_numbers = [i.split('_')[2].lstrip('000') for i in names]
    unique_scan_numbers = np.unique(np.array(scan_numbers))
    unique_scan_numbers_final = [str(item).zfill(2) for item in unique_scan_numbers]
    final_target_scan_numbers = []
    final_cal_scan_numbers = []
    for num in unique_scan_numbers_final:
        obj=kappa.fitsval(data_dir+'s'+waveband+'a'+date+'_000'+num+'_0001.sdf','OBJECT').value
        mode=kappa.fitsval(data_dir+'s'+waveband+'a'+date+'_000'+num+'_0001.sdf','OBS_TYPE').value
        if obj == target and mode == 'science':
            final_target_scan_numbers.append(num)
        if cal_type=='user':
            if obj != target and mode == 'science':
                final_cal_scan_numbers.append(num)
    final_cal_scan_numbers.append('default')
    return final_target_scan_numbers, final_cal_scan_numbers

def check_data():
    '''Checks if there are corrupt files before reducing all the data
    ------------------------------------------------------------------------------
    '''   
    #how do i do this????
            
def FCF(cal_name):
    '''Calculates the FCF value for a given calibrator observation.
    INPUT:  cal_name - name of the calibrator map (str)
    OUTPUT: FCF value (float)
    '''
    result_FCF = picard.scuba2_check_cal([cal_name])
    fcffile=[s for s in result_FCF.logfiles if 'log.fcf' in s][0]
    FCF = float(ascii.read(fcffile,names=('UT','HST','Obs','Source','Mode','Filter','El','Tau225',\
    'Tau','FCF_ARCSEC','FCFerr1','FCF_BEAM','FCFerr2','FCF_BEAMMATCH','FCFerr3','FWHM_eff'))['FCF_BEAM']) # gets the FCF value
    return FCF

def create_list(num,data_dir,output_dir,cal_scan,waveband):
    '''Creates file listing raw data files for a particular scan; for use by mapmaker or pol2map.
    INPUT:  num - target scan number (str or list of strings)
            data_dir - path to directory where the raw data is stored (str)
            output_dir - path to directory where the data products are stored (str)
            cal_scan - calibration scan number or 'default' is auto cal_type selected (str)
    OUTPUT: Nothing returned, but function creates,
            - File listing raw data files
    '''
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Get names
    if isinstance(num,list):
        names=[]
        for scan in num:
            names0 = glob.glob(data_dir+'s'+waveband+'*'+date+'_000'+scan+'_*.sdf')
            names0.sort()
            names.extend(names0)
        label='all_'+waveband
    else:
        names = glob.glob(data_dir+'s'+waveband+'*'+date+'_000'+num+'_*.sdf')
        names.sort()
        label=num+'_'+waveband
    # Write names into a file
    file = open(output_dir2+'mylist'+label+'.lst','w') 
    for name in names:
        file.write('{0}\n'.format(name)) # writes name in each line
    file.close()
    

def create_map_fits(num,configfile,cal_scan,output_dir,data_dir,crop_file,waveband):
    '''Creates a calibrated, cropped map, with and without matched filtering, of a SCUBA-2 scan. The following data reduction
    commands are performed in starlink:
    (1) Checks if in the data_products directory there is an existing folder 
       for each calibrator scan number or 'default' in the case of auto calibration.
       If not, it creates the folder
    (2) Creates a list with the names of raw data by calling create_list function
    (3) Runs smurf.makemap on the raw data
    (4) Gets the FCF value for the given cal_scan number or uses default values in cause of auto calibration.
    (5) Runs kappa.cmult to apply the FCF correction to the map, converting pW to
       Jy/beam. The output is a new file with '_cal' at the end, indicating that
       the FCF was applied.
    (6) Runs picard.crop_scuba2_images on the output of step (5) to crop the map
    (7) Converts the output of step 6 into fits files

    INPUT:  num - target scan number (str)
            configfile - name of configuration file WITH full path (str)
            cal_scan - calibration scan number or 'default' if auto calibration selected (str)
            output_dir -  path to directory where the data products are stored (str)
            data_dir - path to directory where the raw data is stored (str)
            crop_file - cropping parameter file (str)
            waveband - '8' or '4' for 850um or 450um (str)
    OUTPUT: Nothing returned, but function creates,
            - Calibrated and cropped map files (sdf and fits)
            - Calibrated and cropped map files with matched filter applied (sdf and fits)
    '''
    final_outputdir=output_dir+'/data_products/'
    # Check if directory of cal_scan exists already, if not, create directory
    if not os.path.isdir(output_dir+'calibrate_'+cal_scan):
        os.mkdir(output_dir+'calibrate_'+cal_scan)
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Create list
    create_list(num,data_dir,output_dir,cal_scan,waveband)
    # create map
    smurf.makemap(in_='^'+output_dir2+'mylist'+num+'_'+waveband+'.lst',out=output_dir2+'scan_'+num+'_'+waveband,config = '^'+configfile)  
    # get FCF value
    if cal_scan=='default' and waveband=='8':
        fcf_val= 537
    elif cal_scan=='default' and waveband=='4':
        fcf_val= 491
    else:
        fcf_val = FCF(output_dir2+'scan_'+cal_scan+'_'+waveband+'.sdf')
    # run cmult
    kappa.cmult(in_=output_dir2+'scan_'+num+'_'+waveband+'.sdf',scalar = fcf_val, out = output_dir2+'scan_'+num+'_'+waveband+'_cal')
    # cropping map
    cropped_map = picard.crop_scuba2_images([output_dir2+'scan_'+num+'_'+waveband+'_cal.sdf'],recpars = crop_file)
    # move to Data products directory
    os.system('mv '+cropped_map.datafiles[0]+' '+output_dir2)
    #make a matched filter version
    mf_map=picard.scuba2_matched_filter([output_dir2+'scan_'+num+'_'+waveband+'_cal_crop.sdf'])
    os.system('mv '+[item for item in mf_map.datafiles if '_mf' in item][0]+' '+output_dir2+'scan_'+num+'_'+waveband+'_cal_crop_mf.sdf')
    # create fits file
    convert.ndf2fits(in_= output_dir2+'scan_'+num+'_'+waveband+'_cal_crop.sdf',out = output_dir2+'scan_'+num+'_'+waveband+'_cal_crop.fits')
    convert.ndf2fits(in_= output_dir2+'scan_'+num+'_'+waveband+'_cal_crop_mf.sdf',out = output_dir2+'scan_'+num+'_'+waveband+'_cal_crop_mf.fits')
    #copy all final processed maps to the data_products directory
    if cal_scan=='default':
        os.system('cp -r '+output_dir2+'scan_'+num+'_'+waveband+'_cal_crop.fits'+' '+final_outputdir)
        os.system('cp -r '+output_dir2+'scan_'+num+'_'+waveband+'_cal_crop_mf.fits'+' '+final_outputdir)

def create_mosaic(flag,target_scan_list,output_dir,cal_scan,waveband,mf):
    '''Creates a mosaic map/cube that combines individual scan maps/timing cubes for a
    SCUBA-2 observation. The following procedure is used:
    (1) Searches for all the calibrated target scan maps or timing cubes
    (2) Creates a list with the names that meet step (1) criteria, and writes them to a file
    (3) Reads the header keys 'LBOUND1,2,3' and 'NAXIS1,2,3' to obtain the
        size limits of the maps to be mosaiced
    (4) Runs kappa.wcsmosaic, producing a combined version of the individual maps/cubes of each scan
    INPUT:  flag - indicates whether we are applying the function to the cropped map or the timing cube (str; 'scan' or 'cube')
            target_scan_list- list of strings with the names of target scan numbers
            output_dir - path to directory where the data products are stored (str)
            cal_scan - calibration scan number or 'default' if auto calibration selected (str)
            mf - flag to indicate if mosaicing maps with matched filter applied ('y' or 'n'; str)
    OUTPUT: Nothign returned, but function creates,
            - mosaiced map/cube (sdf and fits files)
    '''
    final_outputdir=output_dir+'/data_products/'
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Get the name of cropped files
    lsts=[]
    if flag == 'scan':
        if mf=='n':
            names_mosaic0 = glob.glob(output_dir2+'*'+waveband+'_cal_crop.sdf')
            outs=output_dir2+'mosaic_map_'+flag+'_'+waveband
        elif mf=='y':
            names_mosaic0 = glob.glob(output_dir2+'*'+waveband+'_cal_crop_mf.sdf')
            outs=output_dir2+'mosaic_map_'+flag+'_'+waveband+'_mf'
    if flag == 'cube':
        names_mosaic0 = glob.glob(output_dir2+'*'+waveband+'_shortmp_cube_cal.sdf')
        outs=output_dir2+'mosaic_map_'+flag+'_'+waveband
    names_mosaic = []
        
    for item in names_mosaic0:
        if any(x in item for x in target_scan_list):
            names_mosaic.append(item)
                
    names_mosaic.sort()
    # Write names into a file
    fileo = open(output_dir2+'mylist_mosaic_'+flag+'_'+waveband+'.lst','w') 
    for name in names_mosaic:
        fileo.write('{0}\n'.format(name))
    fileo.close()
    # Get image size, first convert names_mosaic into .fits file
    fitshdr=fits.open(names_mosaic[0].split('.sdf')[0]+'.fits')[0].header
    lbound1 = int(fitshdr['LBOUND1'])
    lbound2 = int(fitshdr['LBOUND2'])
    lbound3 = int(fitshdr['LBOUND3'])
    naxis1 = int(fitshdr['NAXIS1'])
    naxis2 = int(fitshdr['NAXIS2'])
    naxis3 = int(fitshdr['NAXIS3'])

    # Run wcsmosaic
    kappa.wcsmosaic(in_='^'+output_dir2+'mylist_mosaic_'+flag+'_'+waveband+'.lst', out = outs, lbnd = [lbound1,lbound2,lbound3] , ubnd = [lbound1+naxis1,lbound2+naxis2,lbound3+naxis3],ref = names_mosaic[0])
    # Convert to fits
    convert.ndf2fits(in_=outs+'.sdf',out = outs+'.fits')
    #copy all final processed maps to the data_products directory
    if cal_scan=='default':
        os.system('cp -r '+outs+'.fits '+final_outputdir)
    
def timing_cube(num,configfile,cal_scan,output_dir,waveband):
    '''Creates a calibrated timing cube (RA,DEC,time) of a SCUBA-2 scan. The following data reduction
    commands are performed in starlink:
    (1) Checks if in the data_products directory there is an existing folder 
       for each calibrator scan number or 'default' in the case of auto calibration.
       If not, it creates the folder
    (2) Runs smurf.makemap on the raw data for the specified scan (the raw data list file was already
        created in the full scan map function)
    (3) Runs smurf.stackframes to create a cube
    (4) Gets the FCF value for the given cal_scan number or uses default values in cause of auto calibration.
    (5) Runs kappa.cmult to apply the FCF correction to the map, converting pW 
        to Jy/beam. The output is a new file with '_shortmp_cube_cal' at the 
        end,indicating that the FCF was applied.
    (6) Converts the output of step 5 into fits files
    INPUT:  num - target scan number (str)
            configfile - name of configuration file, WITH the full path (str)
            cal_scan - calibration scan number or 'default' if auto calibration selected (str)
            output_dir - path to directory where the data products are stored (str)
    OUTPUT: Nothing returned, but function creates,
            - Timing cube (RA,DEC,time; sdf and fits files)
    '''
    final_outputdir=output_dir+'/data_products/'
    # Check if directory of cal_scan exists already, if not, create directory
    if not os.path.isdir(output_dir+'calibrate_'+cal_scan):
        os.mkdir(output_dir+'calibrate_'+cal_scan)
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # making timing cube 
    smurf.makemap(in_='^'+output_dir2+'mylist'+num+'_'+waveband+'.lst',out=output_dir2+'scan_'+num+'_'+waveband+'_shortmp',config = '^'+configfile)
    # stacking frames
    smurf.stackframes(in_ = output_dir2+'scan_'+num+'_'+waveband+'_shortmp.more.smurf.shortmaps', out = output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube', sort = False, sortby = '')
    # get FCF val
    if cal_scan=='default' and waveband=='8':
        fcf_val= 537
    elif cal_scan=='default' and waveband=='4':
        fcf_val= 491
    else:
        fcf_val = FCF(output_dir2+'scan_'+cal_scan+'_'+waveband+'.sdf')
    # run cmult
    kappa.cmult(in_=output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube.sdf',scalar = fcf_val, out = output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube_cal')
    # create .fits file
    convert.ndf2fits(in_=output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube_cal.sdf',out = output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube_cal.fits')
    #copy all final processed maps to the data_products directory
    if cal_scan=='default':
        os.system('cp -r '+ output_dir2+'scan_'+num+'_'+waveband+'_shortmp_cube_cal.fits '+final_outputdir)
   
def noise_estimate_scuba2(num,output_dir,cal_scan,flag,single,waveband):
    '''Creates a log file that contains the rms noise estimation for a calibrated SCUBA-2 map
    INPUT:  num - target scan number (str)
            output_dir - path to directory where the data products are stored (str)
            cal_scan - calibration scan number or 'default' if auto calibration selected (str)
            flag - label for indicating whether the rms noise is calculated from a calibrator ('cal') or target ('targ') map
            single - indicates if the map is made from a single scan ('scan') or mosaic of scans ('mosaic'; str)
    OUTPUT: Nothing returned, but function creates,
            - A .txt file with 
                column 1: target scan number
                column 2: flag indicating type of file (cal or targ)
                column 3: noise estimation in mJy
                column 4: maximum intensity in mJy
                column 5: x-coordinate of pixel with max intensity
                column 6: y-coordinate of pixel with max intensity
                column 7: beam major axis (fwhm) in arcsec calculted from calibrator map or defualt values
                column 8: beam minor axis (fwhm) in arcsec calculted from calibrator map or defualt values
                column 9: beam position angle in deg calculted from calibrator map or 0.
    ''' 
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    if single=='scan':
        filename = output_dir2+'scan_'+num+'_'+waveband+'_cal_crop.sdf'
    elif single=='mosaic':
        filename = output_dir2+'mosaic_map_scan_'+waveband+'.sdf'
    rms_noise = kappa.stats(filename,comp='err')
    noise_val = rms_noise.mean*1e3 # uncertainty in mJy
    max_intensity = kappa.stats(filename)
    max_val = max_intensity.maximum*1e3
    # read fits file corresponding to output_dir2+'scan_'+num+'_cal_crop
    fits_map = fits.open(filename.strip('.sdf')+'.fits')[0].data
    # find the indices corresponding to the maximum flux
    indices = np.unravel_index(np.argmax(np.nan_to_num(fits_map), axis=None),fits_map.shape)
    x = indices[1]
    y = indices[2]
    # get the coordinates of the maximum intensity in the calibrator scan, to be
    # used later to fit the beam size in that map
    if cal_scan=='default' and waveband=='8':
        majfwhm=14.6
        minfwhm=14.6
        orient=0.0
    elif cal_scan=='default' and waveband=='4':
        majfwhm=9.8
        minfwhm=9.8
        orient=0.0
    else:
        cal_stat = kappa.stats(output_dir2+'scan_'+cal_scan+'_'+waveband+'_cal_crop.sdf').maxwcs.split(',')
        cal_beam = kappa.beamfit(output_dir2+'scan_'+cal_scan+'_'+waveband+'_cal_crop.sdf', pos = '"'+cal_stat[0]+','+cal_stat[1]+'"',mode='interface', gauss = False)
        majfwhm = Angle(cal_beam.majfwhm[0]*u.rad).arcsec
        minfwhm = Angle(cal_beam.minfwhm[0]*u.rad).arcsec
        orient = cal_beam.orient[0]
    # Keep record of noise_val for each target number in a file
    if os.path.isfile(output_dir2+'noise_log'+waveband+'.txt'):       
        fileo = open(output_dir2+'noise_log'+waveband+'.txt','a') 
    elif not os.path.isfile(output_dir2+'noise_log'+waveband+'.txt'): 
        fileo = open(output_dir2+'noise_log'+waveband+'.txt','w')
        fileo.write('#Scan type rms max_I X Y beammaj beammin bpa\n')
        fileo.write('#() () (mJy) (mJy) (pix) (pix) (arcsec) (arcsec) (deg)\n')
    fileo.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(num,flag,noise_val,max_val,x,y,majfwhm,minfwhm,orient))
    fileo.close()

def create_pol_map(nums,data_dir,output_dir,date,configfile,crop_file,waveband):
    '''Creates calibrated, cropped stokes cubes of a POL-2 scan. The following data reduction
    commands are performed in starlink:
    (1) Checks if in the data_products directory there is an existing folder 
       for each calibrator scan number or 'default' in the case of auto calibration.
       If not, it creates the folder
    (2) Creates a list with the names of raw data by calling create_list function
    (3) Runs smurf.pol2map command on the raw polarization data to create a 
    total intensity map in pW. It is run again, with the first I map as input, to produce the
    final I, Q and U maps, and the vector catalogue.
    (3) Runs kappa.cmult to apply defualt FCFs to IQU maps (individual scans and mosaics)
    (4) Runs picard.crop_scuba2_images to crop IQU maps (individual scans and mosaics)
    (5) IQU maps are converted to fits files, and a stokes cube is built (individual scans and mosaics)
    INPUT:  nums - list of strings of target scan numbers
            data_dir - path to directory where raw data is stored (str)
            output_dir - path to directory where the data products are stored (str)
            date - date of observation (YYYMMDD, str)
            configfile - name of configuration file WITH full path (str)
            crop_file - name of cropping parameter file WITH full path (str)
            waveband - '8' or '4' for 850um or 450um (str)
    OUTPUT: Nothing returned, but function creates,
            - Calibrated, cropped stokes cubes (individual scans and mosaics with and without matched filters; fits files)
    '''
    final_outputdir=output_dir+'/data_products/'
    if waveband=='8':
        fcf_val=725
    elif waveband=='4':
        fcf_val=962
    if not os.path.isdir(output_dir+'calibrate_default'):
        os.mkdir(output_dir+'calibrate_default')
    output_dir2 = output_dir+'calibrate_default/'
    create_list(nums,data_dir,output_dir,'default',waveband)
    smurf.pol2map(in_ = '^'+output_dir2+'mylistall_'+waveband+'.lst', iout = output_dir2+'stokes_i_'+waveband+'/'+'mosaic_Imap',qout='!',\
        uout='!', mapdir = output_dir2+'stokes_i_'+waveband, qudir = output_dir2+'stokes_qu_'+waveband,jy=True, fcf='!',skyloop=False,\
        config = '^'+configfile)
    smurf.pol2map(in_=output_dir2+'stokes_qu_'+waveband+'/*',iout = output_dir2+'stokes_i_'+waveband+'/'+'mosaic_Imap_2',
              qout=output_dir2+'stokes_qu_'+waveband+'/'+'mosaic_Qmap', uout=output_dir2+'stokes_qu_'+waveband+'/'+'mosaic_Umap',
              mapdir=output_dir2+'stokes_i_'+waveband+'/',mask=output_dir2+'stokes_i_'+waveband+'/'+'mosaic_Imap',\
              maskout1=output_dir2+'stokes_i_'+waveband+'/'+'astmask',
              maskout2=output_dir2+'stokes_i_'+waveband+'/'+'pcamask',ipref=output_dir2+'stokes_i_'+waveband+'/'+'mosaic_Imap_2',
              cat=output_dir2+'stokes_qu_'+waveband+'/'+'mycat',debias=True,jy=True, fcf='!',skyloop=False,\
              config =  '^'+configfile)
    #calibrate,crop and convert all output I,Q,U maps (indv scans and mosaics, both with matched filter versions) to fits
    for n in nums:
        list_pol=[glob.glob(output_dir2+'stokes_i_'+waveband+'/'+date+'_000'+n+'*_'+i+'map.sdf')[0] for i in ['I','Q','U']]
        for item in list_pol:
            kappa.cmult(in_=item,scalar = fcf_val, out = item.strip('.sdf')+'_cal.sdf')
            cropped_map = picard.crop_scuba2_images([item.strip('.sdf')+'_cal.sdf'],recpars = crop_file)
            os.system('mv '+cropped_map.datafiles[0]+' '+output_dir2+'stokes_i_'+waveband)
            mf_map=picard.scuba2_matched_filter([item.strip('.sdf')+'_cal_crop.sdf'])
            os.system('mv '+[i for i in mf_map.datafiles if '_mf' in i][0]+' '+item.strip('.sdf')+'_cal_crop_mf.sdf')
            convert.ndf2fits(in_=item.strip('.sdf')+'_cal_crop.sdf',out = item.strip('.sdf')+'_cal_crop.fits')
            convert.ndf2fits(in_=item.strip('.sdf')+'_cal_crop_mf.sdf',out = item.strip('.sdf')+'_cal_crop_mf.fits')
        lpf=[glob.glob(output_dir2+'stokes_i_'+waveband+'/'+date+'_000'+n+'*_'+i+'map_cal_crop.fits')[0] for i in ['I','Q','U']]
        stokes_cube=create_stokes_cubes(lpf,output_dir2,waveband,date,n,'n','n')
        lpf_mf=[glob.glob(output_dir2+'stokes_i_'+waveband+'/'+date+'_000'+n+'*_'+i+'map_cal_crop_mf.fits')[0] for i in ['I','Q','U']]
        stokes_cube_mf=create_stokes_cubes(lpf_mf,output_dir2,waveband,date,n,'n','y')
        #copy all final processed maps to the data_products directory
        os.system('cp -r '+stokes_cube+' '+final_outputdir)
        os.system('cp -r '+stokes_cube_mf+' '+final_outputdir)
    list_pol_mosaic=glob.glob(output_dir2+'stokes_*'+waveband+'/mosaic_*map.sdf')
    for item in list_pol_mosaic:
        kappa.cmult(in_=item,scalar = fcf_val, out = item.strip('.sdf')+'_cal.sdf')
        cropped_map = picard.crop_scuba2_images([item.strip('.sdf')+'_cal.sdf'],recpars = crop_file)
        os.system('mv '+cropped_map.datafiles[0]+' '+output_dir2+item.split('/')[-2])
        mf_map=picard.scuba2_matched_filter([item.strip('.sdf')+'_cal_crop.sdf'])
        os.system('mv '+[i for i in mf_map.datafiles if '_mf' in i][0]+' '+item.strip('.sdf')+'_cal_crop_mf.sdf')
        convert.ndf2fits(in_=item.strip('.sdf')+'_cal_crop.sdf',out = item.strip('.sdf')+'_cal_crop.fits')
        convert.ndf2fits(in_=item.strip('.sdf')+'_cal_crop_mf.sdf',out = item.strip('.sdf')+'_cal_crop_mf.fits')
    lpfmos=[glob.glob(output_dir2+'stokes_*'+waveband+'/mosaic_'+str(i)+'map_cal_crop.fits')[0] for i in ['I','Q','U']]
    stokes_cube_mosaic=create_stokes_cubes(lpfmos,output_dir2,waveband,date,'','y','n')
    lpfmos_mf=[glob.glob(output_dir2+'stokes_*'+waveband+'/mosaic_'+str(i)+'map_cal_crop_mf.fits')[0] for i in ['I','Q','U']]
    stokes_cube_mosaic_mf=create_stokes_cubes(lpfmos_mf,output_dir2,waveband,date,'','y','y')
    #copy all final processed maps to the data_products directory
    os.system('cp -r '+stokes_cube_mosaic+' '+final_outputdir)
    os.system('cp -r '+stokes_cube_mosaic_mf+' '+final_outputdir)

def create_stokes_cubes(lists,output_dir2,waveband,date,n,mos,mf):
    '''Builds a stokes cube from a list of I,Q,U maps
    INPUT:  lists - list of strings of I,Q,U fits files
            output_dir2 - output directory (str)
            waveband - '8' or '4' for 850um or 450um (str)
            date - date of observation (YYYMMDD, str)
            n - scan number in case of individual maps (str)
            mos - flag indicating if the I,Q,U maps are mosaics ('y' or 'n'; str)
            mf - flag indicating if the I,Q,U maps are matched filtered ('y' or 'n'; str)
    OUTPUT: Function returns filename of output stokes cube, and it creates,
            - Stokes cube of input maps (planes I,Q,U, extensions DATA,VARIANCE)
    '''
    list_pol_fits=[]
    list_pol_fitsvar=[]
    for fname in lists:
        list_pol_fits.append(getdata(fname,0))
        list_pol_fitsvar.append(getdata(fname,1))
    list_pol_fits_arr=np.array(list_pol_fits)
    list_pol_fits_arrvar=np.array(list_pol_fitsvar)
    if mos=='n':
        mname=output_dir2+'stokes_i_'+waveband+'/'+date+'_000'+n+'_stokes_cube'
    elif mos=='y':
        mname=output_dir2+'stokes_i_'+waveband+'/mosaic_map_stokes_cube'
    if mf=='n':
        label='.fits'
    elif mf=='y':
        label='_mf.fits'
    fits.writeto(mname+label, list_pol_fits_arr,header=getheader(lists[0]))
    fits.writeto(mname+'_var'+label, list_pol_fits_arrvar,header=getheader(lists[0]))
    hdu1 = fits.PrimaryHDU(data=list_pol_fits_arr,header=getheader(lists[0]))
    hdu2 = fits.ImageHDU(data=list_pol_fits_arrvar)
    new_hdu = fits.HDUList([hdu1, hdu2])
    new_hdu.writeto(mname+'_full'+label, overwrite=True)
    return(mname+'_full'+label)

def create_pol_timing(nums,data_dir,output_dir,date,configfile,waveband):
    '''Creates calibrated, cropped I,Q,U timing cubes of a POL-2 scan. The following data reduction
    commands are performed in starlink:
    (1) Checks if in the data_products directory there is an existing folder 
       for each calibrator scan number or 'default' in the case of auto calibration.
       If not, it creates the folder
    (2) Creates a list with the names of raw data by calling create_list function
    (3) Runs smurf.pol2map command on the raw polarization data to create a 
    total intensity cube in pW. It is run again, with the first I cube as input, to produce the
    final I, Q and U cubes.
    (3) Runs kappa.cmult to apply defualt FCFs to IQU cubes (individual scans only)
    (4) IQU cubes are converted to fits files (individual scans only)
    INPUT:  nums - list of strings of target scan numbers
            data_dir - path to directory where raw data is stored (str)
            output_dir - path to directory where the data products are stored (str)
            data - date of observation (YYYMMDD, str)
            configfile - name of configuration file WITH full path (str)
            crop_file - name of cropping parameter file WITH full path (str)
            waveband - '8' or '4' for 850um or 450um (str)
    OUTPUT: Nothing returned, but function creates,
            - Calibrated I,Q,U timing cubes (individual scans only) 
    '''
    final_outputdir=output_dir+'/data_products/'
    if waveband=='8':
        fcf_val=725
    elif waveband=='4':
        fcf_val=962
    output_dir2 = output_dir+'calibrate_default/'
    smurf.pol2map(in_ = '^'+output_dir2+'mylistall_'+waveband+'.lst', iout = output_dir2+'stokes_icube_'+waveband+'/'+'mosaic_Imap',qout='!',\
        uout='!', mapdir = output_dir2+'stokes_icube_'+waveband, qudir = output_dir2+'stokes_qucube_'+waveband,jy=True, fcf='!',skyloop=False,\
        config = '^'+configfile)
    smurf.pol2map(in_=output_dir2+'stokes_qucube_'+waveband+'/*',iout = output_dir2+'stokes_icube_'+waveband+'/'+'mosaic_Imap_2',
              qout=output_dir2+'stokes_qucube_'+waveband+'/'+'mosaic_Qmap', uout=output_dir2+'stokes_qucube_'+waveband+'/'+'mosaic_Umap',
              mapdir=output_dir2+'stokes_icube_'+waveband+'/',mask=output_dir2+'stokes_icube_'+waveband+'/'+'mosaic_Imap',\
              maskout1=output_dir2+'stokes_icube_'+waveband+'/'+'astmask',
              maskout2=output_dir2+'stokes_icube_'+waveband+'/'+'pcamask',ipref=output_dir2+'stokes_icube_'+waveband+'/'+'mosaic_Imap_2',
              cat=output_dir2+'stokes_qucube_'+waveband+'/'+'mycat',debias=True,jy=True, fcf='!',skyloop=False,\
              config =  '^'+configfile)
    #calibrate and convert all output I,Q,U maps (only indv scans have cubes) to fits
    for n in nums:
        #need to cut name length down as starlink gives truncation errors for stackframes task!
        for i in ['I','Q','U']:
            os.system('mv '+output_dir2+'stokes_icube_'+waveband+'/'+date+'_000'+n+'*_'+i+'map.sdf'+' '+output_dir2+'stokes_icube_'+waveband+'/'+'scan_'+n+'_'+i+'map.sdf')
        list_pol=[glob.glob(output_dir2+'stokes_icube_'+waveband+'/'+'scan_'+n+'_'+i+'map.sdf')[0] for i in ['I','Q','U']]
        for item in list_pol:
            smurf.stackframes(in_ = item.strip('.sdf')+'.more.smurf.shortmaps',out = item.strip('.sdf')+'_cube', sort = False, sortby = '')
            kappa.cmult(in_=item.strip('.sdf')+'_cube.sdf',scalar = fcf_val, out = item.strip('.sdf')+'_cube_cal.sdf')
            convert.ndf2fits(in_=item.strip('.sdf')+'_cube_cal.sdf',out = item.strip('.sdf')+'_cube_cal.fits')
            #copy all final processed maps to the data_products directory
            os.system('cp -r '+item.strip('.sdf')+'_cube_cal.fits '+final_outputdir)

def noise_estimate_pol2(nums,output_dir,waveband):
    '''Creates a log file that contains the rms noise estimation for calibrated POL-2 I, Q, U target
    scan maps (individual scans and mosaics).
    INPUT:  nums - list of strings of target scan numbers
            output_dir - path to directory where the data products are stored (str)
            waveband - '8' or '4' for 850um or 450um
    OUTPUT: Nothing returned, but function creates,
            - A .txt file with 
                column 1: target scan number
                column 2: flag indicating type of file (cal or targ)
                column 3: noise estimation in mJy
                column 4: maximum intensity in mJy
                column 5: x-coordinate of pixel with max intensity
                column 6: y-coordinate of pixel with max intensity
                column 7: beam major axis (fwhm) in arcsec calculted from calibrator map or defualt values
                column 8: beam minor axis (fwhm) in arcsec calculted from calibrator map or defualt values
                column 9: beam position angle in deg calculted from calibrator map or 0.
    ''' 
    output_dir2 = output_dir+'calibrate_default/'
    files=[]
    flag=[]
    scans=[]
    for n in nums:
        files.extend([glob.glob(output_dir2+'stokes_i_'+waveband+'/'+date+'_000'+n+'*_'+i+'map_cal_crop.sdf')[0] for i in ['I','Q','U']])
        flag.extend(['I','Q','U'])
        scans.extend([n,n,n])
    files.append(output_dir2+'stokes_i_'+waveband+'/mosaic_Imap_cal_crop.sdf')
    files.append(output_dir2+'stokes_qu_'+waveband+'/mosaic_Qmap_cal_crop.sdf')
    files.append(output_dir2+'stokes_qu_'+waveband+'/mosaic_Umap_cal_crop.sdf')
    flag.extend(['I','Q','U'])
    scans.extend(['mosaic','mosaic','mosaic'])
    for i in range(0,len(files)):
        rms_noise = kappa.stats(files[i],comp='err')
        noise_val = rms_noise.mean*1e3 # uncertainty in mJy
        max_intensity = kappa.stats(files[i])
        max_val = max_intensity.maximum*1e3
        fits_map = fits.open(files[i].strip('.sdf')+'.fits')[0].data
        indices = np.where(fits_map==max_val/1000) # returns z,x,y
        x = indices[1][0]
        y = indices[2][0]
        if waveband=='8':
            majfwhm=14.6
            minfwhm=14.6
            orient=0.0
        elif waveband=='4':
            majfwhm=9.8
            minfwhm=9.8
            orient=0.0
        if os.path.isfile(output_dir2+'noise_log'+waveband+'.txt'):       
            fileo = open(output_dir2+'noise_log'+waveband+'.txt','a') 
        elif not os.path.isfile(output_dir2+'noise_log'+waveband+'.txt'): 
            fileo = open(output_dir2+'noise_log'+waveband+'.txt','w')
            fileo.write('#Scan type rms_IQU max_IQU X Y beammaj beammin bpa\n')
            fileo.write('#() () (mJy) (mJy) (pix) (pix) (arcsec) (arcsec) (deg)\n')
        fileo.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(scans[i],flag[i],noise_val,max_val,x,y,majfwhm,minfwhm,orient))
        fileo.close()

def scuba2_reduce(data_dir,output_dir,cal_scan_numbers,target_scan_numbers,cal_type,file_crop,file_config,file_config_timing,waveband,integ):
    '''SCUBA-2 reduction procedure
    INPUT:  data_dir - path to directory where the raw data files are stored (str)
            output_dir - path to directory where the data products are stored (str)
            cal_scan_numbers - list of strings of calibration scan numbers or ['default'] if auto calibration selected
            target_scan_numbers - list of strings of target scan numbers
            cal_type - 'user' (outputs list of all cal scans present for manual FCF determination) or
            'auto' (outputs ['default'] for application of standard FCFs; str)
            file_crop - name of cropping parameter file WITH full path (str)
            file_config - name of configuration file WITH full path for mapmaker single scan map creation (str)
            file_config_timing - name of configuration file WITH full path for mapmaker cube creation (str)
            waveband - '8' or '4' for 850um or 450um (str)
            integ - timebins for lightcurves (0 for no timing cubes created, float)
    OUTPUT: Nothing returned, but function creates,
            - Calibrated, cropped maps for each scan and mosaic of scans, with and without matched filters (sdf and fits files)
            - Calibrated timing cubes for each scan and mosaic of scan cubes (if integ>0; sdf and fits files)
            - rms noise estimation log
    '''
    for item in cal_scan_numbers:
        if cal_type=='user' and item !='default':
            create_map_fits(item,file_config,item,output_dir,data_dir,file_crop,waveband)
            noise_estimate_scuba2(item,output_dir,item,'cal','scan',waveband)
        for targ in target_scan_numbers:
            create_map_fits(targ,file_config,item,output_dir,data_dir,file_crop,waveband)
            noise_estimate_scuba2(targ,output_dir,item,'targ','scan',waveband)
            if integ != 0:
                timing_cube(targ,file_config_timing,item,output_dir,waveband)
        create_mosaic('scan',target_scan_numbers,output_dir,item,waveband,'n')
        create_mosaic('scan',target_scan_numbers,output_dir,item,waveband,'y')
        noise_estimate_scuba2('all',output_dir,item,'mosaic','mosaic',waveband)
        if integ != 0:
            create_mosaic('cube',target_scan_numbers,output_dir,item,waveband,'n')
def pol2_reduce(data_dir,output_dir,target_scan_numbers,file_crop,file_config,file_config_timing,waveband,integ):
    '''POL-2 reduction procedure
    INPUT:  data_dir - path to directory where the raw data files are stored (str)
            output_dir - path to directory where the data products are stored (str)
            target_scan_numbers - list of strings of target scan numbers
            file_crop - name of cropping parameter file WITH full path (str)
            file_config - name of configuration file WITH full path for mapmaker single scan map creation (str)
            file_config_timing - name of configuration file WITH full path for mapmaker cube creation (str)
            waveband - '8' or '4' for 850um or 450um (str)
            integ - timebins for lightcurves (0 for no timing cubes created, float)
    OUTPUT: Nothing returned, but function creates,
            - Calibrated, cropped stokes cubes for each scan and mosaic of scans, with and without matched filters (sdf and fits files)
            - Calibrated timing I, Q, U cubes for each scan (if integ>0; sdf and fits files)
            - rms noise estimation log 
    '''
    create_pol_map(target_scan_numbers,data_dir,output_dir,date,file_config,file_crop,waveband)
    noise_estimate_pol2(target_scan_numbers,output_dir,waveband)
    if integ!=0:
        create_pol_timing(target_scan_numbers,data_dir,output_dir,date,file_config_timing,waveband)
def inputs2file(output_dir,plist,plist_names):
    '''Dumps all user input parameters into a file for use by the analysis script later.
    INPUT:  output_dir - directory to save file in (str)
            plist - list of parameters
            plist_names - list of strings of the names of the parameters in plist
    OUTPUT: Returns nothing, but creates,
            - A file (shared.pkl) containing all parameters in plist
    '''
    shared={}
    for i in range(0,len(plist)):
        shared[plist_names[i]]=plist[i]
    with open(output_dir+'shared.pkl','wb') as fp:
        pickle.dump(shared, fp)

#########################  
if __name__=="__main__":
    ######################### 
    #Parameter section
    ######################### 
    date = '20160125'#20150622'#
    target = '3c273'#'BHXRB V404 Cyg'
    cal_type= 'user'#'user' or auto', always auto for pol2
    obs_type='pol2'#'scuba2 or pol2'
    waveband='8'#'4'
    integ=120#time bins in seconds, a 0 indicates no timing cubes created
    crop_params=200.0#radius in arcsec to crop output maps of full scans
    wout='/export/data2/atetarenko/PB_test/tutorial/'
    #########################

    #define input/output directories
    data_dir =wout+'raw/'+date+'/'
    output_dir = wout+'results/'+date+'/'

    #make/clear out the output diretory for the observation date
    if os.path.isdir(output_dir):
        os.system('rm -rf '+output_dir+'*')
    else:
        os.mkdir(output_dir)
    #create a directory within output_dir where all processed data products will go
    if not os.path.isdir(output_dir+'/data_products'):
        os.mkdir(output_dir+'/data_products')

    #write script inputs to a file for use by the analysis script later
    plist=[date,target,cal_type,obs_type,waveband,integ,crop_params]
    plist_names=['date','target','cal_type','obs_type','waveband','integ','crop_params']
    inputs2file(output_dir,plist,plist_names)

    #grab and write config files for map cropping recipe and mapmaker/pol2map
    file_crop=create_crop_config(output_dir,crop_params)
    file_config,file_config_timing=create_configs(wrapper.starpath,output_dir,obs_type,integ)

    #create lists of target and calibrator scan numbers
    target_scan_numbers,cal_scan_numbers = get_scan_numbers(data_dir,date,target,cal_type,waveband)

    if obs_type=='scuba2':
        scuba2_reduce(data_dir,output_dir,cal_scan_numbers,target_scan_numbers,cal_type,file_crop,file_config,file_config_timing,waveband,integ)
    elif obs_type=='pol2':
        pol2_reduce(data_dir,output_dir,target_scan_numbers,file_crop,file_config,file_config_timing,waveband,integ)
        #move pol2map logfiles out of working directory
        os.system('mv pol2map.log* '+output_dir)
    #move the noise log (from auto calibration) into the final processing directory
    os.system('cp -r '+output_dir+'calibrate_default/noise_log'+waveband+'.txt '+output_dir+'/data_products')
    #delete the temporary directory(s) created in the current working directory by starlink
    os.system('rm -rf PICARD*')
    os.system('rm -rf tmp*')
