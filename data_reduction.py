#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 10:23:55 2020

@author: constanza
"""

from starlink import wrapper
from starlink import kappa
from starlink import smurf
from starlink import picard
from starlink import convert
wrapper.change_starpath('/Users/constanza/star-2018A')
import glob
import numpy as np
from astropy.io import ascii
import os
from astropy.io import fits

######################### Defining intensity functions
def get_scan_numbers(data_dir,date,target,cal_name_list):
    '''
    This function gets the scan numbers from the header of 'science' and 'scan' files
    ---------------------------------------------------------------------------------
    Parameters: - data_dir: path to folder where the data is stored
                - date: date of observation
                - target: string with the name of the target
                - cal_name_list: list of strings with the names of calibrators
    Output: - list of strings with target scan numbers
            - list of strings with calibration scan numbers
    '''
    names = glob.glob(data_dir+'*.sdf')
    scan_numbers = [i.split('_')[1].lstrip('000') for i in names]
    unique_scan_numbers = np.unique(np.array(scan_numbers))
    unique_scan_numbers_final = [str(item).zfill(2) for item in unique_scan_numbers]
    final_target_scan_numbers = []
    final_cal_scan_numbers = []
    # convert to fits to read header
    for num in unique_scan_numbers_final:
        os.system('rm -rf '+data_dir+'s8a'+date+'_000'+num+'_0001.fits')
        convert.ndf2fits(in_=data_dir+'s8a'+date+'_000'+num+'_0001.sdf',out = data_dir+'s8a'+date+'_000'+num+'_0001.fits')
        obj = fits.open(data_dir+'s8a'+date+'_000'+num+'_0001.fits')[0].header['OBJECT']
        mode = fits.open(data_dir+'s8a'+date+'_000'+num+'_0001.fits')[0].header['OBS_TYPE']
        if obj == target and mode == 'science':
            final_target_scan_numbers.append(num)
        elif obj in cal_name_list and mode == 'science':
            final_cal_scan_numbers.append(num)
    return final_target_scan_numbers, final_cal_scan_numbers

def check_data():
    '''
    This function checks if there are 'damaged' files before reducing all the data
    ------------------------------------------------------------------------------
    '''   
    foo=2
            
def FCF(cal_name):
    '''
    This function calculates the FCF value for a given calibrator (cal_name)
    ------------------------------------------------------------------------
    Parameters: - cal_name: string with the name of the calibrator
    Output: - a number corresponding to the FCF
    '''
    result_FCF = picard.scuba2_check_cal([cal_name])
    fcffile=result_FCF.logfiles[2]
    FCF8 = float(ascii.read(fcffile)['col12']) # gets the FCF value
    return FCF8

def create_list(num,data_dir,output_dir,cal_scan):
    '''
    This function creates a text file named 'mylist'+target scan number
    and saves it into the Data_products > 'calibrate_cal scan number' directory.
    Inside the text file there is a list with the names of raw data.
    ----------------------------------------------------------------------------
    Parameters: - num: string with target scan number
                - data_dir: path to folder where the data is stored
                - output_dir: path to folder where the data products are stored
                - cal_scan: string with calibration scan number
    Output: - creates a text file named 'mylist'+target scan number
    '''
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Get names
    names = glob.glob(data_dir+'s8*'+date+'_000'+num+'_*.sdf')
    names.sort()
    # Write names into a file
    file = open(output_dir2+'mylist'+num+'.lst','w') 
    for name in names:
        file.write('{0}\n'.format(name)) # writes name in each line
    file.close()
    

def create_map_fits(num,configfile,cal_scan,output_dir,data_dir):
    '''
    This function creates a cropped map of the raw data by running a series of
    starlink commands. It goes as follows:
    (1) Checks if in the Data_products directory, there is an existing folder 
       for each calibrator scan number. If not, it creates the folder
    (2) Creates a list with the names of raw data by calling create_list function
    (3) Runs smurf.makemap on the raw data
    (4) Gets the FCF value for the given cal_scan number
    (5) Runs kappa.cmult to apply the FCF correction to the map, converting pW to
       Jy/beam. The output is a new file with '_cal' at the end, indicating that
       the FCF was applied.
    (6) Runs picard.crop_scuba2_images on the output of step (5) to crop the map
    (7) Moves the file created in step (6) to the Data_products > 'calibrate_cal 
    scan number' directory.
    ----------------------------------------------------------------------------
    Parameters: - num: string with target scan number
                - configfile: name of configuration file
                - cal_scan: string with tcalibration scan number
                - output_dir: path to folder where the data products are stored
                - data_dir: path to folder where the data is stored
    Output: - calibrated and cropped map files
    '''
    # Check if directory of cal_scan exists already, if not, create directory
    if not os.path.isdir(output_dir+'calibrate_'+cal_scan):
        os.mkdir(output_dir+'calibrate_'+cal_scan)
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Create list
    create_list(num,data_dir,output_dir,cal_scan)
    # create map
    smurf.makemap(in_='^'+output_dir2+'mylist'+num+'.lst',out=output_dir2+'scan_'+num,config = '^'+data_dir+'dimmconfig_bright_compact.lis')  
    # get FCF value
    fcf_val = FCF(output_dir2+'scan_'+cal_scan+'.sdf')
    # run cmult
    kappa.cmult(in_=output_dir2+'scan_'+num+'.sdf',scalar = fcf_val, out = output_dir2+'scan_'+num+'_cal')
    # cropping map
    cropped_map = picard.crop_scuba2_images([output_dir2+'scan_'+num+'_cal.sdf'],recpars = data_dir+'crop_parameters.lis')
    # move to Data products directory
    os.system('mv '+cropped_map.datafiles[0]+' '+output_dir2)
    # create fits file
    convert.ndf2fits(in_= output_dir2+'scan_'+num+'_cal_crop.sdf',out = output_dir2+'scan_'+num+'_cal.fits')

def create_mosaic(flag,target_scan_list,output_dir,cal_scan):
    '''
    This function creates a map that combines the individual maps of each scan,
    by running a series of starlink commands. It goes as follows:
    (1) Goes into the calibrate_'+cal_scan number directory
    (2) Finds the names of the relevant files according to the flag
    (3) Searches for all the target scan numbers in the names of step (2)
    (4) Creates a list with the names that meet step (2) and step (3) criteria
    (5) Writes the names of step (5) into a file named 'mylist_mosaic_'+flag
    (6) Reads the header keys 'LBOUND1,2,3' and 'NAXIS1,2,3' to obtain the
        size limits
    (7) Runs kappa.wcsmosaic, the output is named 'mosaic_map_'+flag, which is
        a combined version of the individual maps of each scan
    ----------------------------------------------------------------------------
    Parameters: - flag: string indicating whether we are applying the function to
                  the cropped map or the timing cube
                - target_scan_list: list of strings with the names of target scan
                  numbers
                - output_dir: path to folder where the data products are stored
                - cal_scan: string with tcalibration scan number
    Output: - a file with combined maps of individual scans
    '''
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # Get the name of cropped files
    if flag == 'crop':
        names_mosaic0 = glob.glob(output_dir2+'*_crop.sdf')
                
    if flag == 'shortmp_cube':
        names_mosaic0 = glob.glob(output_dir2+'*_shortmp_cube_cal.fits')
    names_mosaic = []
        
    for item in names_mosaic0:
        if any(x in item for x in target_scan_list):
            names_mosaic.append(item)
                
    names_mosaic.sort()
    # Write names into a file
    file = open(output_dir2+'mylist_mosaic_'+flag+'.lst','w') 
    for name in names_mosaic:
        file.write('{0}\n'.format(name))
    file.close()
    # Get image size, first convert names_mosaic into .fits file
    lbound1 = int(kappa.fitsval(names_mosaic[0],'LBOUND1').value)
    lbound2 = int(kappa.fitsval(names_mosaic[0],'LBOUND2').value)
    lbound3 = int(kappa.fitsval(names_mosaic[0],'LBOUND3').value)
    naxis1 = int(kappa.fitsval(names_mosaic[0],'NAXIS1').value)
    naxis2 = int(kappa.fitsval(names_mosaic[0],'NAXIS2').value)
    naxis3 = int(kappa.fitsval(names_mosaic[0],'NAXIS3').value)

    # Run wcsmosaic
    kappa.wcsmosaic(in_='^'+output_dir2+'mylist_mosaic_'+flag+'.lst', out = output_dir2+'mosaic_map_'+flag, lbnd = [lbound1,lbound2,lbound3] , ubnd = [lbound1+naxis1,lbound2+naxis2,lbound3+naxis3],ref = names_mosaic[0])
    
def timing_cube(num,configfile,cal_scan,output_dir):
    '''
    This function stacks the scans to build a timing cube (RA,DEC,time). It 
    goes as follows:
    (1) Checks if the directory of cal_scan exists already, if not, it creates
        a new directory
    (2) Runs smurf.makemap on the raw data
    (3) Runs smurf.stackframes (joints frames)
    (4) Gets the FCF value for the given cal_scan number
    (5) Runs kappa.cmult to apply the FCF correction to the map, converting pW 
        to Jy/beam. The output is a new file with '_shortmp_cube_cal' at the 
        end,indicating that the FCF was applied.
    (6) Converts the output of step 5 in .fits files
    ----------------------------------------------------------------------------
    Parameters: - num: string with target scan number
                - configfile: name of configuration file
                - cal_scan: string with tcalibration scan number
                - output_dir: path to folder where the data products are stored
    Output: - A file corresponding to the timing cube (RA,DEC,time)
    '''
    # Check if directory of cal_scan exists already, if not, create directory
    if not os.path.isdir(output_dir+'calibrate_'+cal_scan):
        os.mkdir(output_dir+'calibrate_'+cal_scan)
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    # making timing cube 
    smurf.makemap(in_='^'+output_dir2+'mylist'+num+'.lst',out=output_dir2+'scan_'+num+'_shortmp',config = '^'+data_dir+'dimmconfig_bright_compact_shortmaps.lis')
    # stacking frames
    smurf.stackframes(in_ = output_dir2+'scan_'+num+'_shortmp.more.smurf.shortmaps', out = output_dir2+'scan_'+num+'_shortmp_cube', sort = False, sortby = '')
    # get FCF val
    fcf_val = FCF(output_dir2+'scan_'+cal_scan+'.sdf')
    # run cmult
    kappa.cmult(in_=output_dir2+'scan_'+num+'_shortmp_cube.sdf',scalar = fcf_val, out = output_dir2+'scan_'+num+'_shortmp_cube_cal')
    # create .fits file
    convert.ndf2fits(in_=output_dir2+'scan_'+num+'_shortmp_cube_cal.sdf',out = output_dir2+'scan_'+num+'_shortmp_cube_cal.fits')
    # Obtain start date
    #mjd_start = fits.open(output_dir+'scan_'+num+'_cal_shortmp_cube.fits')[0].header['MJD-OBS']
   
def noise_estimate(num,output_dir,cal_scan):
    '''
    This function creates a file that contain the noise estimation for each
    target scan number.
    ----------------------------------------------------------------------------
    Parameters: - num: string with target scan number
                - output_dir: path to folder where the data products are stored
                - cal_scan: string with tcalibration scan number
    Output: - A .txt file with the target scan number in column 1, and the noise
              estimation in column 2
    '''
    output_dir2 = output_dir+'calibrate_'+cal_scan+'/'
    rms_noise = kappa.stats(output_dir2+'scan_'+num+'_cal.sdf')
    noise_val = rms_noise.sigma*1e3 # uncertainty in mJy
    # Keep record of noise_val for each target number in a file
    if os.path.isfile(output_dir2+'noise_log.txt'):       
        file = open(output_dir2+'noise_log.txt','a') 
    elif not os.path.isfile(output_dir2+'noise_log.txt'): 
        file = open(output_dir2+'noise_log.txt','w') 
    file.write('{0} {1}\n'.format(num,noise_val)) # writes name in each line
    file.close()
######################### End of Defining Intensity functions   
    
######################### Defining Polarization functions

def create_pol_map(num,output_dir):
    '''
    This function runs pol2map command on the raw polarization data to create a 
    total intensity map in pW. It is ran again to produce the final I, Q and U 
    maps, and the vector catalogue
    ----------------------------------------------------------------------------
    Parameters: - num: string with target scan number
                - output_dir: path to folder where the data products are stored
    Output: 
    '''
    smurf.pol2map(in_ = '^'+output_dir+'mylist'+num+'.lst', iout = output_dir+'stokes_i/'+'total_intensity_'+num,qout='!', uout='!', mapdir = output_dir+'stokes_i', qudir = output_dir+'stokes_qu')
    smurf.pol2map(in_=output_dir+'stokes_qu'+'/*',iout = output_dir+'stokes_i/'+'total_intensity_'+num+'_part2',
              qout=output_dir+'stokes_qu/'+'q_intensity_'+num, uout=output_dir+'stokes_qu/'+'u_intensity_'+num,
              mapdir=output_dir+'stokes_i/',mask=output_dir+'stokes_i/'+'total_intensity_'+num,maskout1=output_dir+'stokes_i/'+'astmask',
              maskout2=output_dir+'stokes_i/'+'pcamask',ipref=output_dir+'stokes_i/'+'total_intensity_'+num+'_part2',
              cat=output_dir+'stokes_qu/'+'mycat',debias=True)

######################### End of Defining Polarization functions
    
    
    
######################### Parameter section

data_dir ='/Users/constanza/Desktop/Research/Data/'
output_dir = '/Users/constanza/Desktop/Research/Data_products/'
date = '20150622'
# NOTE: numbers must have 2 digits!
target = 'BHXRB V404 Cyg'
cal_name_list = ['Arp220','CRL2688']

config_file ='dimmconfig_bright_compact.lis'
config_file2 ='dimmconfig_bright_compact_shortmaps.lis'

target_scan_numbers,cal_scan_numbers = get_scan_numbers(data_dir,date,target,cal_name_list)

######################### End of Parameter section


#os.system('rm -rf '+output_dir+'*')

cal_scan_numbers = ['08']

for item in cal_scan_numbers:
    #create_map_fits(item,config_file,item,output_dir,data_dir)
    for targ in target_scan_numbers:
        #create_map_fits(targ,config_file,item,output_dir,data_dir)
        #timing_cube(targ,config_file,item,output_dir)
        noise_estimate(targ,output_dir,item)
    #create_mosaic('crop',target_scan_numbers,output_dir,item)
    #create_mosaic('shortmp_cube',target_scan_numbers,output_dir,item)

