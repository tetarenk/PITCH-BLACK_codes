# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from starlink import wrapper
from starlink import kappa
from starlink import smurf
from starlink import picard
from starlink import convert
from starlink.utilities import starhelp
wrapper.change_starpath('/Users/constanza/star-2018A')
import glob
import numpy as np
from astropy.io import ascii
import os
from astropy.io import fits

'''
Parameters
'''

data_dir ='/Users/constanza/Desktop/Research/Data/'
output_dir = '/Users/constanza/Desktop/Research/Data_products/'
date = '20150622'
scan8 = [36,37,41,42,46,47,52]
cal_scan_numbers = [8,32,54]
target = 'V404Cyg'

'''
end of parameter section
'''

'''
Create a file with the names of raw data
'''

names = glob.glob(data_dir+'s8*'+date+'_00037_*.sdf')
names.sort()

file1 = open(output_dir+'mylist37.lst','w') 
for name in names:
    file1.write('{0}\n'.format(name)) # writes name in each line
file1.close()


# Step 2
smurf.makemap(in_='^'+output_dir+'mylist37.lst',out=output_dir+'calibrator_scan_37',config = '^'+data_dir+'dimmconfig_bright_compact.lis')


# Step 4
result_FCF = picard.scuba2_check_cal([output_dir+'calibrator_scan_8.sdf'])
fcffile=result_FCF.logfiles[2]

FCF8 = float(ascii.read(fcffile)['col12']) # gets the FCF value

# Step 5 
kappa.cmult(in_=output_dir+'calibrator_scan_37.sdf',scalar = FCF8, out = output_dir+'calibrator_scan_37_cal')


'''
Cropping target maps
'''

cropped_map = picard.crop_scuba2_images([output_dir+'calibrator_scan_37_cal.sdf'],recpars = output_dir+'crop_parameters.lis')
os.system('rm -rf '+output_dir+cropped_map.datafiles[0])
os.system('mv '+cropped_map.datafiles[0]+' '+output_dir)

# Step 6
convert.ndf2fits(in_=output_dir+'calibrator_scan_37_cal.sdf',out = output_dir+'calibrator_scan_37_cal.fits')

# This is for step 7
names_mosaic = glob.glob(output_dir+'*_crop.sdf')
names_mosaic.sort()

file2 = open(output_dir+'mylist_mosaic.lst','w') 
for name in names_mosaic:
    file2.write('{0}\n'.format(name)) # writes name in each line
file2.close()

# in step 7 we need to specify the size of the image, so we need to read the header
convert.ndf2fits(in_=names_mosaic[0],out = names_mosaic[0]+'.fits')
lbound1 = fits.open(names_mosaic[0]+'.fits')[0].header['LBOUND1']
lbound2 = fits.open(names_mosaic[0]+'.fits')[0].header['LBOUND2']
lbound3 = fits.open(names_mosaic[0]+'.fits')[0].header['LBOUND3']
naxis1 = fits.open(names_mosaic[0]+'.fits')[0].header['NAXIS1']
naxis2 = fits.open(names_mosaic[0]+'.fits')[0].header['NAXIS2']
naxis3 = fits.open(names_mosaic[0]+'.fits')[0].header['NAXIS3']

# Step 7
kappa.wcsmosaic(in_='^'+output_dir+'mylist_mosaic.lst', out = output_dir+'mosaic_map', lbnd = [lbound1,lbound2,lbound3] , ubnd = [lbound1+naxis1,lbound2+naxis2,lbound3+naxis3],ref = names_mosaic[0])

'''
Noise estimation in map. Gives the uncertainty in the flux, to be used later
'''

rms_noise = kappa.stats(output_dir+'calibrator_scan_36_cal.sdf')
noise_val = rms_noise.mean*1e3 # unceratinty in mJy

'''
Making timing cube
'''
smurf.makemap(in_='^'+output_dir+'mylist36.lst',out=output_dir+'calibrator_scan_36_shortmp',config = '^'+data_dir+'dimmconfig_bright_compact_shortmaps.lis')

'''
Processing timing cube
'''

smurf.stackframes(in_ = output_dir+'calibrator_scan_36_shortmp.more.smurf.shortmaps', out = output_dir+'calibrator_scan_36_shortmp_cube', sort = False, sortby = '')
kappa.cmult(in_=output_dir+'calibrator_scan_36_shortmp_cube.sdf',scalar = FCF8, out = output_dir+'calibrator_scan_36_shortmp_cube_cal')

#convert to fits
convert.ndf2fits(in_=output_dir+'calibrator_scan_36_shortmp_cube.sdf',out = output_dir+'calibrator_scan_36_cal_shortmp_cube.fits')

mjd_start = fits.open(output_dir+'calibrator_scan_36_cal_shortmp_cube.fits')[0].header['MJD-OBS']


######## POLARIZATION DATA REDUCTION

'''
Parameters
'''

data_dir ='/Users/constanza/Desktop/Research/Data_pol/tutorial/raw/'
output_dir = '/Users/constanza/Desktop/Research/Data_products_pol/'
date = '20160125'
scan8 = [43]



'''
Create a file with the names of raw data
'''

names = glob.glob(data_dir+'s8*'+date+'_00043_*.sdf')
names.sort()

file1 = open(output_dir+'mylist43.lst','w') 
for name in names:
    file1.write('{0}\n'.format(name)) # writes name in each line
file1.close()

'''
create pol map - total intensity map in pW (calibrate later!)
'''

smurf.pol2map(in_ = '^'+output_dir+'mylist43.lst', iout = output_dir+'stokes_i/'+'total_intensity_43',
              qout='!', uout='!', mapdir = output_dir+'stokes_i', qudir = output_dir+'stokes_qu')
              #jy = True, fcf = '!' (!=default value)


'''
produce the final I, Q and U maps, and the vector catalogue
takes the product in stokes_i map to remove the polarization signal
produced by the instrument 
'''
smurf.pol2map(in_=output_dir+'stokes_qu'+'/*',iout = output_dir+'stokes_i/'+'total_intensity_43_part2',
              qout=output_dir+'stokes_qu/'+'q_intensity_43', uout=output_dir+'stokes_qu/'+'u_intensity_43',
              mapdir=output_dir+'stokes_i/',mask=output_dir+'stokes_i/'+'total_intensity_43',maskout1=output_dir+'stokes_i/'+'astmask',
              maskout2=output_dir+'stokes_i/'+'pcamask',ipref=output_dir+'stokes_i/'+'total_intensity_43_part2',cat=output_dir+'stokes_qu/'+'mycat',
              debias=True)



######## POLARIZATION TIMING -- ASK SOMEONE

smurf.pol2map(in_ = '^'+output_dir+'mylist43.lst', iout = output_dir+'stokes_i/'+'total_intensity_43',
              qout='!', uout='!', mapdir = output_dir+'stokes_i', qudir = output_dir+'stokes_qu', 
              config = data_dir+'dimmconfig_pol2_compact_shortmap.lis',jy=True)



smurf.pol2map(in_=output_dir+'stokes_qu'+'/*',iout = output_dir+'stokes_i/'+'total_intensity_43_part2',
              qout=output_dir+'stokes_qu/'+'q_intensity_43', uout=output_dir+'stokes_qu/'+'u_intensity_43',
              mapdir=output_dir+'stokes_i/',mask=output_dir+'stokes_i/'+'total_intensity_43',maskout1=output_dir+'stokes_i/'+'astmask',
              maskout2=output_dir+'stokes_i/'+'pcamask',ipref=output_dir+'stokes_i/'+'total_intensity_43_part2',cat=output_dir+'stokes_qu/'+'mycat',
              debias=True,config = data_dir+'dimmconfig_pol2_compact_shortmap.lis')



