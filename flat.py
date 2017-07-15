import numpy as np
from numpy import ndarray
import math
import pylab
from astropy.io import fits
import astropy.io.fits.hdu.hdulist
import os

print('number slit:')
number = int(input())
print('with normalization - 1, without normalization - 0:')
norm = int(input()) 
hdus = fits.open('flat_slit'+repr(number)+'_jos_g-1.fits',memmap=True)
master_flat = np.zeros((hdus[0].data.shape[0],hdus[0].data.shape[1]),dtype=np.float)
n = hdus[0].header['NAXIS1']
centr = hdus[0].data[n/4:3*n/4,n/4:3*n/4]
flux = np.sum(centr,0)
flux_centr = np.sum(flux,0)
hdus.close()

name = 'flat_slit'+repr(number)+'_jos_g-'
print('using fits:')
for i in range(5):
	hdu = fits.open(name+repr(i+1)+'.fits',memmap=True)
	print(name+repr(i+1)+'.fits')
	flat = hdu[0].data
	centr_flat = hdu[0].data[n/4:3*n/4,n/4:3*n/4]
	flux_flat = np.sum(centr_flat,0)
	flux_centr_flat = np.sum(flux_flat,0)
	if (norm==1):
		k = flux_centr/flux_centr_flat
	else:
		k=1
	master_flat = master_flat+flat*k
	hdu.close()
print('results:')
stat = np.std(master_flat[100:200,100:200]/5)
mean = np.mean(master_flat[100:200,100:200]/5)
print('stat = '+repr(stat)+'	mean = '+ repr(mean))

try:
	os.mkdir("done_res")
except OSError:
	pass

hdu = fits.PrimaryHDU()
hdu.data = master_flat/5
hdu.writeto('done_res/master_flat-'+repr(number)+'.fits', clobber=True)

