import numpy as np
from numpy import ndarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import pylab
from astropy.io import fits
import astropy.io.fits.hdu.hdulist
import os
namber = 11
hdus = fits.open('flat_slit'+repr(namber)+'_jos_g-1.fits',memmap=True)
master_flat = np.zeros((hdus[0].data.shape[0],hdus[0].data.shape[1]),dtype=np.float)
centr = hdus[0].data[0:512,0:512]
flux = np.sum(centr,0)
flux_centr = np.sum(flux,0)
hdus.close()

name = 'flat_slit'+repr(namber)+'_jos_g-'
for i in range(5):
	hdu = fits.open(name+repr(i+1)+'.fits',memmap=True)
	flat = hdu[0].data
	centr_flat = hdu[0].data[0:512,0:512]
	flux_flat = np.sum(centr_flat,0)
	flux_centr_flat = np.sum(flux_flat,0)
	k = flux_centr/flux_centr_flat
	master_flat = master_flat+flat*k
	hdu.close()

try:
	os.mkdir("done_res")
except OSError:
	pass

hdu = fits.PrimaryHDU()
hdu.data = master_flat/5
hdu.writeto('done_res/master_flat-'+repr(namber)+'.fits', clobber=True)

