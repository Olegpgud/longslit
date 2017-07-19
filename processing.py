import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import matplotlib.pyplot as plt
from pylab import * # This program requires the matplotlib python library for the plots.
from scipy import ndimage
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.interpolate import griddata
import pylab
font = {'family': 'Droid Sans',
        'weight': 'normal',
        'size': 20}


cmap = plt.get_cmap('jet_r')
calib = fits.open('slit7_jos_g_4s-1.fits')
calib_dat = calib[0].data[0:4000,530:1532]
flat = fits.open('slit7_jos_g_4s-1.1499362072.fits')
flat_dat = flat[0].data[0:4000,530:1532]
flat_dat[flat_dat==0] = -1
devided = np.array(calib_dat) / np.array(flat_dat)

hdu = fits.PrimaryHDU(data=devided)
hdulist = fits.HDUList([hdu])
hdulist.writeto('devided.fits', clobber=True)

regions = [[200,250],[330,370],[480,530],[600,650],[700,770],[830,900],[917,970],[1000,1120],[1117,1230],[1330,1370],[1455,1490],[1560,1660],[1810,1850]]
q = 0
poly_arr = []
poly_arr_x = []
for region in regions:
    q = q + 1
    x_corr = []
    y_corr = []
    y_min = region[0]
    y_max = region[1]
    strip0 = devided[y_min:y_max,500].copy()
    for i in range(1000)[10:]:
        strip = devided[y_min:y_max,i].copy()
        corr = np.correlate(strip0, strip, mode='same')
        ind = corr.tolist().index(max(corr)) - (y_max-y_min)/2
        z = np.polyfit([ind-1,ind,ind+1],[corr[ind-1+(y_max-y_min)/2],corr[ind+(y_max-y_min)/2],corr[ind+1+(y_max-y_min)/2]],2)
        x_max = -0.5*z[1]/z[0]
        diff = 0
        if i>13:
            diff = abs(x_max - y_corr[-1])+abs(x_max - y_corr[-2])
        print x_max, ind
        if not (diff>20 and i>10):
            x_corr.append(i)
            y_corr.append(x_max)
    z = np.polyfit(np.array(x_corr), np.array(y_corr), 4)

    p = np.poly1d(z)
    poly_arr.append(z)
    poly_arr_x.append((y_max+y_min)/2)
    final_shifted=devided[y_min:y_max,:].copy()
    for i in range(1000):
        strip = devided[y_min:y_max,i].copy()
        output_arr = ndimage.interpolation.shift(np.array(strip), p(i))
        final_shifted[:,i] = output_arr
        #hdu = fits.PrimaryHDU(data=final_shifted)
        #hdulist = fits.HDUList([hdu])
        #hdulist.writeto('shifted-'+str(q)+'.fits', clobber=True)
    print z
    color = cmap(float(q)/20)
    plot(x_corr,y_corr-p(x_corr)+q*0.3,c=color)
plt.show()
poly_arr = np.array(poly_arr)

f0 = interp1d(poly_arr_x, poly_arr[:,0], kind='cubic')
f1 = interp1d(poly_arr_x, poly_arr[:,1], kind='cubic')
f2 = interp1d(poly_arr_x, poly_arr[:,2], kind='cubic')

y_arr = np.array(range(1800))[250:]
f0_coeffs = f0(y_arr)
f1_coeffs = f1(y_arr)
f2_coeffs = f2(y_arr)

shifts_arr = np.zeros((1800,1000), dtype=np.float)

for i in range(len(y_arr)):
    p_interpoled = np.poly1d([f0_coeffs[i],f1_coeffs[i],f2_coeffs[i]])
    diff_y = np.polyval(p_interpoled, range(1000))
    shifts_arr[y_arr[i]] = diff_y


print devided.shape

points_x = []
points_y = []
values = []
for i in range(1800):
    for j in range(1000):
        points_y.append(i+shifts_arr[i,j]+0.)
        points_x.append(j)
        values.append(devided[i,j])
'''
f = interpolate.interp2d(np.arange(0, 1000, 1), np.arange(0,1800 , 1), np.array(devided)[0:1800,0:1000], kind='linear')
x_new = np.arange(0, 1000, 1)
y_new = np.arange(0, 1800, 1)

final_res = f(x_new, y_new)
'''
grid_z2 = griddata((np.array(points_x),np.array(points_y)), np.array(values), (np.linspace(0,1000,1000)[None,:], np.linspace(0,1800,1800)[:,None]), method='cubic')
#hdu = fits.PrimaryHDU(data=grid_z2)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('shifts.fits', clobber=True)

cutted = grid_z2[270:1700,:]
print cutted
summed_result = np.zeros((cutted.shape[0],1000),dtype=np.float)
for i in range(cutted.shape[0]):
    summed_result[i,0:1000] = np.nansum(cutted[i])
hdu = fits.PrimaryHDU(data=summed_result)
hdulist = fits.HDUList([hdu])
hdulist.writeto('summed.fits', clobber=True)   