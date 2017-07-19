from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
from astropy.io import fits
import math
from sklearn import linear_model

x = np.linspace(0,1400, 1400)
data_poly = []
x_poly = []
# create data to be fitted
datafile = fits.open('summed_2.fits')
data_in = datafile[0].data
data = np.zeros(1400,dtype=np.float)
for j in range(1400):
    data[j] = data_in[j,0]
    
params = []    
for j in range(1400)[3:]:
    diff = data_in[j,0] - data_in[j-1,0]
    diff2 = data_in[j,0] - data_in[j-2,0]
    diff3 = data_in[j,0] - data_in[j-3,0]
    if (abs(diff)+abs(diff2)+abs(diff3))<30:
        data_poly.append(data_in[j,0])
        x_poly.append(j)
        params.append([1.])

z = np.polyfit(x_poly, data_poly, 13)
p = np.poly1d(z)
reg = linear_model.LinearRegression()
reg.fit(params, data_poly)
print reg.coef_
import pylab
y_poly = np.dot(np.array([reg.coef_]), np.array(params).transpose())
print y_poly.shape
pylab.plot(x_poly, y_poly[0], 'k')
pylab.plot(x_poly, data_poly, 'r+')
pylab.show()


hdu = fits.PrimaryHDU(data=data-p(x))
hdulist = fits.HDUList([hdu])
hdulist.writeto('summed_non_cont_2.fits', clobber=True)   