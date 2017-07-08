from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
from astropy.io import fits

x = np.linspace(0,1400, 1400)
data_poly = []
x_poly = []
# create data to be fitted
datafile = fits.open('summed.fits')
data_in = datafile[0].data
data = np.zeros(1400,dtype=np.float)
for j in range(1400):
    data[j] = data_in[j,0]
    
    
for j in range(1400)[2:]:
    diff = data_in[j,0] - data_in[j-1,0]
    diff2 = data_in[j,0] - data_in[j-2,0]
    if (abs(diff)+abs(diff2))<40:
        data_poly.append(data_in[j,0])
        x_poly.append(j)
    

z = np.polyfit(x_poly, data_poly, 14)
p = np.poly1d(z)



try:
    import pylab
    pylab.plot(x, data-p(x), 'k')
    pylab.plot(x, p(x), 'r')
    pylab.show()
except:
    pass
