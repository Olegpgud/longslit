from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
from astropy.io import fits

x = np.linspace(0,1400, 1400)

def profile_function(l_c, sigma,flux,x):
    return flux * np.exp(-np.power(x - l_c, 2.) / (2 * np.power(sigma, 2.)))
    
def get_lambd(x,a,b,c):
    return a + b*x +c*x**2
    
    
linelist=[
    [5051,1],
    [5370,0.8],
    [5590,0.6],
    [6240,0.7]
    ]

# create data to be fitted
datafile = fits.open('summed.fits')
data_in = datafile[0].data
data = np.zeros(1400,dtype=np.float)
for j in range(1400):
    data[j] = data_in[j,0]
'''
data = (5. * np.sin(2 * x - 0.1) * np.exp(-x*x*0.025) +
        np.random.normal(size=len(x), scale=0.2) )
'''
# define objective function: returns the array to be minimized
model_ = np.zeros(1400,dtype=np.float)
def fcn2min(params, x,data):
    """ model decaying sine wave, subtract data"""
    a = params['a']
    b = params['b']
    c = params['c']
    Q = params['Q']
    sigma = params['sigma']
    model = np.zeros(1400,dtype=np.float)
    for line in linelist:
        model = model + Q*profile_function(line[0],sigma,line[1],get_lambd(x,a,b,c))
        #print get_lambd(x,a,b,c)
    model_ = model
    return (model-data)

# create a set of Parameters
params = Parameters()
params.add('a',   value= 9285)
params.add('b', value= -4.5)
params.add('c', value= 0.0)
params.add('Q', value= 1700.0)
params.add('sigma', value= 7.0)
#print get_lambd(x,7000,-0.3,0)


# do fit, here with leastsq model
minner = Minimizer(fcn2min, params, fcn_args=(x, data))
result = minner.minimize( method='nelder')

# calculate final result
final = data + result.residual

# write error report
report_fit(result)

# try to plot results



try:
    import pylab
    pylab.plot(x, data, 'k+')
    pylab.plot(x, final, 'r')
    #pylab.plot(x, fcn2min(params, x,data), 'r')
    pylab.show()
except:
    pass
