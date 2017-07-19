from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
from astropy.io import fits

x = np.linspace(0,1400, 1400)

def H (x,n):
	if n==3:
		return(-12*x+8*x**3)
	else:
		return(12-48*x**2+16*x**4)

def profile_function(l_c, sigma,flux,x,h3,h4):
    return flux * np.exp(-np.power(x - l_c, 2.) / (2 * np.power(sigma, 2.)))*(1+H((x-l_c)/(sigma*(2**0.5)),3)*h3+h4*H((x-l_c)/(sigma*(2**0.5)),4))
    
def get_lambd(x,a,b,c,d,e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4
    
'''    
linelist=[
    [5051,1],
    [5370,0.8],
    [5590,0.6],
    [6240,0.7]
    ]
'''
linelist=[
    [14254.351,0.017],
    [14097.5,0.13],
    [13829.51,0.04],
    [13914.4,0.05],
    [13722.34,1],
    [13682.29,0.20],
    [13626.38,0.3],
    [13603.05,0.05],
    [13507.88,0.87],
    [13410.26,0.05],
    [13370.77,0.6],
    [13316.85,0.3],
    [13276.27,0.4],
    [13231.72,0.23],
    [13217.16,0.27],
    [13012.0,0.15],
    [12960.1,0.65],
    [12936.9,0.08],
    [12806.2,0.2],
    [12705.9,0.2],
    [12491.08,0.4],
    [12459.53,0.17],
    [12442.73,0.58],
    [12406.22,0.36],
    [12359.67,0.02],
    [12346.77,0.09],
    [12143.06,0.10],
    [12115.64,0.15]
    ]

# create data to be fitted
datafile = fits.open('summed_non_cont_2.fits')
data_in = datafile[0].data
data = np.zeros(1400,dtype=np.float)
mask = np.zeros(1400,dtype=np.float)
for j in range(1400):
    data[j] = data_in[j]
    if (j>300) & (j<1000):
        mask[j]=1.
'''
data = (5. * np.sin(2 * x - 0.1) * np.exp(-x*x*0.025) +
        np.random.normal(size=len(x), scale=0.2) )
'''
# define objective function: returns the array to be minimized
model = np.zeros(1400,dtype=np.float)
def fcn2min(params, x,data):
    """ model decaying sine wave, subtract data"""
    a = params['a']
    b = params['b']
    c = params['c']
    d = params['d']
    e = params['e']
    Q = params['Q']
    sigma = params['sigma']
    h3 = params['h3']
    h4 = params['h4']
    model = np.zeros(1400,dtype=np.float)
    for line in linelist:
        model = model + Q*profile_function(line[0],sigma,line[1],get_lambd(x,a,b,c,d,e),h3,h4)
        #print get_lambd(x,a,b,c)
    
    return mask*(model-data)
# create a set of Parameters
params = Parameters()
params.add('a',   value= 11005.3960, vary=True)
params.add('b', value= 2.6, vary=True)
params.add('c', value= 0., vary=True)
params.add('d', value= 0., vary=True)
params.add('e', value= 0., vary=False)
params.add('Q', value= 10000.0, vary=True)
params.add('sigma', value= 4.0, vary=True, min=3.)
params.add('h3',value= 0,vary=False, min=-0.0005, max=0.0005)
params.add('h4',value= 0,vary=False, min=-0.0005, max=0.0005)
#print get_lambd(x,7000,-0.3,0)


# do fit, here with leastsq model
minner = Minimizer(fcn2min, params, fcn_args=(x, data))
result = minner.minimize()

# calculate final result
final = data + result.residual

# write error report
report_fit(result)

# try to plot results

model = np.zeros(1400,dtype=np.float)
for line in linelist:
    model = model + result.params['Q'].value*profile_function(line[0],result.params['sigma'].value,line[1],get_lambd(x,result.params['a'].value,result.params['b'].value,result.params['c'].value,result.params['d'].value,result.params['e'].value),0,0)


try:
    import pylab
    pylab.plot(x, data, 'k')
    pylab.plot(x, model, 'r')
    pylab.plot(x, model-data, 'g')
    #pylab.plot(x, fcn2min(params, x,data), 'r')
    pylab.show()
except:
    pass

