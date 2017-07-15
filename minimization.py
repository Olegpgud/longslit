from scipy.optimize import minimize, rosen, rosen_der
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
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

datafile = fits.open('summed_non_cont_2.fits')
data_in = datafile[0].data
data = np.zeros(1400,dtype=np.float)
mask = np.zeros(1400,dtype=np.float)
for j in range(1400):
    data[j] = data_in[j]
    if (j>300) & (j<1000):
        mask[j]=1.

def x_to_l(x):
    return 11205.9595+2.47028877*x+0.00028425*x**2-7.8345e-08*x**3
    
Q = 9714.89573
def profile_function(lsf, l_c,flux):
    model = np.zeros(1400,dtype=np.float)
    for i in range(model.shape[0]):
        arg = int((1./2.47)*(x_to_l(i) - l_c) + len(lsf)/2.)
        if (arg>0) & (arg<len(lsf)):
            model[i] = flux * lsf[int(arg)]
    return model





model = np.zeros(1400,dtype=np.float)
def func2min(a):
    model = np.zeros(1400,dtype=np.float)
    for line in linelist:
        model = model + Q*profile_function(a,line[0],line[1])
    return sum((model-data)**2/data) + 0*(sum(a**2)+sum(np.diff(a)**2))


x0 = np.ones(40,dtype=np.float)

x = np.linspace(0,20, 20)
x0 = np.exp(-np.power(x - len(x)/2., 2.) / (2 * np.power(2.0, 2.)))

'''
print func2min(x0)
import pylab
pylab.plot(range(1400), func2min(x0), 'r')
#pylab.plot(x, fcn2min(params, x,data), 'r')
pylab.show()
'''


model = np.zeros(1400,dtype=np.float)
for line in linelist:
    model = model + Q*profile_function(x0,line[0],line[1])

'''
pylab.plot(range(len(model)), model, 'r')
pylab.plot(range(len(model)), data, 'k')
pylab.show()
'''
res = minimize(func2min, x0, method='Nelder-Mead', tol=1e-6)
print res.x

pylab.plot(range(len(res.x)), res.x, 'k')
pylab.show()
#pylab.plot(range(len(res.x)), res.x, 'r')




'''
[-5.76525358e-06  -2.65140479e-06   5.33515941e-04   9.07608756e-03
   3.69647727e-02   5.01099431e-02   1.87579224e-01   5.43860736e-01
   8.52775570e-01   9.82225679e-01   7.95348852e-01   4.98583669e-01
   1.87909277e-01   7.67543622e-02   3.71614220e-02   1.82079836e-02
   1.02745861e-02   5.28429251e-04  -1.16494738e-04   8.23001517e-06]



'''