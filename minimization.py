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
a = 11205.9595
b = 2.47028877
c = 0.00028425
d = -7.8345e-08
def x_to_l(x):
    return a+b*x+c*x**2+d*x**3
    
Q = 9714.89573
def profile_function(lsf, l_c,flux):
    model = np.zeros(1400,dtype=np.float)
    for i in range(model.shape[0]):
        x_ = x_to_l(i)
        step = b + 2*c*i+3*d*i**2
        arg = int((2./step)*(x_ - l_c) + len(lsf)/2.)
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

x = np.linspace(0,40, 40)
x0 = np.exp(-np.power(x - len(x)/2., 2.) / (2 * np.power(4.5, 2.)))

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


pylab.step(range(len(model)), model, 'r')
pylab.step(range(len(model)), data, 'k')
pylab.show()
''''''
res = minimize(func2min, x0, method='Nelder-Mead', tol=1e-6)
print res.x

pylab.step(range(len(res.x)), res.x, 'k')
pylab.show()
#pylab.plot(range(len(res.x)), res.x, 'r')




'''
[-5.76525358e-06  -2.65140479e-06   5.33515941e-04   9.07608756e-03
   3.69647727e-02   5.01099431e-02   1.87579224e-01   5.43860736e-01
   8.52775570e-01   9.82225679e-01   7.95348852e-01   4.98583669e-01
   1.87909277e-01   7.67543622e-02   3.71614220e-02   1.82079836e-02
   1.02745861e-02   5.28429251e-04  -1.16494738e-04   8.23001517e-06]

[ -6.67369321e-21  -2.19018261e-19   3.02186212e-16   1.75483163e-14
  -1.52359484e-12  -1.27493080e-11  -6.32976562e-10  -2.25594616e-08
   8.00759817e-08  -7.35494381e-06   1.56417795e-04   9.23454621e-04
   1.30453454e-02   2.06351948e-01   3.82127948e-01   5.87146918e-01
   7.55316280e-01   9.26857928e-01   9.94030663e-01   9.66210850e-01
   9.07713531e-01   6.98211637e-01   5.89245352e-01   3.86324002e-01
   2.33214267e-01   1.27986462e-01   1.03583912e-01   5.63316968e-03
   1.51066070e-04  -4.15629312e-05   3.04638534e-06  -2.76986592e-07
  -2.53194207e-08  -2.16183630e-09  -3.57726113e-11  -3.64990253e-13
  -3.61063999e-15  -7.84492379e-17   1.57516457e-18  -1.18038669e-20]



[ -1.74296856e-21   2.63493615e-19   6.18838083e-18   4.21812986e-15
  -1.52443273e-13   3.66697345e-12  -7.31576970e-10  -1.56480923e-08
  -6.13783532e-08  -9.24880972e-06   3.72089877e-04  -1.46355670e-03
   1.51887328e-02   1.43176736e-01   2.97138563e-01   4.35013336e-01
   7.55814337e-01   8.14643391e-01   1.02428515e+00   9.45312244e-01
   9.00067166e-01   6.68346249e-01   5.76047205e-01   3.24695739e-01
   2.05317912e-01   7.09126866e-02   5.99446198e-02  -9.21096032e-03
  -3.77839789e-03  -8.98754569e-07  -4.67631769e-06   1.00748197e-07
  -4.39610537e-08  -3.70041512e-10   1.81881117e-11   1.42956849e-13
   1.43531805e-15   4.97593577e-17   1.52830539e-19   1.75654304e-22]
'''