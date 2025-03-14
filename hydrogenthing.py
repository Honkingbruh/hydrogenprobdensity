import numpy as np
import scipy as spy
import math
from scipy import stats, special
from scipy.constants import find, physical_constants
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from mayavi import mlab
import time
import pyvista as pv

slice = True
start = time.time()
n = 3
l = 1
m = 0

a_0arr = np.array(physical_constants['Bohr radius'])
a_0 = (float)(a_0arr[0])

maxCoord = 2.8*a_0
if(n == 1):
    maxCoord = 2.8*a_0
if(n == 2):
    maxCoord = 19.2*a_0
if(n == 3):
    maxCoord = 32*a_0
if(n == 4):
    maxCoord = 100*a_0
maxPoints = 10000000
halfMaxPoints = maxPoints/2

hbar = spy.constants.hbar
mu, sigma = 0, 5*hbar
x = maxCoord*(np.random.random(maxPoints)-0.5)
y = maxCoord*(np.random.random(maxPoints)-0.5)
z = maxCoord*(np.random.random(maxPoints)-0.5)

print("1: " + str(time.time() - start))
#if(slice):
#    x = maxCoord*(np.random.random(halfMaxPoints)-0.5)
#    y = maxCoord*(np.random.random(halfMaxPoints)-0.5)
#    z = (maxCoord/2)*(np.random.random(halfMaxPoints))
    


#scaling = 2.8447564687401866*(10**31)*maxCoord
scaling = (10**31)*maxCoord


#density = stats.gaussian_kde(xyz)(xyz) 



#earr = np.array(physical_constants['elementary charge'])
#e = (float)(earr[0])
print(a_0)

r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arctan2(y, x)
phi = np.arccos(z/r)
print("2: " + str(time.time() - start))
#print("sharty ", spy.special.assoc_laguerre((2*r[0]/(n*a_0)), n - l - 1, 2*l + 1) - 10, " ajdoi ", (2*r[0]/(n*a_0)))
#np.sqrt(np.exp(x, 2) + np.exp(y, 2) + np.exp(z, 2))
#theta = np.array([])
#np.arccos(z/(np.sqrt(np.exp(x, 2) + np.exp(y, 2) + np.exp(z, 2))))
#phi = np.array([])
#np.sgn(y)*np.arccos(x/np.sqrt(np.exp(x, 2) + np.exp(y, 2)))
#for i in range(0, 1000):
    #r = np.append(r, (float)(np.sqrt((x[i])**2 + (y[i])**2 + (z[i])**2)))
    #theta = np.append(theta, (float)(np.arccos(z[i]/(np.sqrt((x[i])**2 + (y[i])**2 + (z[i])**2)))))
    #phi = np.append(phi, (float)(np.sign(y[i])*np.arccos(x[i]/np.sqrt((x[i])**2 + (y[i])**2))))
#print(spy.special.eval_genlaguerre(n - l - 1, 2*l + 1, x[0])*spy.special.sph_harm(m, n, theta[0], phi[0]))
#print(np.exp(spy.constants.e, -r[0]/(n*a_0[0]))*spy.special.eval_genlaguerre(n - l - 1, 2*l + 1, x[0])*spy.special.sph_harm(m, n, theta[0], phi[0]))
#print(2*(float)(a_0[0]))
#print(type(2*r[i]/(n*((float)(a_0[0])))))
rtp = np.vstack([r, theta, phi])
print("3: " + str(time.time() - start))
sqrtPart = math.sqrt(((2/n*a_0)**3)*math.factorial(n - l - 1)/(2*n*(math.factorial(n + l))))
print("4: " + str(time.time() - start))
def eChargePart(r):
    return np.exp(-r/(n*a_0))
def expTerm(r):
    return (2*r/(n*a_0))**l
def laguerreTerm(r):
    return spy.special.assoc_laguerre((2*r[0]/(n*a_0)), n - l - 1, 2*l + 1)
def harmPart(phi, theta):
    return spy.special.sph_harm(m, l, phi, theta)
#print("aiosdjw", r[i], " ", (-r[i]/(n*a_0)), " ", e**(-r[0]/(n*a_0)))
def psi(r, theta, phi):
    return sqrtPart*eChargePart(r)*expTerm(r)*laguerreTerm(r)*harmPart(phi, theta)

psiValues = psi(r, theta, phi)
print("5: " + str(time.time() - start))
probDensity = np.abs(psiValues)**2
#np.divide(probDensity, np.max(probDensity))
probDensity = (probDensity - np.min(probDensity)) / (np.max(probDensity) - np.min(probDensity))
print("6: " + str(time.time() - start))
threshold = 0.4
print(np.size(probDensity))
x, y, z = x[probDensity >= threshold], y[probDensity >= threshold], z[probDensity >= threshold]
probDensity = probDensity[probDensity >= threshold]
print(np.size(probDensity))
newMaxPoints = np.size(probDensity)
print("7: " + str(time.time() - start))

def getProbDensity(counter):
    return probDensity[counter]

#figure = mlab.figure('DensityPlot')
print(x[0], " ", y[0], " ", z[0], " ", probDensity[0])
#x = np.append(x, 0)
#y = np.append(y, 0)
#z = np.append(z, 0)
xyz = np.vstack([x,y,z])
pc = pv.PolyData(np.transpose(xyz))
print("8: " + str(time.time() - start))
#pointSizes = np.full(x.shape, 5.0)
#probDensity = np.append(probDensity, 0)
#pointSizes[pointSizes.size - 1] = 5.0 * 1835
#counter = 0
    
pc.plot(render_points_as_spheres=True, scalars=probDensity, point_size=5.0, cmap='magma')
print("9: " + str(time.time() - start))
pc = pv.PolyData(np.transpose(np.array([[0],[0],[0]])))
#pts = mlab.points3d(0.1*x/hbar, 0.1*y/hbar, 0.1*z/hbar, probDensity, scale_mode='none', scale_factor=scaling, colormap='magma')
#mlab.colorbar(pts, title='probability')
#mlab.axes()
print("program took ", time.time() - start, " s to run.")
#mlab.show()
