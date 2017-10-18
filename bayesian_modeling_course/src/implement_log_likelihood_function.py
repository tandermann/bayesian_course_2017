#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 11:47:48 2017

@author: tobias
"""

import numpy as np
from numpy import *
from scipy.special import gamma
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from matplotlib import cm

#%matplotlib qt

s = 5.
x = 4.
u = 7.
pi
e
# manually log-transforming a function:
y = (1/(sqrt(2*pi*s**2)))*e**(-((x-u)**2)/(2*s**2))
y

# now let's log transform the function by hand
log_transformed_y = -0.5*log(2*pi*s**2)-((x-u)**2)/(2*s**2)
log_transformed_y
log_transformed_y == log(y)


# probability density functions
def pdf_normal(x,u,s):
    lik = (1/(sqrt(2*pi*s**2)))*e**(-((x-u)**2)/(2*s**2))
    return lik

def log_pdf_normal(x,u,s):
    lik = -0.5*log(2*pi*s**2)-((x-u)**2)/(2*s**2)
    return lik

def pdf_cauchy(x,x0,s):
    lik = 1/(pi*s*(1+((x-x0)/s)**2))
    return lik

def pdf_exp(x,l):
    lik = l*e**(-l*x)
    return lik

def pdf_gamma(x,a,b):
    lik = (b**a/gamma(a)) * x**(a-1) * e**(-b*x)
    return lik

def log_pdf_gamma(x,a,b):
    lik = a*log(b)-log(gamma(a))+a*log(x)-log(x)-b*x
    return lik



#______________________________DATA_________________________________________
# the exercise with panda body mass data
# here is the data as a numpy array:
pandas_bm =np.array([84.74, 84.88, 94.60, 96.37, 102.93, 109.11, 125.76])
u = 100 #mean
s = 5 #stdev
# what's the likelihood of our panda data under these parameters?
sum(log_pdf_normal(pandas_bm,u,s))



#____________________PLOT CAUCHY DISTRIBUTION____________________
# pdf_cauchy
x = np.arange(min(pandas_bm)-10,max(pandas_bm)+10,0.01)
plt.plot(x,pdf_cauchy(x,100,s))
plt.plot(pandas_bm,pdf_cauchy(pandas_bm,100,s),'ro')
#print('Global likelihood with chosen parameters:',sum(pdf_cauchy(pandas_bm,100,s)))



#_____________________________NORMAL____________________________________________
#____________________PLOT NORMAL DISTRIBUTION____________________
# plot the panda data under normal distribution
x = np.arange(min(pandas_bm)-10,max(pandas_bm)+10,0.01)
plt.plot(x,pdf_normal(x,u,s))
plt.plot(pandas_bm,pdf_normal(pandas_bm,u,s),'ro')
#print('Global likelihood with chosen parameters:',sum(pdf_normal(pandas_bm,u,s)))
#____________________PLOT LOGTRANSFORMED NORMAL DISTRIBUTION____________________
# plot under log transformed normal distribution
plt.plot(x,log_pdf_normal(x,u,s))
plt.plot(pandas_bm,log_pdf_normal(pandas_bm,u,s),'ro')
print('Global likelihood with chosen parameters:',sum(log_pdf_normal(pandas_bm,u,s)))
#____________________MAKE MATRIX (LOGTRANSFORMED NORMAL)____________________
# make a matrix with likelihood values for a range of u and s values
n_steps=100
x_vec = np.linspace(84,125,n_steps)
y_vec = np.linspace(6,50,n_steps)
z_2d  = np.zeros((n_steps,n_steps))
for i in range(0,n_steps):
    for j in range(0,n_steps):
        x= x_vec[i]
        y= y_vec[j]
        z_2d[j,i] = sum(log_pdf_normal(pandas_bm,x,y))
# get the best point from the matrix (best combo of u and s)
i,j = np.unravel_index(z_2d.argmax(), z_2d.shape)
best_point = [x_vec[j],y_vec[i]]
print ("max likelihood:",z_2d[i,j], "mu:", x_vec[j], "sig:",y_vec[i])
#____________________PLOT LIKELIHOOD SURFACE (LOGTRANSFORMED NORMAL)___________
# plot contour plot
f=plt.figure()
levs = np.linspace(z_2d.min(),z_2d.max(),200)      
CS = plt.contour(x_vec, y_vec, z_2d,levs)
plt.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
#plt.axis([min(pandas_bm), max(pandas_bm), 0, 30])
#__________insert: make up some data__________
x_test = np.random.normal(100,15,[500])
y_test = np.random.normal(15,3,[500])
plt.plot(x_test,y_test,'bo')
#_________________end insert__________________
plt.plot(best_point[0], best_point[1], 'rX',label="maximum likelihood")
plt.annotate('maximum likelihood', xy=(best_point[0], best_point[1]), xytext=(best_point[0] + 10, best_point[1] + 10), arrowprops=dict(facecolor='red',width=1.5, headwidth=8))
plt.show()




#_____________________________GAMMA____________________________________________
#____________________PLOT GAMMA DISTRIBUTION____________________
# gamma distribution
a = 58.
b = .58
plt.plot(x,pdf_gamma(x,a,b))
plt.plot(pandas_bm,pdf_gamma(pandas_bm,a,b),'ro')
#____________________PLOT LOGTRANSFORMED GAMMA DISTRIBUTION____________________
# log transformed gamma distirbution
a = 58.
b = .58
x = np.arange(min(pandas_bm)-10,max(pandas_bm)+10,0.01)
# check if function is correct
#log(pdf_gamma(x,a,b)) == log_pdf_gamma(x,a,b)
plt.plot(x,log_pdf_gamma(x,a,b))
plt.plot(pandas_bm,log_pdf_gamma(pandas_bm,a,b),'ro')
#____________________MAKE MATRIX (LOGTRANSFORMED GAMMA)____________________
n_steps=100
x_vec = np.linspace(50.,100,n_steps)
y_vec = np.linspace(0.05,2,n_steps)
z_2d  = np.zeros((n_steps,n_steps))
for i in range(0,n_steps):
    for j in range(0,n_steps):
        x= x_vec[i]
        y= y_vec[j]
        z_2d[j,i] = sum(log_pdf_gamma(pandas_bm,x,y))
i,j = np.unravel_index(z_2d.argmax(), z_2d.shape)
best_point = [x_vec[j],y_vec[i]]
print ("max likelihood:",z_2d[i,j], "a:", x_vec[j], "b:",y_vec[i])
#____________________PLOT LIKELIHOOD SURFACE (LOGTRANSFORMED GAMMA)___________
# plot contour plot
f=plt.figure()
levs = np.linspace(z_2d.min(),z_2d.max(),200)      
CS = plt.contour(x_vec, y_vec, z_2d,levs)
plt.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
#plt.axis([min(pandas_bm), max(pandas_bm), 0, 30])
#plt.plot(best_point[0], best_point[1], 'rX',label="maximum likelihood")
#plt.annotate('maximum likelihood', xy=(best_point[0], best_point[1]), xytext=(best_point[0] + 10, best_point[1] + 10), arrowprops=dict(facecolor='red',width=1.5, headwidth=8))
plt.show()




# plot interactive 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.view_init(azim=0)
surf = ax.plot_surface(x_vec, y_vec, z_2d , rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()


