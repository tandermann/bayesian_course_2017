#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 13:40:10 2017

@author: tobias
"""

import csv
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

#___________________________ LIKELIHOOD FUNCTION ______________________________
def log_pdf_normal(x,m,s):
    lik = -0.5*np.log(2*np.pi*s**2)-((x-m)**2)/(2*s**2)
    return lik

#__________________________________ PRIORS ____________________________________
# std parameter
def log_pdf_gamma(x,a,b):
    lik = a*np.log(b)-np.log(gamma(a))+a*np.log(x)-np.log(x)-b*x
    return lik

def pdf_gamma(x,a,b):
    lik = (b**a/gamma(a)) * x**(a-1) * np.e**(-b*x)
    return lik

# mu parameter
def log_pdf_uniform(x,a,b):
    if x<a or x>b:
        lik = -np.inf
    else:
        lik = -np.log(b-a)
    return lik

# std parameter
def pdf_normal(x,u,s):
    lik = (1/(np.sqrt(2*np.pi*s**2)))*np.e**(-((x-u)**2)/(2*s**2))
    return lik

# _____________________________ PROPOSAL FUNCTIONS ____________________________ 
# draw a random float from a normal distribution and store it as the new x value
def update_normal(x,d=0.25): # d is the standard deviation
    x_prime = np.random.normal(x,d)
    return x_prime

def sliding_window(x,d=2): #d is the size of the normal window
    x_prime = x+(np.random.random()-.5)*d
    return x_prime

# _____________________________________ DATA __________________________________
pandas_bm =np.array([84.74, 84.88, 94.60, 96.37, 102.93, 109.11, 125.76])
commute_length = np.array([63.5,81.,61.5,120.5,64.5,66.3,62.3,65.6,66.5,70.5,62.3,66.5,68.9])
data = pandas_bm

#___________________________________MCMC SETTINGS______________________________
# initiate log file
logfile = open("/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples.txt","w")
wlog = csv.writer(logfile,delimiter='\t')
wlog.writerow(["it","post","prior","lik","mu","sig"])
logfile.flush() #this will tell python to write to the file whenever it is being stated in the code

log_every_n_lines = 100

# initial values for the MCMC
mu_init = 100
sd_init = 10

# parameters defining prior distribution
# prior probability of mu will be calculated from a normal distribution between a and b
mu_prior_a_unif = 1
mu_prior_b_unif = 1000
# prior probability of sd is calculated from a gamma distribution with the parameters a and b
sd_prior_shape_gamma = 1
sd_prior_rate_gamma = .1

#MCMC settings
n_iterations = 100000
# 'window size' for new proposals
d_mu = 2
d_sd = 1

# calculating the initial likelihood and prior probability
current_lik = sum(log_pdf_normal(data,mu_init,sd_init))
current_prior_mu = log_pdf_uniform(mu_init,mu_prior_a_unif,mu_prior_b_unif)
current_prior_sd = pdf_gamma(sd_init,sd_prior_shape_gamma,sd_prior_rate_gamma)
current_prior = current_prior_mu + current_prior_sd
current_posterior = current_lik + current_prior

# count how many steps we accepted
acceptance_counter = 0

# this ratio defines how often parameter 1 (mu) is sampled in relation to parameter 2 (sd)
update_frequency = .5

for i in range(0,n_iterations):
    # proposing new parameter values  
    if np.random.random() < update_frequency:
        proposed_mu = update_normal(mu_init,d_mu)
        proposed_sd = sd_init
    else:
        proposed_mu = mu_init
        proposed_sd = update_normal(sd_init,d_sd) 
    #proposed_mu = update_normal(mu_init,d_mu)
    #proposed_sd = update_normal(sd_init,d_sd)
    # calculate likelihood ratio
    new_likelihood = sum(log_pdf_normal(data,proposed_mu,proposed_sd))
    #likelihood_ratio = new_likelihood/current_lik
    
    # calculate prior ratio
    new_prior_mu = log_pdf_uniform(proposed_mu,mu_prior_a_unif,mu_prior_b_unif)
    new_prior_sd = pdf_gamma(proposed_sd,sd_prior_shape_gamma,sd_prior_rate_gamma)
    new_prior = new_prior_mu + new_prior_sd
    #prior_ratio = new_prior/current_prior_sd
    
    # calculate new posterior
    new_posterior = new_likelihood + new_prior
    posterior_ratio = new_posterior - current_posterior
    
    # draw a random number between 0 and 1
    random_draw = np.random.random()
    transformed_draw = np.log(random_draw)
    # see if that random draw lies within your acceptance ratio
    if transformed_draw < posterior_ratio:
        accept = True
        mu_init = proposed_mu
        sd_init = proposed_sd
        current_lik = new_likelihood
        current_prior = new_prior
        current_posterior = new_likelihood+new_prior
        acceptance_counter += 1
    else:
        accept = False

    # if the iteration can be divided by the user-input logging frequency, print to file
    if i % log_every_n_lines == 0:
        wlog.writerow([i,current_posterior,current_prior,current_lik,mu_init,sd_init])
        logfile.flush() #this will tell python to write to the file whenever it is being stated in the code


print("Acceptance rate is", float(acceptance_counter)/float(n_iterations))   



















#___________________________MESSING AROUND__________________________

def normal_dist(x,m,s):
    y = (1/(np.sqrt(2*np.pi*s**2)))*np.e**(-((x-m)**2)/(2*s**2))
    return y

import pandas as pd
mcmc_log = pd.read_csv('/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples.txt', sep = '\t')

plt.plot(mcmc_log.iloc[:,0],mcmc_log.iloc[:,1]), plt.title('post'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[:,0],mcmc_log.iloc[:,2]), plt.title('prior'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[:,0],mcmc_log.iloc[:,3]), plt.title('lik'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[:,0],mcmc_log.iloc[:,4]), plt.title('mu'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[:,0],mcmc_log.iloc[:,5]), plt.title('sig'),plt.show()#, plt.axis([0, 10000,90,110])


#_______________________________________DENSITY PLOT___________________________
import seaborn as sns

data = mcmc_log.mu
sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw=1,color="green", shade=True)
plt.show()

data = mcmc_log.sig
sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw=1,color="blue", shade=True)
plt.show()

## plto a normal distribution
#plt.hist(mcmc_log.mu,bins=100)
#plt.plot(density(mcmc_log.mu))
#x_test = np.arange(50,95,1)
#y_test = normal_dist(x_test,70.883,6.371)
#plt.plot(x_test,y_test)
#

