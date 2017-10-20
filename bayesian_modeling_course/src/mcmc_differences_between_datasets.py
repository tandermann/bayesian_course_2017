#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:56:47 2017

@author: tobias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:15:35 2017

@author: tobias
"""

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
def log_pdf_gamma(x,a,b):
    lik = a*np.log(b)-np.log(gamma(a))+a*np.log(x)-np.log(x)-b*x
    return lik

#__________________________________ PRIORS ____________________________________
# shape and rate parameter
def pdf_gamma(x,a,b):
    lik = (b**a/gamma(a)) * x**(a-1) * np.e**(-b*x)
    return lik

# _____________________________ PROPOSAL FUNCTIONS ____________________________ 
# draw a random float from a normal distribution and store it as the new x value
def update_normal(x,d=0.25): # d is the standard deviation
    x_prime = np.random.normal(x,d)
    return x_prime

def sliding_window(x,d=2): #d is the size of the normal window
    x_prime = x+(np.random.random()-.5)*d
    return x_prime

# this updating function is useful because it updates the new value in percent rather than in absolute values. it is however an assymetric function, which is why we need to consider the Hastings ratio, whcih is exported by this function and needs to be added to the posterior ratio
def update_multiplier_proposal(i,d=1.2):
    u = np.random.uniform(0,1,np.shape(i))
    l = 2*np.log(d)
    m = np.exp(l*(u-.5))
    ii = i * m
    U=sum([np.log(m)])
    return ii, U

# _____________________________________ DATA __________________________________
pandas_bm = np.array([84.74, 84.88, 94.60, 96.37, 102.93, 109.11, 125.76])
brown_bears = np.array([105.12,121.15,130.23,140.33,147.3,150.7])
joined_data = np.append(pandas_bm,brown_bears)
# create an array with species indeces
species_ind = np.append(np.repeat(0,len(pandas_bm)) , np.repeat(1,len(brown_bears)))
data = joined_data

#___________________________________MCMC SETTINGS______________________________
# initiate log file
logfile = open("/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples_comparing_means_gamma.txt","w")
wlog = csv.writer(logfile,delimiter='\t')
wlog.writerow(["it","joined_post","joined_prior","joined_lik","shape_sp0","rate_sp0","shape_sp1","rate_sp1","mean_bm_sp0","mean_bm_sp1","sp1-sp0"])
logfile.flush() #this will tell python to write to the file whenever it is being stated in the code

log_every_n_lines = 10

# initial values for the MCMC
shape_init = 5
rate_init = 5

# parameters defining prior distribution
# prior probability of shape will be calculated from a normal distribution between a and b
shape_prior_mu = 1
shape_prior_sd = .1
# prior probability of rate is calculated from a gamma distribution with the parameters a and b
rate_prior_mu = 1
rate_prior_sd = .1

#MCMC settings
n_iterations = 100000
# 'window size' for new proposals
d_shape = 1
d_rate = .1
update_factor = 1.5

# define shape and rate for both species
shape_init_sp0 = shape_init
shape_init_sp1 = shape_init
rate_init_sp0 = rate_init
rate_init_sp1 = rate_init

# calculating the initial likelihood and prior probability
current_lik_sp0 = sum(log_pdf_gamma(data[species_ind ==0],shape_init_sp0,rate_init_sp0))
current_lik_sp1 = sum(log_pdf_gamma(data[species_ind ==1],shape_init_sp1,rate_init_sp1))
current_lik_joined = current_lik_sp0 + current_lik_sp1

# calculating the initial prior probability of the parameters
current_prior_shape_sp0 = log_pdf_gamma(shape_init,shape_prior_mu,shape_prior_sd)
current_prior_rate_sp0 = log_pdf_gamma(rate_init,rate_prior_mu,rate_prior_sd)
current_prior_sp0 = current_prior_shape_sp0 + current_prior_rate_sp0
current_prior_sp1 = current_prior_sp0

# calculating the initial posterior of the data
current_posterior_joined = current_lik_sp0 + current_lik_sp1 + current_prior_sp0 + current_prior_sp1

# count how many steps we accepted
acceptance_counter = 0

# this ratio defines how often parameter 1 (shape) is sampled in relation to parameter 2 (rate)
update_frequency = .5



for i in range(0,n_iterations):
    # initiate a hastings ratio of 0
    hasting_ratio = 0
    # proposing new parameter values
    if np.random.random() < update_frequency:
        # using the update_multiplier_proposal function we will get a new proposed parameter and a hastings ratio which we have to add to the posterior ratio in a later step
        proposed_shape_sp0,hast = update_multiplier_proposal(shape_init_sp0,update_factor)
        hasting_ratio += hast
        proposed_rate_sp0,hast = update_multiplier_proposal(rate_init_sp0,update_factor)
        hasting_ratio += hast        
        proposed_shape_sp1 = shape_init_sp1
        proposed_rate_sp1 = rate_init_sp1
    else:
        proposed_shape_sp0 = shape_init_sp0
        proposed_rate_sp0 = rate_init_sp0
        proposed_shape_sp1,hast = update_multiplier_proposal(shape_init_sp1,update_factor)
        hasting_ratio += hast
        proposed_rate_sp1,hast = update_multiplier_proposal(rate_init_sp1,update_factor)
        hasting_ratio += hast

    #print(hasting_ratio)
    #proposed_shape = update_normal(shape_init,d_shape)
    #proposed_rate = update_normal(rate_init,d_rate)
    # calculate likelihood ratio
    new_lik_sp0 = sum(log_pdf_gamma(data[species_ind ==0],proposed_shape_sp0,proposed_rate_sp0))
    new_lik_sp1 = sum(log_pdf_gamma(data[species_ind ==1],proposed_shape_sp1,proposed_rate_sp1))
    new_lik_joined = new_lik_sp0 + new_lik_sp1
    #likelihood_ratio = new_likelihood/current_lik
    
    # calculate new prior
    new_prior_shape_sp0 = log_pdf_gamma(proposed_shape_sp0,shape_prior_mu,shape_prior_sd)
    new_prior_rate_sp0 = log_pdf_gamma(proposed_rate_sp0,rate_prior_mu,rate_prior_sd)
    new_prior_shape_sp1 = log_pdf_gamma(proposed_shape_sp1,shape_prior_mu,shape_prior_sd)
    new_prior_rate_sp1 = log_pdf_gamma(proposed_rate_sp1,rate_prior_mu,rate_prior_sd)
    new_prior_joined = new_prior_shape_sp0 + new_prior_rate_sp0 + new_prior_shape_sp1 + new_prior_rate_sp1
    
    #prior_ratio = new_prior/current_prior_rate
    
    # calculate new posterior
    new_posterior_joined = new_lik_joined + new_prior_joined
    posterior_ratio_joined = new_posterior_joined - current_posterior_joined
    
    # draw a random number between 0 and 1
    random_draw = np.random.random()
    transformed_draw = np.log(random_draw)
    # see if that random draw lies within your acceptance ratio
    #print(proposed_rate_sp0,proposed_rate_sp1,new_lik_joined)
    if transformed_draw < posterior_ratio_joined + hasting_ratio:
        accept = True
        #print('accepted:',proposed_rate_sp0,proposed_rate_sp1,new_lik_joined)
        shape_init_sp0 = proposed_shape_sp0
        shape_init_sp1 = proposed_shape_sp1
        rate_init_sp0 = proposed_rate_sp0
        rate_init_sp1 = proposed_rate_sp1
        current_lik_joined = new_lik_joined
        current_prior_joined = new_prior_shape_sp0 + new_prior_rate_sp0 + new_prior_shape_sp1 + new_prior_rate_sp1
        current_posterior_joined = new_lik_joined + new_prior_joined
        acceptance_counter += 1
    else:
        accept = False

    # if the iteration can be divided by the user-input logging frequency, print to file
    if i % log_every_n_lines == 0:
        wlog.writerow([i,current_posterior_joined,current_prior_joined,current_lik_joined,shape_init_sp0,rate_init_sp0,shape_init_sp1,rate_init_sp1,shape_init_sp0/rate_init_sp0,shape_init_sp1/rate_init_sp1,shape_init_sp1/rate_init_sp1-shape_init_sp0/rate_init_sp0])
        logfile.flush() #this will tell python to write to the file whenever it is being stated in the code


print("Acceptance rate is", float(acceptance_counter)/float(n_iterations))   



# plotting the prior
#x = np.arange(0,100,0.1)
#plt.plot(x,pdf_gamma(x,1,.1))














#___________________________MESSING AROUND__________________________

def normal_dist(x,m,s):
    y = (1/(np.sqrt(2*np.pi*s**2)))*np.e**(-((x-m)**2)/(2*s**2))
    return y

import pandas as pd
mcmc_log = pd.read_csv("/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples_comparing_means_gamma.txt", sep = '\t')

plt.plot(mcmc_log.iloc[10:,0],mcmc_log.iloc[10:,1]), plt.title('post'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[10:,0],mcmc_log.iloc[10:,2]), plt.title('prior'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[10:,0],mcmc_log.iloc[10:,3]), plt.title('lik'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[10:,0],mcmc_log.iloc[10:,4]), plt.title('shape'),plt.show()#, plt.axis([0, 10000,90,110])
plt.plot(mcmc_log.iloc[10:,0],mcmc_log.iloc[10:,5]), plt.title('sig'),plt.show()#, plt.axis([0, 10000,90,110])


#_______________________________________DENSITY PLOT___________________________
import seaborn as sns

data_log = mcmc_log['shape'][1:]
sns.set_style('whitegrid')
sns.kdeplot(np.array(data_log), bw=4,color="green", shade=True)
plt.show()

data_log = mcmc_log.rate[1:]
sns.set_style('whitegrid')
sns.kdeplot(np.array(data_log), bw=0.05,color="blue", shade=True)
plt.show()

## plto a normal distribution
#plt.hist(mcmc_log.shape,bins=100)
#plt.plot(density(mcmc_log.shape))
#x_test = np.arange(50,95,1)
#y_test = normal_dist(x_test,70.883,6.371)
#plt.plot(x_test,y_test)
#

