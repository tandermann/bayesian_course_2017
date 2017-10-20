#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 15:21:56 2017

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
pandas_bm =np.array([84.74, 84.88, 94.60, 96.37, 102.93, 109.11, 125.76])
commute_length = np.array([63.5,81.,61.5,120.5,64.5,66.3,62.3,65.6,66.5,70.5,62.3,66.5,68.9])
data = pandas_bm

#___________________________________MCMC SETTINGS______________________________
# initiate log file
logfile = open("/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples_gamma_model_testing.txt","w")
wlog = csv.writer(logfile,delimiter='\t')
wlog.writerow(["it","post","prior","lik","shape","rate","mean_bm"])
logfile.flush() #this will tell python to write to the file whenever it is being stated in the code

log_every_n_lines = 10

# initial values for the MCMC
shape_init = 5.
rate_init = 5.

# parameters defining prior distribution
# prior probability of shape will be calculated from a normal distribution between a and b
shape_prior_mu = 1
shape_prior_sd = .1
# prior probability of rate is calculated from a gamma distribution with the parameters a and b
rate_prior_mu = 1
rate_prior_sd = .1

#MCMC settings
n_iterations = 10000
# 'window size' for new proposals
d_shape = 1
d_rate = .1
update_factor = 1.5

# calculating the initial likelihood and prior probability
current_lik = sum(log_pdf_gamma(data,shape_init,rate_init))
current_prior_shape = log_pdf_gamma(shape_init,shape_prior_mu,shape_prior_sd)
current_prior_rate = log_pdf_gamma(rate_init,rate_prior_mu,rate_prior_sd)
current_prior = current_prior_shape + current_prior_rate
current_posterior = current_lik + current_prior

# count how many steps we accepted
acceptance_counter = 0

# this ratio defines how often parameter 1 (shape) is sampled in relation to parameter 2 (rate)
update_frequency = .5

# turn the factor proposal function on or off
multiplier_proposal = True

# set the steps (temperature) for the model testing
temperatures = np.array([1.,0.64,0.32,0.16,0.08,0.04,0.02,0.01,0.005,0.001])
temp_ind = -1

posterior_mean_dict = {}

for i in range(0,n_iterations*len(temperatures)):
    if i % n_iterations == 0:
        temp_ind += 1
        if temp_ind != 0:
            average_posterior_heatstep = sum(posterior_list)/len(posterior_list)
            posterior_mean_dict.setdefault(temperature,average_posterior_heatstep)
        posterior_list = []
    #print(temp_ind)
    temperature = temperatures[temp_ind]
    # initiate a hastings ratio of 0
    hasting_ratio = 0
    hasting_ratio_add = 0
    # proposing new parameter values
    if multiplier_proposal:
        if np.random.random() < update_frequency:
            # using the update_multiplier_proposal function we will get a new proposed parameter and a hastings ratio which we have to add to the posterior ratio in a later step
            proposed_shape,hasting_ratio_add = update_multiplier_proposal(shape_init,update_factor)
            proposed_rate = rate_init
        else:
            proposed_shape = shape_init
            proposed_rate,hasting_ratio_add = update_multiplier_proposal(rate_init,update_factor)
    # if multiplier_proposal==False, the new parameters should be updated by a fixed normal distirbution
    else:
        # proposing new parameter values  
        if np.random.random() < update_frequency:
            proposed_shape = update_normal(shape_init,d_shape)
            proposed_rate = rate_init
        else:
            proposed_shape = shape_init
            proposed_rate = update_normal(rate_init,d_rate)
    # update the hastings ratio
    hasting_ratio = hasting_ratio + hasting_ratio_add
    #print(hasting_ratio)
    #proposed_shape = update_normal(shape_init,d_shape)
    #proposed_rate = update_normal(rate_init,d_rate)
    # calculate likelihood ratio
    new_likelihood = sum(log_pdf_gamma(data,proposed_shape,proposed_rate))
    likelihood_ratio = new_likelihood-current_lik
    
    # calculate prior ratio
    new_prior_shape = log_pdf_gamma(proposed_shape,shape_prior_mu,shape_prior_sd)
    new_prior_rate = log_pdf_gamma(proposed_rate,rate_prior_mu,rate_prior_sd)
    new_prior = new_prior_shape + new_prior_rate
    prior_ratio = new_prior - current_prior
    
    # calculate new posterior
    #new_posterior = new_likelihood + new_prior
    posterior_ratio = temperature*(likelihood_ratio) + prior_ratio
    
    # draw a random number between 0 and 1
    random_draw = np.random.random()
    transformed_draw = np.log(random_draw)
    # see if that random draw lies within your acceptance ratio
    if transformed_draw < (posterior_ratio + hasting_ratio):
        accept = True
        shape_init = proposed_shape
        rate_init = proposed_rate
        current_lik = new_likelihood
        current_prior = new_prior
        current_posterior = new_likelihood+new_prior
        acceptance_counter += 1
    else:
        accept = False

    # if the iteration can be divided by the user-input logging frequency, print to file
    if i % log_every_n_lines == 0:
        posterior_list.append(current_posterior)
        #print (i,current_posterior,current_prior,current_lik,shape_init,rate_init,shape_init/rate_init)
        wlog.writerow([i,current_posterior,current_prior,current_lik,shape_init,rate_init,shape_init/rate_init])
        logfile.flush() #this will tell python to write to the file whenever it is being stated in the code
# get the last mean of posterior for the temperature at the end of the loop
last_temp_mean = sum(posterior_list)/len(posterior_list)
posterior_mean_dict.setdefault(temperature,last_temp_mean)

print("Acceptance rate is", (float(acceptance_counter)/float(n_iterations))/len(temperatures))   



# plotting the prior
#x = np.arange(0,100,0.1)
#plt.plot(x,pdf_gamma(x,1,.1))




# plot the heat steps and the eman of their posterior distribution
x_heat = list(posterior_mean_dict.keys())
y_heat = list(posterior_mean_dict.values())
f=plt.figure()
plt.plot(x_heat,y_heat),plt.plot(x_heat,y_heat,'ro'),[plt.axvline(x=i,color='grey',linestyle='--') for i in x_heat]
f.savefig('/Users/tobias/GitHub/bayesian_modeling_course/data/heating_chains.png', dpi = 500)
plt.show()




#___________________________MESSING AROUND__________________________

def normal_dist(x,m,s):
    y = (1/(np.sqrt(2*np.pi*s**2)))*np.e**(-((x-m)**2)/(2*s**2))
    return y

import pandas as pd
mcmc_log = pd.read_csv("/Users/tobias/GitHub/bayesian_modeling_course/data/py_mcmc_samples_gamma_model_testing.txt", sep = '\t')

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

