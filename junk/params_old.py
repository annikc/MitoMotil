from __future__ import division
import numpy as np
import seaborn as sns

np.random.seed(1234) 								# set seed for pseudorandom number generation in python

# ===============================
#      SIMULATION PARAMETERS
# ===============================
num_runs    = 10                                    # number of runs to average over 

day_seconds = 1500									# number of seconds per epoch
days        = 5										# number of epochs

density     = 1/4                                   # density of mitochondria (number/um)
dend_L      = 70                                    # length of dendrite (um)

event_range = [1,2]                                 # um affected by synaptic events
event_min   = (event_range[0]/dend_L)               # min percentage of dendrite affected by synaptic event
event_max   = (event_range[1]/dend_L)               # max percentage of dendrite affected by synaptic event
event_mean  = event_min + (event_max - event_min)/2 # mean percentage of dendrite affected by synaptic event
event_sd    = (event_max - event_mean)/3            # 99% of distr covered by 3sd

mito_pop    = 500 								    # population of mitos 
MM_pct_init = 1                                     # proportion of mobile mitochondria at t = 0

mean_stop   = 7*60									# mean for distribution of mito immobilization times 
sd_stop     = (2.5)*60 								# sd for distribution of mito immobilization times


# PLOTTING STYLE PARAMETERS
sns.set_palette("Spectral",num_runs)
# Set color palette for plots
color_list = [tuple((np.array(x)+0.85)%1) for x in sns.color_palette("Spectral", num_runs)]
