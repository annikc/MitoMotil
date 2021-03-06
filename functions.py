# import necessary packages
import numpy as np
from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms
import util

pd = util.params_dict

def get_spike_trains(num_runs, event_hz, freq_block=pd['day_seconds']):
    ### Create Spike trains from poisson process
    spiketrain_list = []
    for i in range(num_runs):
        spike_holder = []
        for j in range(len(event_hz)):
            spikes = homogeneous_poisson_process(rate=event_hz[j]*Hz,
                                                   t_start=(j*freq_block)*s+0.01*s,
                                                   t_stop=(j+1)*freq_block*s, as_array=True)
            spike_holder = np.concatenate([spike_holder,spikes], axis=None)
        spiketrain_list.append(spike_holder)
    return spiketrain_list

# normal distribution
def populate_list_n(mito_pop, MM_pct_init, mean_stop, sd_stop, noisy_recovery=True):
    # includes noisy_recovery option for giving different mitochondria different recovery times (drawn from normal distribution around mean = fast_stop and sd = 0.1*fast_stop)
    mito_list = []         # initialize list storing instances of mitochondrion class

    for i in range(mito_pop):
        if noisy_recovery == True:
            #normal
            recovery_time = int(np.round(np.random.normal(loc=mean_stop, scale=sd_stop))) # distribution for determining recovery for mobile mitos (~seconds or minutes)
        else:
            recovery_time = mean_stop
        mobility = np.random.choice([0,1],1,p=[1-MM_pct_init,MM_pct_init])
        if recovery_time<0:
            recovery_time = 5
        mito_list.append(util.mitochondrion(recovery_time, mobility))

    return mito_list

# calculate fraction of mitochondria that are mobile
def calc_frac_mm(mito_list):
    pop = len(mito_list)
    mobile_num = 0
    for k in mito_list:
        if k.mobile:
            mobile_num += 1
    mobile_frac = mobile_num/pop
    return mobile_frac

def immobilize_mitos(pd, mito_list, synaptic_event_times, experiment_time=pd['day_seconds']*len(pd['event_hz'])):
    pct_mm = []
    for t in range(experiment_time):
        if t in [int(np.round(event)) for event in synaptic_event_times]:
            #generate percentage affected by syn event
            frac_immobilized = np.random.normal(loc=pd['event_mean'], scale=pd['event_sd'])
            #select mitochondria to immobilize
            immobilized_ind = np.random.choice(pd['mito_pop'], int(np.round(pd['mito_pop']*frac_immobilized)))

            for m in immobilized_ind:
                mito_list[m].freeze(t)
        for m in range(pd['mito_pop']):
            mito_list[m].release(t)

        mobile_frac = calc_frac_mm(mito_list)
        pct_mm.append(mobile_frac)

    return np.array(pct_mm)