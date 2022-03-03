# import necessary packages
import numpy as np
from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms
from util import params_dict as pd


def get_spike_trains(num_runs, event_hz, freq_block=pd['day_seconds']):
    ### Create Spike trains from poisson process
    spiketrain_list = []
    for i in range(num_runs):
        spike_holder = []
        for j in range(len(event_hz)):
            spikes = homogeneous_poisson_process(rate=event_hz[j]*Hz,
                                                   t_start=(j*freq_block)*s+0.01*s,
                                                   t_stop=(j+1)*freq_block*s, as_array=True)
            spike_holder = np.concatenate([spikes, spike_holder])

        spiketrain_list.append(spike_holder)
    return spiketrain_list