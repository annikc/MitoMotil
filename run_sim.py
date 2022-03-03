from __future__ import division
import os
import numpy as np
from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms
import util
pd = util.params_dict

np.random.seed(1234)

load_data = False
if load_data:
	dose_response = np.load('condition_1/doseresponse.npy').item()
else:
	dose_response = {}

def run_sim(pd, num_runs, **kwargs):
	event_hz = kwargs.get('event_hz', pd['event_hz']) ## must be list
	recov_means = kwargs.get('recov_means', pd['recov_means']) ## must be list
	savedata = kwargs.get('savedata', True)

	track_pct_mm_dict = {}
	for i in range(len(recov_means)):
		if savedata:
			directory = './data/MeanRecovTimes/{}sec/'.format(recov_means[i])
			if not os.path.exists(directory):
				os.makedirs(directory)
			savedir = kwargs.get('savedir', directory)
		print( "     Series of {} runs for mean mito reocvery of {}".format(num_runs*len(event_hz), recov_means[i]) )
		# set up dictionary to store days data for each frequency
		track_pct_mm_frq = {}
		# step over different frequencies
		for j in range(len(event_hz)):
			event_freq = event_hz[j]
			print('Beginning {} runs at {} Hz:'.format(num_runs, event_freq))
			spiketrain_list = []
			
			# set up list to store mobile fraction data for each day
			track_pct_mm = []
			#step over number of example runs
			for k in range(num_runs):
				# Create list of mitochondria objects with freezing times drawn from distribution 
				mito_list = util.populate_list_n(pd['mito_pop'], pd['MM_pct_init'], recov_means[i], pd['recov_sd'])
		
				pct_mm = []
				day = homogeneous_poisson_process(rate    = event_freq*Hz,
												  t_start = 1*s,
												  t_stop  = pd['day_seconds']*s, 
												  as_array= True)
				spiketrain_list.append(day)

				synaptic_event_times = spiketrain_list[k]
				#step over seconds of the day 
				for t in range(pd['day_seconds']):
					if t in [int(np.round(x)) for x in synaptic_event_times]:
						#generate percentage affected by syn event
						frac_immobilized = np.random.normal(loc=pd['event_mean'], scale=pd['event_sd'])
						#select mitochondria to immobilize
						immobilized_ind = np.random.choice(pd['mito_pop'], int(np.round(pd['mito_pop']*frac_immobilized)))
						for m in immobilized_ind:
							mito_list[m].freeze(t)
					for m in range(pd['mito_pop']):
						mito_list[m].release(t)

					mobile_frac = util.calc_frac_mm(mito_list)
					pct_mm.append(mobile_frac)

				# save mobile fraction data for the day 
				track_pct_mm.append(pct_mm)
				if k == 0:
					print( "Finished 1...", end=" ")
				elif k < num_runs-1:
					print ("{}...".format(k+1), end=" ")
				else:
					print('{}'.format(k+1))

			# store data for each frequency 
			track_pct_mm_frq[str(event_hz[j])] = track_pct_mm


		#store data for each mito recovery average condition
		track_pct_mm_dict[str(recov_means[i])] = track_pct_mm_frq
		if savedata: 
			np.save(savedir+f'{recov_means[i]}_{event_hz[0]}to{event_hz[-1]}_data.npy',track_pct_mm_frq)
		print( "Done runs for mean mito recovery of {} \n==================".format(recov_means[i]))

	return(track_pct_mm_dict)        