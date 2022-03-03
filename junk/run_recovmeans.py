from __future__ import division, print_function
import os, sys, time

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms

import util as mm

from params import *


#### Run full event hz for dictionary population
recov_means = np.array([5])*60 #np.array([1,2,3,4,5,6,7,8,9,10,15,20,25])*60

event_hzs = [0.001, 0.012, 0.035, 0.058, 0.081, 0.104, 0.127, 0.15, 0.25, 0.35, 0.5]

dose_response = {}

for mean_stop in recov_means:
	mean_stop_min = int(mean_stop/60)
	print ("Mean Mito Immobilization = {0} min (Avg {1} runs/ stim freq)".format(mean_stop_min, num_runs), '\n========')
	mainpath = 'data/MeanRecovTimes/condition_'
	if not os.path.isdir(mainpath+'{}'.format(mean_stop_min)):
		os.makedirs(mainpath+'{}'.format(mean_stop_min))
	savedir = mainpath+'{}/'.format(mean_stop_min)
	
	for ev_hz in event_hzs:
		event_hz = [ev_hz]

		### Create Spike trains from poisson process
		spiketrain_list = []
		plot_rasters = []

		for i in range(num_runs):
			for j in range(len(event_hz)):
				day = homogeneous_poisson_process(rate=event_hz[j]*Hz,
													   t_start=0.1*s,
													   t_stop=day_seconds*s, as_array=True)
				c1 = np.concatenate([[0],day])
				plot_rasters.append(np.concatenate([c1, [day_seconds-1]]))
				if event_hz[j] == 0:
					day= np.concatenate([day,[3*day_seconds - 1]])
				a = day+(j*day_seconds)
				if j == 0: 
					day_holder = np.concatenate([[0,1],a])
				else:
					day_holder = np.concatenate([day_holder,a])

			spiketrain_list.append(day_holder)

		### Spike Raster Plot
		if mean_stop == recov_means[0]:
			extra_height = 2
			plt.figure()
			plt.plot(np.zeros(day_seconds), 'o', markersize=0)
			for i, spiketrain in enumerate(spiketrain_list):
				#if i == 0: 
				#    plt.plot((i)*np.ones(len(event_hz)*day_seconds), 'o', markersize=0)
				t = spiketrain #.rescale(s)
				plt.plot(t, (i+1)* np.ones_like(t), marker='o', markersize=10, color = color_list[i])
				#if i == len(spiketrain_list)-1:
				#    plt.plot((len(spiketrain_list)+1)*np.ones(len(event_hz)*day_seconds), 'o', markersize=0)

			plt.plot((num_runs+extra_height)*np.ones(day_seconds), 'o', markersize=0)
			for i in range(len(event_hz)-1):
				plt.axvline(x=(i+1)*day_seconds, color='k', alpha=0.2)

			for i in range(len(event_hz)):
				x = ((i)*day_seconds + day_seconds/3)
				y = num_runs + extra_height
				plt.annotate('{0:0.3f} Hz'.format(event_hz[i]),(x,y), rotation = 90)
			plt.axis('tight')
			#plt.xlim([0,10])
			plt.xlabel('Time (sec)')
			plt.ylabel('Spike Train Index')
			plt.gca().tick_params(axis='both', which='major', labelsize=10)
			#plt.suptitle('Spike Raster Plot', fontsize=12)
			plt.savefig('data/SpikeRasters/SpikeRaster{}.png'.format(event_hz[0]), format='png')
			plt.savefig('data/SpikeRasters/SpikeRaster{}.svg'.format(event_hz[0]), format='svg')
			plt.close()

		mito_list = mm.populate_list_n(mito_pop, MM_pct_init, mean_stop, sd_stop)
		recov_times = []
		for i in range(mito_pop):
			recov_times.append(mito_list[i].MM_stop)

		plt.hist(recov_times, bins=25, color = color_list[0])
		#plt.xlim([0, days*day_seconds])
		plt.ylabel('Number of Mitochondria')
		plt.xlabel('Recovery Time (sec)')
		#plt.suptitle('Distribution of Mitochondria \nRecovery Times (sec)')
		plt.savefig(savedir+'/recov_distribution.png'.format(mean_stop_min), format='png')
		plt.savefig(savedir+'/recov_distribution.svg'.format(mean_stop_min), format='svg')
		plt.close()




		#### RUN SIMULATION
		recovery_times = np.zeros(num_runs)
		track_pct_mm = []
		print ("@{} Hz".format(event_hz[0]))
		for i in range(num_runs):
			recovery_time=0
			pct_mm = []
			synaptic_event_times = spiketrain_list[i]
			#normally distributed stopping times
			mito_list = mm.populate_list_n(mito_pop, MM_pct_init, mean_stop, sd_stop)
			#bimodally distributed stopping times
			#mito_list = mm.populate_list_bimodal(mito_pop, MM_pct_init, fast_stop, slow_stop)
			for j in range(len(event_hz)*day_seconds):
				if j in [int(np.round(x)) for x in synaptic_event_times]:
					#generate percentage affected by syn event
					frac_immobilized = np.random.normal(loc=event_mean, scale=event_sd)
					#select mitochondria to immobilize
					immobilized_ind = np.random.choice(mito_pop, int(np.round(mito_pop*frac_immobilized)))
					for k in immobilized_ind:
						mito_list[k].freeze(j)

				for k in range(mito_pop):
					mito_list[k].release(j)
				# track fraction of mitochondria that are mobile
				mobile_frac = mm.calc_frac_mm(mito_list)  
				pct_mm.append(mobile_frac)
				if event_hz[int(j/day_seconds)] == 0:
					rec_recov = True
					if rec_recov:
						recovery_time +=1
					if mobile_frac == pct_mm[0]:
						recovery_times[i]=recovery_time
						rec_recov = False

			track_pct_mm.append(pct_mm)
			print ("{0}...".format(i+1),end='')
			sys.stdout.flush()
			time.sleep(0.3)
		print ("Done")

		average_pct_mm = []
		for i in range(len(track_pct_mm[0])):
			avg_sum = 0
			for j in range(num_runs):
				avg_sum += track_pct_mm[j][i]
			avg = avg_sum/num_runs
			average_pct_mm.append(avg)

		dose_response['{}'.format(event_hz[0])] = mm.running_mean(average_pct_mm, 25)[-1]
		np.save(savedir+'doseresponse.npy',dose_response)
		
		plt.figure()
		for i in range(num_runs):
			plt.plot(track_pct_mm[i], alpha = 0.5)
		plt.plot(average_pct_mm, 'k', label = 'Average')
		for i in range(len(event_hz)):
			plt.axvline(x=(i+1)*day_seconds, color='k', alpha=0.2)
		plt.legend(loc=0)
		for i in range(len(event_hz)):
			x = ((i)*day_seconds + day_seconds/2)
			y = 0.08
			plt.annotate('{0:.3f}Hz'.format(event_hz[i]),(x,y), rotation=90)
		#plt.ylim([0,0.1])
		#plt.xlim([0,days*day_seconds+10])
		plt.ylabel('Mobile Fraction of Mitochondria')
		plt.xlabel('Time (sec)')
		plt.suptitle('Average Mitochondrial Motility Over {} Runs'.format(num_runs), fontsize = 12)
		plt.savefig(savedir+'AvgMitoMotility{}.png'.format(event_hz[0]), format='png')
		#plt.show()
		plt.close()