from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from elephant.spike_train_generation import homogeneous_poisson_process
from quantities import Hz, s, ms
from functions import get_spike_trains
import util

pd = util.params_dict


def spike_plot(pd, **kwargs):
	'''
	Spike Raster Plot
	Create Spike trains from poisson process
	'''
	event_hz = kwargs.get('event_hz', pd['event_hz'])
	num_runs = kwargs.get('num_runs', 10)
	savedir  = kwargs.get('savedir', './data/SpikeRasters/')
	savefig  = kwargs.get('savefig', True)
	showfig  = kwargs.get('showfig', False)
	format   = kwargs.get('format','png')

	sns.set_palette("Spectral",num_runs) # Set color palette for plots

	if savefig:
		if not os.path.exists(savedir):
			os.makedirs(savedir)

	print(f'Plotting Spike Rasters for {event_hz} Hz Inputs')

	spiketrain_list = get_spike_trains(num_runs, event_hz, freq_block=pd['day_seconds'])


	color_list = [tuple((np.array(x)+0.85)%1) for x in sns.color_palette("Spectral", num_runs)]

	plt.figure()
	plt.plot(np.zeros(pd['day_seconds']), 'o', markersize=0)
	for i, spiketrain in enumerate(spiketrain_list):
		t = spiketrain
		plt.vlines(t, i+0.6, i+1.4, color = color_list[i])

	
	plt.xlabel('Time (sec)')
	plt.ylabel('Example Spike Trains         ')
	plt.gca().tick_params(axis='both', which='major', labelsize=10)

	if len(event_hz) > 1: 
		plt.plot((num_runs+3.5)*np.ones(pd['day_seconds']), 'o', markersize=0)
		for i in range(len(event_hz)-1):
			plt.axvline(x=(i+1)*pd['day_seconds'], color='k', alpha=0.2)

		for i in range(len(event_hz)):
			x = ((i)*pd['day_seconds'] + pd['day_seconds']/3)
			y = num_runs + 3.25
			plt.annotate('{0:0.3f} Hz'.format(event_hz[i]),(x,y), rotation = 80)
			plt.yticks(np.arange(num_runs)+1, tuple(['']*num_runs) )
		
		if savefig:
			plt.savefig(savedir+f'Raster{len(event_hz)}_{event_hz[0]}-{event_hz[-1]}.{format}', format=format)
	else:
		plt.plot((num_runs+1)*np.ones(pd['day_seconds']), 'o', markersize=0)
		for i in range(len(event_hz)-1):
			plt.axvline(x=(i+1)*pd['day_seconds'], color='k', alpha=0.2)

		plt.suptitle('Spike Raster Plot for {} Hz'.format(event_hz[0]), fontsize=12)
	
		if savefig:
			plt.savefig(f'./data/SpikeRasters/Raster{event_hz[0]}.{format}', format=format)

	if showfig:
		plt.show()
	plt.close()

def mito_stop_hist(pd, n_bins=25, **kwargs):
	# make example histogram of recovery times
	recov_means = kwargs.get('recov_times', pd['recov_means'])
	savedir     = kwargs.get('savedir', './data/Histograms/')
	savefig     = kwargs.get('savefig', True)
	showfig     = kwargs.get('showfig', False)
	format      = kwargs.get('format','png')
	fixed_scale = kwargs.get('fixed_scale', True)

	if savefig:
		if not os.path.exists(savedir):
			os.makedirs(savedir)

	sns.set_palette("Spectral",10) # Set color palette for plots
	color_list = [tuple((np.array(x)+0.85)%1) for x in sns.color_palette("Spectral", 10)]

	print('Plotting Mitochondria Stopping Time Histograms...')
	for i in range(len(recov_means)):
		mean_stop = pd['recov_means'][i]
		mito_list = util.populate_list_n(pd['mito_pop'], pd['MM_pct_init'], mean_stop, pd['recov_sd'])
		recov_times = []
		for i in range(pd['mito_pop']):
			recov_times.append(mito_list[i].MM_stop)

		quick_mitos = [x for x in recov_times if x<10]
		
		plt.figure()

		plt.hist(recov_times, bins=n_bins, color = color_list[-1])
		
		if fixed_scale:
			plt.ylim([0,200])

		plt.ylabel('Number of Mitochondria')
		plt.xlabel('Freezing Time (sec)')
		plt.suptitle('Mitochondria Freezing Times ($\mu$={}, $\sigma^2$={} sec) '.format(mean_stop, int(pd['recov_sd'])))

		if savefig:
			plt.savefig(savedir+f'Mito_freezing_histogram_{mean_stop}sec.{format}', format=format)
		if showfig:
			plt.show()
		plt.close()

def plot_pct_mm(data, **kwargs):
	savefig = kwargs.get('savefig', True)
	showfig = kwargs.get('showfig', False)
	for mean_freeze in data.keys():
		if savefig:
			directory = './data/MeanRecovTimes/{}sec/'.format(mean_freeze)
			if not os.path.exists(directory):
				os.makedirs(directory)
			savedir = kwargs.get('savedir', directory)

		num_runs = len(data)
		for event_hz in data[str(mean_freeze)].keys():
			track_pct_mm = data[mean_freeze][event_hz]
			num_runs = len(track_pct_mm)

			average_pct_mm = []
			for i in range(len(track_pct_mm[0])):
				avg_sum = 0
				for j in range(num_runs):
					avg_sum += track_pct_mm[j][i]
				avg = avg_sum/num_runs
				average_pct_mm.append(avg)

			plt.figure()
			sns.set_palette("Spectral",num_runs) # Set color palette for plots
			color_list = [tuple((np.array(x)+0.85)%1) for x in sns.color_palette("Spectral", num_runs)]

			for i in range(num_runs):
				plt.plot(track_pct_mm[i], color=color_list[i],alpha = 0.5)
			plt.plot(average_pct_mm, 'k', label = 'Average')

			plt.suptitle('Mitochondrial Motility with {} Hz Stimulation \n(Mean Freeze Time {} sec)'.format(event_hz, mean_freeze), fontsize = 12)
			plt.ylabel('Mobile Fraction of Mitochondria')
			plt.xlabel('Time (sec)')
			plt.legend(loc=0)
			plt.ylim([0,1.05])
			if savefig:
				plt.savefig(savedir+'AvgMitoMotility{}.png'.format(event_hz), format='png')
			if showfig:
				plt.show()
			
			plt.close()


## for specific use case
def plot_spike_raster(spiketrain_list, day_seconds, event_hz, **kwargs):
	### Spike Raster Plot
	num_runs = len(spiketrain_list)
	sns.set_palette("Spectral",num_runs) # Set color palette for plots
	color_list = [tuple((np.array(x)+0.85)%1) for x in sns.color_palette("Spectral", num_runs)]

	savefig = kwargs.get('savefig', False)
	file_format  = kwargs.get('fileformat','png')

	plt.figure()
	plt.plot(np.zeros(day_seconds), 'o', markersize=0)
	for i, spiketrain in enumerate(spiketrain_list):
		t = spiketrain
		plt.vlines(t, i+0.6, i+1.4, color = color_list[i])

	plt.plot((num_runs+1)*np.ones(day_seconds), 'o', markersize=0)

	for i in range(len(event_hz)):
		x = ((i)*day_seconds + day_seconds/3)
		y = num_runs + 1
		plt.annotate('{0:0.3f} Hz'.format(event_hz[i]),(x,y+0.25), rotation = 0)
		rect = mpatches.Rectangle((i*day_seconds, y),day_seconds,1, color='gray', alpha=0.2)
		plt.gca().add_patch(rect)

	plt.xlabel('Time (sec)')
	plt.ylabel('Spike Train Index')
	plt.gca().tick_params(axis='both', which='major', labelsize=10)

	if savefig:
		plt.savefig(f'SpikeRaster.{file_format}', format=file_format)
	plt.show()
	plt.close()
