import util
import run_sim as run
import plotting as pt

pd = util.params_dict


def main(pd, plot_stop_times=False,plot_all_spike_freq=False,plot_each_spike_freq=False):
	print('Running with Parameters:')
	for k,v in pd.items():
		print(f'{k}: {v}')

	if plot_stop_times:
		### Make mito freeze time histograms
		pt.mito_stop_hist(pd)

	if plot_all_spike_freq:
		### Make Raster plot for all Frequencies in one plot
		pt.spike_plot(pd)
	if plot_each_spike_freq:
		### Make Raster plots for each frequency
		for i in range((len(pd['event_hz']))):
			pt.spike_plot(pd, event_hz=[pd['event_hz'][i]])

	# collect simulation data
	data = run.run_sim(pd, num_runs=10)

	# plot simulation data
	pt.plot_pct_mm(data)

if __name__ == '__main__':
	main(pd, plot_stop_times=True,plot_all_spike_freq=False,plot_each_spike_freq=False)
