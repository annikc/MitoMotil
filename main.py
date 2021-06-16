import util
import run_sim as run
import plotting as pt 

pd = util.params_dict
print(pd)

### Make mito freeze time histograms
pt.mito_stop_hist(pd)

### Make Scatterplot for all Frequencies in one plot 
pt.spike_plot(pd)

### Make Scatterplots for each frequency
for i in range((len(pd['event_hz']))):
	pt.spike_plot(pd, event_hz=[pd['event_hz'][i]])

data = run.run_sim(pd)

pt.plot_pct_mm(data)