import numpy as np 

dend_L 		= 70
event_range = [1,2]                                 # um affected by synaptic events
event_min   = event_range[0]/dend_L                 # min percentage of dendrite affected by synaptic event
event_max   = event_range[1]/dend_L                 # max percentage of dendrite affected by synaptic event
event_mean  = event_min + (event_max - event_min)/2 # mean percentage of dendrite affected by synaptic event
event_sd    = (event_max - event_mean)/3   

density 	= 1/4
mito_pop 	= 500
MM_pct_init = 1

days 		= 5
day_seconds = 1500 

event_hz 	= [0.001, 0.012, 0.035, 0.058, 0.081, 0.104, 0.127, 0.15, 0.25, 0.35, 0.5]
recov_means = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25])*60
recov_sd 	= 2.5*60

params_dict = {'day_seconds':  day_seconds,
				'days':        days,
				'density':     density,      # density of mitochondria (number/um)
				'dend_L':      dend_L,       # length of dendrite (um)
				'event_range': event_range,  # um affected by synaptic events
				'event_min':   event_min,    # min percentage of dendrite affected by synaptic event
				'event_max':   event_max,    # max percentage of dendrite affected by synaptic event
				'event_mean':  event_mean,   # mean percentage of dendrite affected by synaptic event
				'event_sd':    event_sd,     # 99% of distr covered by 3sd
				'mito_pop':    mito_pop, 		 
				'MM_pct_init': MM_pct_init,  # proportion of mobile mitochondria at t = 0
				'recov_means': recov_means,  # Mean and sd of recovery times for normally distributed mitos
				'recov_sd':    recov_sd,
				'event_hz':    event_hz
}



## Define Mitochondrion Objects to populate list
class mitochondrion(object):
    def __init__(self, MM_stop, mobile):
        self.mobile 			= mobile
        self.time_immobilized 	= 0
        self.MM_stop 			= MM_stop
        self.freeze_counter 	= 0 
        
    def freeze(self, time):
        '''immobilize mitochondria 
        each object gets time at which it receives synaptic signal
        if already frozen, continues to freeze
        '''
        self.mobile 			= False
        self.time_immobilized 	= time
        self.freeze_counter		+=1
        
    def release(self,time):
        if time == self.time_immobilized+self.MM_stop+(self.freeze_counter):
            self.mobile = True
            self.time_immobilized = np.nan

# normal distribution
def populate_list_n(mito_pop, MM_pct_init, mean_stop, sd_stop):
    mito_list = []         # initialize list storing instances of mitochondrion class
    noisy_recovery = True # option for giving different mitochondria different recovery times (drawn from normal distribution around mean = fast_stop and sd = 0.1*fast_stop)

    for i in range(mito_pop):
        if noisy_recovery == True: 
            #normal 
            recovery_time = int(np.round(np.random.normal(loc=mean_stop, scale=sd_stop))) # distribution for determining recovery for mobile mitos (~seconds or minutes)
        else:
            recovery_time = mean_stop
        mobility = np.random.choice([0,1],1,p=[1-MM_pct_init,MM_pct_init])
        if recovery_time<0:
            recovery_time = 5
        mito_list.append(mitochondrion(recovery_time, mobility))
        
    return mito_list

# calculate fraction of mitochondria that are mobile            
def calc_frac_mm(mito_list):
    pop = len(mito_list)
    mobile_num = 0
    for k in mito_list:
        if k.mobile == True:
            mobile_num += 1
    mobile_frac = mobile_num/pop   
    return mobile_frac

# calculate a running mean over a list or array (effectively smooths noisy data)
def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result




### Junkyard

# bimodal distribution
def populate_list_bimodal(mito_pop, MM_pct_init, fast_stop, slow_stop):
    mito_list = []         # initialize list storing instances of mitochondrion class
    noisy_recovery = True # option for giving different mitochondria different recovery times (drawn from normal distribution around mean = fast_stop and sd = 0.1*fast_stop)

    for i in range(mito_pop):
        if noisy_recovery == True: 
            fast_recovery_time = int(np.round(np.random.normal(loc=fast_stop, scale=0.15*fast_stop))) # distribution for determining recovery for mobile mitos (~seconds or minutes)
            slow_recovery_time = int(np.round(np.random.normal(loc=slow_stop, scale=0.1*slow_stop))) # distribution for determining recovery time for immobile mitos (~day)

        else:
            fast_recovery_time = fast_stop
            slow_recovery_time = slow_stop
        mobility = np.random.choice([0,1],1,p=[1-MM_pct_init,MM_pct_init])
        if mobility == 0:
            recovery_time = slow_recovery_time
        else:
            recovery_time = fast_recovery_time
        mito_list.append(mitochondrion(recovery_time, mobility))
    return mito_list