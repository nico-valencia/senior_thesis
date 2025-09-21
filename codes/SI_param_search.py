import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy.optimize as sci 
# from astropy.stats import sigma_clip
from multiprocessing import Pool 
import itertools 

# REMEMBER I DISABLED THE WARNING REGARDING A value is trying to be set on a copy of a slice from a DataFrame
pd.options.mode.chained_assignment = None  # Disables the warning



# set up parameter space to explore
param_grid = {
                'snr'           : np.around(np.linspace(4,10,7),3), 
                # 'snrbefore'     : np.around(np.logspace(0,10,10),3), 
                # 'snrafter'      : np.around(np.logspace(0,10,10),3), 
                'stalta'        : np.around(np.linspace(2,10,9),3), 
                # 'sierrmax'      : np.around(np.logspace(0,1,10),3), 
                # 'quality'       : np.around(np.logspace(0,10,10),3), 
                'sksweight'     : np.around(np.linspace(0.1,1,6),3), 
                'sweight'       : np.around(np.linspace(0.1,1,6),3), 
                'sdepthweight'  : np.around(np.linspace(0.1,1,6),3), 
                    }
param_combination = list(itertools.product(*param_grid.values()))

# station and method selection
# station = 'FORT'
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_mar12.txt',sep='|')
# ravel = False
# single = False
# sierrmax = 0.5
# s_depth = 400

station = 'STKA'
eventdataa = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_apr1.txt',sep='|')


def SI_finder(p,eventdata,station,mode):
    try:
        mode = mode
        snr_threshold = p[0] 
        sta_lta_threshold = p[1]
        sks_weight = p[2]
        s_depth_weight = p[3]
        s_weight = p[4]

        # station = 'DAV'
        # eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_apr1.txt',sep='|')
        ravel = False
        single = False
        sierrmax = 1
        s_depth = 400

        # # actually truncate the dataset lmao 
        eventdata = eventdataa.copy()
        eventdata = eventdata[eventdata['snr'] > snr_threshold]
        eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]


        if ravel: 
            eventdata['calcbaz'] = eventdata['calcbaz'] % 180
            degs = np.arange(0,180,1)
        else:
            degs = np.arange(0,360,1)

        rads = np.radians(degs)




        eventdata = eventdata.sort_values(by='calcbaz', ascending=True)

        # event si errors so they can be filtered out

        eventdata['SI_act_error'] = abs(eventdata['wiggleSI_act']-eventdata['wiggleSIlow_act'])
        eventdata['SI_calc_error'] = abs(eventdata['wiggleSI_calc']-eventdata['wiggleSIlow_calc'])

        eventdata = eventdata[eventdata['SI_act_error'] < sierrmax]
        eventdata = eventdata[eventdata['SI_calc_error'] < sierrmax]


        # eventdata = eventdata.dropna()[
        # eventdata = eventdata[eventdata['wiggleSI_calc'] != 0]
        eventdata = eventdata[eventdata['SI_calc_error'] != 0]


        ### quality factor criteria (mar 12 2025) ###
        '''
        curve_fit sigma parameter: inverse weight (closer to 0 is better)
        terms to consider: SI error , SNR , phase type , STA/LTA ... 
        '''

        inv_error = (1 / eventdata['SI_calc_error'])                # smaller error is better (so inverted)
        snr = eventdata['snrbefore'] + eventdata['snrafter']    # larger value is better 
        stalta = eventdata['sta_lta']


        eventdata['phase_weight'] = 0 
        eventdata['quality'] = 0 
        for i in range(0,len(eventdata)):
            if eventdata['evphase'].iloc[i] == 'SKS':
                eventdata['phase_weight'].iloc[i] = sks_weight
            if eventdata['evphase'].iloc[i] == 'S':
                if eventdata['evdp'].iloc[i] >= s_depth:
                    eventdata['phase_weight'].iloc[i] = s_depth_weight
                else:
                    eventdata['phase_weight'].iloc[i] = s_weight


        eventdata['quality'] = np.sqrt((inv_error**2)+(snr**2)+(stalta**2)) * eventdata['phase_weight']

        # quality_threshold = np.percentile(eventdata['quality'],99)
        # eventdata = eventdata[eventdata['quality'] <= quality_threshold]

        eventdata = eventdata.dropna(subset=['quality'])
        # eventdata = eventdata[eventdata['quality']>3]

        inv_weight = 1 /eventdata['quality']


        if single:
            def SI_function(baz,fast,lag):
                SI = lag * 2 * np.sin(np.deg2rad(baz-fast)) * np.cos(np.deg2rad(baz-fast))
                return SI
            bounds = ([-360, 0], [360,3])
            # p0 = [90.0 , 1.0]

            # raw measurements 
            sine_fit = sci.curve_fit(SI_function,eventdata['calcbaz'],eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
            params,cov,fitdata,trash1,trash2 = sine_fit
            phi,dt = params[0],params[1]


        ### USE MULTI DIMENSIONAL INPUT ARRAY FOR VARIABLE TO GET DATA IN 
        else:
            def SI_function(baz_incident,fast,lag,tilt):
                baz,incident = baz_incident
                SI = lag * ((np.cos(np.deg2rad(tilt))**2 * np.sin(2*np.deg2rad(baz-fast))) 
                            + (np.sin(2*np.deg2rad(tilt)) * np.sin(np.deg2rad(baz-fast)) * np.tan(np.deg2rad(incident))))
                return SI
            bounds = ([-360, 0, 0], [360, 3, 90])

            # raw measurements 
            baz_and_incident = np.vstack((eventdata['calcbaz'], eventdata['evincident']))
            sine_fit = sci.curve_fit(SI_function,baz_and_incident,eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
            params,cov,fitdata,trash1,trash2 = sine_fit
            phi,dt,tilt = params[0],params[1],params[2]


        # print(params)
        print(f"Total number of waveforms used: {len(eventdata)}")
        print(f"Fast Angle: {round(phi,4)} Delay Time {round(dt,4)}")
        if single != True:
            print(f"Tilt Angle: {round(tilt,4)}")

        # print("Covariance Matrix:")
        # print(cov)

        
        # print(len(eventdata))
        # print(len(fitdata['fvec']))
        if mode == 'covariance':
            print("Covariance Matrix Diagonal Elements:")
            # print(np.diag(cov))
            print(np.sum(np.diag(cov)**2))
            return np.sum(np.diag(cov)**2)
        if mode == 'SI':
            # print(f"Fast Angle: {round(phi,4)} Delay Time {round(dt,4)}")
            return round(phi,4),round(dt,4),round(tilt,4)
    except Exception as exception: 
        return np.inf



def parallel_grid_search(param_combination):
    return param_combination, SI_finder(param_combination,eventdataa,station,'covariance')

def parallel_grid_split(param_combination):
    return param_combination, SI_finder(param_combination,eventdataa,station,'SI')

if __name__ == '__main__':
    # Use multiprocessing Pool to distribute tasks across multiple processes
    with Pool(processes=64) as pool:
        # Perform parallel grid search
        results = pool.map(parallel_grid_search, param_combination)
        si_results = pool.map(parallel_grid_split, param_combination )

        # # Aggregate results
        best_params, best_performance = None, float('inf')
        data = np.zeros([len(results),6])
        si_data = np.zeros([len(si_results),8])
        # print(results)

        for i in range(0,len(results)):
            params , performance = results[i] 
            
            if performance < best_performance:
                best_params, best_performance = params, performance

            data[i] = params[0],params[1],params[2],params[3],params[4],performance

        for i in range(0,len(si_results)):
            try:
                params, vals = si_results[i]
                print(vals)
                si_data[i] = params[0],params[1],params[2],params[3],params[4], vals[0],vals[1],vals[2]
            except Exception as e:
                si_data[i] = params[0],params[1],params[2],params[3],params[4], np.inf,np.inf,np.inf

    # Print or return the best parameters and performance
    print("Best parameters:", best_params)
    print("Best performance:", best_performance)
    print(len(results))
    print(len(data))

si_data = si_data[~np.isnan(si_data).any(axis=1) & ~np.isinf(si_data).any(axis=1)]


si_values = np.array([r[5] for r in si_data]) % 180
dt_values = np.array([r[6] for r in si_data])
tilt_values = np.array([r[7] for r in si_data])


# stability = True 

# # if stability: 

# # Best parameter set from previous grid search
# # best_params = (snr_opt, stalta_opt, sksweight_opt, sweight_opt, sdepthweight_opt)

# # Define percentage perturbations (small changes)
# perturbations = np.array([-0.1, -0.05, -0.01, 0.01, 0.05, 0.1])  # ±10%, ±5%, ±1%

# # Generate perturbed parameter sets
# perturbed_param_sets = []
# for i, param in enumerate(best_params):
#     for delta in perturbations:
#         new_params = list(best_params)  # Copy the original best params
#         new_params[i] = param * (1 + delta)  # Apply perturbation
#         perturbed_param_sets.append(tuple(new_params))

# # Run SI_finder for each perturbed parameter set
# results = []
# for params in perturbed_param_sets:
#     si_value = SI_finder(params,eventdataa,station,'SI')  # Compute splitting intensity
#     results.append((params, (si_value)))
# print(results)

# si_values = np.array([r[1][0] for r in results]) % 180
# dt_values = np.array([r[1][1] for r in results])
# tilt_values = np.array([r[1][2] for r in results])

print(si_values)

# plt.hist(tilt_values)
# plt.savefig('../data/sample/'+station+'/SI_sensitivity_apr1.png',dpi=200)# # Compute variability
# # si_values = np.array([r[1] for r in results])
# # mean_si = np.mean(si_values)
# # std_si = np.std(si_values)

# # print(f"Mean Splitting Intensity: {mean_si:.4f}")
# # print(f"Standard Deviation: {std_si:.4f}")

# # # Optional: Plot results
# # import matplotlib.pyplot as plt

# # perturbation_labels = [f"{p}%" for p in (perturbations * 100)]
# # plt.figure(figsize=(8, 5))
# # plt.plot(perturbation_labels, si_values, marker='o', linestyle='-')
# # plt.axhline(mean_si, color='r', linestyle='--', label="Mean SI")
# # plt.xlabel("Parameter Perturbation (%)")
# # plt.ylabel("Splitting Intensity")
# # plt.legend()
# # plt.title("Sensitivity of Splitting Intensity to Parameter Perturbations")
# # plt.savefig('../data/sample/'+station+'/SI_sensitivity_apr1.png',dpi=200)