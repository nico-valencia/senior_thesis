import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy.optimize as sci 
import scipy.stats as stats
import sys
import traceback
import time
# from astropy.stats import sigma_clip

'''
all data analysis techniques looped across the entire network


'''



# station = 'STKA'
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_apr1.txt',sep='|')
# print("Unfiltered measurement catalog length:", len(eventdata))
data_type = 'complete'


data_tag = 'AU'
station = 'australia'
eventdataa = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_measurements_apr21.txt',sep='|')

# data_tag = 'swp'
# station = 'swp'
# eventdataa = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_measurements_may11.txt',sep='|')

# eventdataa = pd.concat([eventdataa,eventdataa2])
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_apr16.txt',sep='|')
dir_stations  = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
stations_df  = pd.read_csv(dir_stations,sep='|',header=0)
output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag


stat = 'CMSA'
station_names = stations_df['Station']
print(station_names)


# output file for saved measurements
columns = ['Station','Latitude','Longitude','phi', 'dt'] 
saved_df = pd.DataFrame(columns=columns)

# sys.exit()


for station_name in station_names:
    print(station_name)

    try:
        if station_name != stat:
            continue

        # eventdata = eventdata[eventdata['evstation'].str.contains('PSA',na=False)]
        eventdata = eventdataa[eventdataa['evstation']==station_name].copy()
        print("Filtered measurement catalog length:", len(eventdata))
        print(stat)
        print()

        if len(eventdata) == 0:
            continue

        # eventdata = eventdata[(eventdata['evdp']>=200) | (eventdata['evphase']=='SKS')]
        '''
        notes from meeting:
        look at deep events and slowly introduce more 
        make much smaller bins for azimuth and maybe ravel measurements
        manually inspect waveforms and see wtf is going on

        '''

        ravel = False

        single = True



        snr_threshold = 4
        snrbefore_threshold = 1.0
        snrafter_threshold = 1.0
        sta_lta_threshold = 6
        ellipticity_threshold = 20
        mag_threshold = 7.0
        sierrmax = 0.3

        s_depth = 200 

        # ### truncate sta, meansnr before and after, ellip

        # # # actually truncate the dataset lmao 
        # eventdata = eventdata[eventdata['snr'] > snr_threshold]
        eventdata = eventdata[eventdata['snrbefore'] > snrbefore_threshold]
        eventdata = eventdata[eventdata['snrafter'] > snrafter_threshold]
        # eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]
        eventdata = eventdata[eventdata['evmg'] < mag_threshold]

        # print(eventdata['snrbefore'] > snrbefore_threshold)
        # print(eventdata[eventdata['snrbefore'] > snrbefore_threshold])

        # eventdata = eventdata[eventdata['ellipticity'] > ellipticity_threshold]

        if ravel: 
            eventdata['calcbaz'] = eventdata['calcbaz'] % 180
            degs = np.arange(0,180,1)
        else:
            degs = np.arange(0,360,1)
        rads = np.radians(degs)



        # event si errors so they can be filtered out

        # eventdata['SI_act_error'] = abs(eventdata['wiggleSI_act']-eventdata['wiggleSIlow_act'])
        # eventdata['SI_calc_error'] = abs(eventdata['wiggleSI_calc']-eventdata['wiggleSIlow_calc'])

        # eventdata.loc[:, 'SI_act_error'] = abs(eventdata['wiggleSI_act'] - eventdata['wiggleSIlow_act'])
        eventdata.loc[:, 'SI_calc_error'] = abs(eventdata['wiggleSI_calc'] - eventdata['wiggleSIlow_calc'])



        # eventdata = eventdata[eventdata['SI_act_error'] < sierrmax]
        eventdata = eventdata[eventdata['SI_calc_error'] < sierrmax]


        # eventdata = eventdata.dropna()[
        # eventdata = eventdata[eventdata['wiggleSI_calc'] != 0]
        eventdata = eventdata[eventdata['SI_calc_error'] != 0]

        # sys.exit()

        ### quality factor criteria (mar 12 2025) ###
        '''
        curve_fit sigma parameter: inverse weight (closer to 0 is better)
        terms to consider: SI error , SNR , phase type , STA/LTA ... 
        '''

        # inv_error = (1 / eventdata['SI_calc_error'])                # smaller error is better (so inverted)
        # snr = eventdata['snrbefore'] + eventdata['snrafter']    # larger value is better 
        # stalta = eventdata['sta_lta']

        print("Filtered measurement catalog length:", len(eventdata))

        if len(eventdata) == 0:
            continue

        # eventdata['phase_weight'] = 0 
        # eventdata['quality'] = 0 
        # for i in range(0,len(eventdata)):
        #     if eventdata['evphase'].iloc[i] == 'SKS':
        #         eventdata['phase_weight'].iloc[i] = 0.46
        #     if eventdata['evphase'].iloc[i] == 'S':
        #         if eventdata['evdp'].iloc[i] >= s_depth:
        #             eventdata['phase_weight'].iloc[i] = 1.0
        #         else:
        #             eventdata['phase_weight'].iloc[i] = 0.1




        # eventdata['quality'] = np.sqrt((inv_error**2)+(snr**2)+(stalta**2)) * eventdata['phase_weight']


        # quality_threshold = np.percentile(eventdata['quality'],99)
        # eventdata = eventdata[eventdata['quality'] <= quality_threshold]


        # eventdata = eventdata.dropna(subset=['quality'])
        # eventdata = eventdata[eventdata['quality']>3]

        # inv_weight = 1 /eventdata['quality']

        # print("Filtered measurement catalog length:", len(eventdata))
        # if len(eventdata) == 0:
        #     print('after filtering for criteria, no data was left')
        #     continue

        ####### april 7 2025 ######
        '''
        gonna attempt to make binned splitting intensity based on baz and then calculate normality or uniformity
        after that then calculate mean/median with a confidence interval on that parameter
        and then these will be used as weights for the model fit 
        '''

        eventdata = eventdata.sort_values(by='calcbaz', ascending=True)

        # Define bins and categorize data
        bin_width = 10
        bin_minimum = 10
        array = np.arange(0, 360, bin_width)
        eventdata['chunk'] = np.digitize(eventdata['calcbaz'], array, right=True)

        # Initialize an empty dictionary to store test results
        normality_results = {}
        uniformity_results = {}

        binned_sks_results = np.zeros((len(array),5))

        # Loop over each bin (chunks)
        for chunk in range(0, len(array)):
            # Extract SI values for the current bin
            all_si_values_in_bin = eventdata[eventdata['chunk'] == chunk+1]
            # sks_phases = all_si_values_in_bin['evphase'] == 'SKS'
            # weights = np.where(sks_phases,1.0,0.5)
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']
            # print(weights)
            ### turned off for batch L
            si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evphase'] == 'SKS']['wiggleSI_calc']
            #### turn back on for regular processing
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']

            # print(f"Azimuth Bin {chunk+1} ({array[chunk%len(array)]}° to {array[(chunk+1)%len(array)]}°):")
            # print(f"Numberof waveforms: {len(si_values_in_bin)}")

            ### calculate sample mean/median and confidence interval assuming normal distribution ###


            n = len(si_values_in_bin)

            if n <= bin_minimum:
                binned_sks_results[chunk,:] = chunk,0,0,0,0
                continue

            # mean = np.mean(si_values_in_bin)
            mean = np.median(si_values_in_bin)
            # mean=np.average(si_values_in_bin,weights=weights)
            std = np.std(si_values_in_bin, ddof=1)
            # variance = np.average((si_values_in_bin - mean) ** 2, weights=weights)
            # std = np.sqrt(variance)
            conf = 0.95

            # Get t-critical value
            t_crit = stats.t.ppf((1 + conf) / 2, df=n-1)

            # Margin of error
            margin = t_crit * (std / np.sqrt(n))

            # Confidence interval
            ci_lower = mean - margin
            ci_upper = mean + margin

            # weight = 1 / ((2*margin)**2)
            weight = margin

            # print(f"margin of error for SI: {margin}")
            # print(weight)

            # print(f"95% Confidence Interval for the Mean: ({ci_lower:.3f}, {ci_upper:.3f})")

            # if n <= 1:
            #     binned_results[chunk,:] = chunk,mean,0,array[chunk]+bin_width/2,margin
            # else:
            #     binned_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin

            binned_sks_results[chunk,:] = chunk,mean,0.05*weight,array[chunk]+bin_width/2,margin 


            # print('-' * 50)

        binned_s_results = np.zeros((len(array),5))

        # Loop over each bin (chunks)
        for chunk in range(0, len(array)):
            # Extract SI values for the current bin
            all_si_values_in_bin = eventdata[eventdata['chunk'] == chunk+1]
            # sks_phases = all_si_values_in_bin['evphase'] == 'SKS'
            # weights = np.where(sks_phases,1.0,0.5)
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']
            # print(weights)
            all_si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evdp'] < s_depth]
            si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evphase'] == 'S']['wiggleSI_calc']
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']

            # print(f"Azimuth Bin {chunk+1} ({array[chunk%len(array)]}° to {array[(chunk+1)%len(array)]}°):")
            # print(f"Numberof waveforms: {len(si_values_in_bin)}")

            ### calculate sample mean/median and confidence interval assuming normal distribution ###


            n = len(si_values_in_bin)

            if n <= bin_minimum:
                binned_s_results[chunk,:] = chunk,0,0,0,0
                continue

            # mean = np.mean(si_values_in_bin)
            mean = np.median(si_values_in_bin)
            # mean=np.average(si_values_in_bin,weights=weights)
            std = np.std(si_values_in_bin, ddof=1)
            # variance = np.average((si_values_in_bin - mean) ** 2, weights=weights)
            # std = np.sqrt(variance)
            conf = 0.95

            # Get t-critical value
            t_crit = stats.t.ppf((1 + conf) / 2, df=n-1)

            # Margin of error
            margin = t_crit * (std / np.sqrt(n))

            # Confidence interval
            ci_lower = mean - margin
            ci_upper = mean + margin

            # weight = 1 / ((2*margin)**2)
            weight = margin

            # print(f"margin of error for SI: {margin}")
            # print(weight)

            # print(f"95% Confidence Interval for the Mean: ({ci_lower:.3f}, {ci_upper:.3f})")

            # if n <= 1:
            #     binned_results[chunk,:] = chunk,mean,0,array[chunk]+bin_width/2,margin
            # else:
            #     binned_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin

            binned_s_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin 


            # print('-' * 50)

        binned_s_depth_results = np.zeros((len(array),5))

        # Loop over each bin (chunks)
        for chunk in range(0, len(array)):
            # Extract SI values for the current bin
            all_si_values_in_bin = eventdata[eventdata['chunk'] == chunk+1]
            # sks_phases = all_si_values_in_bin['evphase'] == 'SKS'
            # weights = np.where(sks_phases,1.0,0.5)
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']
            # print(weights)
            all_si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evdp'] >= s_depth]
            si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evphase'] == 'S']['wiggleSI_calc']
            # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']

            # print(f"Azimuth Bin {chunk+1} ({array[chunk%len(array)]}° to {array[(chunk+1)%len(array)]}°):")
            # print(f"Numberof waveforms: {len(si_values_in_bin)}")

            ### calculate sample mean/median and confidence interval assuming normal distribution ###


            n = len(si_values_in_bin)

            if n <= bin_minimum:
                binned_s_depth_results[chunk,:] = chunk,0,0,0,0
                continue

            # mean = np.mean(si_values_in_bin)
            mean = np.median(si_values_in_bin)
            # mean=np.average(si_values_in_bin,weights=weights)
            std = np.std(si_values_in_bin, ddof=1)
            # variance = np.average((si_values_in_bin - mean) ** 2, weights=weights)
            # std = np.sqrt(variance)
            conf = 0.95

            # Get t-critical value
            t_crit = stats.t.ppf((1 + conf) / 2, df=n-1)

            # Margin of error
            margin = t_crit * (std / np.sqrt(n))

            # Confidence interval
            ci_lower = mean - margin
            ci_upper = mean + margin

            # weight = 1 / ((2*margin)**2)
            weight = margin

            # print(f"margin of error for SI: {margin}")
            # print(weight)

            # print(f"95% Confidence Interval for the Mean: ({ci_lower:.3f}, {ci_upper:.3f})")

            # if n <= 1:
            #     binned_results[chunk,:] = chunk,mean,0,array[chunk]+bin_width/2,margin
            # else:
            #     binned_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin

            binned_s_depth_results[chunk,:] = chunk,mean,0.25*weight,array[chunk]+bin_width/2,margin 


            # print('-' * 50)


        binned_s_results = binned_s_results[binned_s_results[:,2] != 0, :]
        binned_sks_results = binned_sks_results[binned_sks_results[:,2] != 0, :]
        binned_s_depth_results = binned_s_depth_results[binned_s_depth_results[:,2] != 0, :]

        # turned off for now for L testing
        binned_results  = np.vstack((binned_s_results, binned_sks_results,binned_s_depth_results))
        binned_truncated_results = np.vstack((binned_sks_results,binned_s_depth_results))
        # binned_truncated_results = np.vstack((binned_sks_results))

        # binned_results = binned_sks_results




        print(binned_results)
        print(array)

        # binned_results = binned_results[np.where(binned_results[:,2]==0),:]

        # sys.exit()

    





        plt.figure(figsize=(30,10))
        # plt.ylim(-2,2)









        print('testing binned results')
        print(binned_results[:,3])
        print(binned_results[:,1])
        print(binned_results[:,2])





        # ####################
        if single:
            def SI_function(baz,fast,lag):
                SI = lag * 2 * np.sin(np.deg2rad(baz-fast)) * np.cos(np.deg2rad(baz-fast))
                return SI
            bounds = ([-360, 0], [360,3])
            # p0 = [90.0 , 1.0]

            # raw measurements 
            # sine_fit = sci.curve_fit(SI_function,eventdata['calcbaz'],eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
            sine_fit = sci.curve_fit(SI_function,binned_results[:,3],binned_results[:,1],sigma=binned_results[:,2],bounds=bounds,full_output=True,nan_policy='omit') #,p0=p0)

            # sine_fit = sci.curve_fit(SI_function,binned_results[:,3],binned_results[:,1],sigma=binned_results[:,2],bounds=bounds,full_output=True,nan_policy='omit') #,p0=p0)
            params,cov,fitdata,trash1,trash2 = sine_fit
            phi,dt = params[0],params[1]

            # sine_fit_2 = sci.curve_fit(SI_function,binned_truncated_results[:,3],binned_truncated_results[:,1],sigma=binned_truncated_results[:,2],bounds=bounds,full_output=True,nan_policy='omit') #,p0=p0)
            # params_2,cov_2,fitdata_2,trash1_2,trash2_2 = sine_fit_2
            # phi_2,dt_2 = params_2[0],params_2[1]



        ### USE MULTI DIMENSIONAL INPUT ARRAY FOR VARIABLE TO GET DATA IN 
        else:
            # x = eventdata['calcbaz']
            # xx = eventdata['evincident']
            # X = ...


            def SI_function(baz_incident,fast,lag,tilt):
                baz,incident = baz_incident
                SI = lag * ((np.cos(np.deg2rad(tilt))**2 * np.sin(2*np.deg2rad(baz-fast))) 
                            + (np.sin(2*np.deg2rad(tilt)) * np.sin(np.deg2rad(baz-fast)) * np.tan(np.deg2rad(incident))))
                return SI
            bounds = ([-360, 0, 0], [360, 3, 45])

            # raw measurements 
            # baz_and_incident = np.array(np.meshgrid(eventdata['calcbaz'], eventdata['evincident']))
            # baz_and_incident = np.vstack((eventdata['calcbaz'], eventdata['evincident']))
            baz_and_incident = np.vstack((binned_results[:,3],np.full_like(binned_results[:,3],20)))
            # baz_and_incident_2 = np.vstack((binned_truncated_results[:,3],np.full_like(binned_truncated_results[:,3],20)))



            # print("baz_and_incident shape:", np.shape(baz_and_incident))  # Should be (2, N)
            # print("eventdata['wiggleSI_calc'] shape:", np.shape(eventdata['wiggleSI_calc']))  # Should be (N,)
            # print("inv_weight shape:", np.shape(inv_weight))  # Should match `eventdata['wiggleSI_calc']`

            sine_fit = sci.curve_fit(SI_function,baz_and_incident,binned_results[:,1],sigma=binned_results[:,2],bounds=bounds,full_output=True,nan_policy='omit') #,p0=p0)
            # sine_fit = sci.curve_fit(SI_function,baz_and_incident,eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
            params,cov,fitdata,trash1,trash2 = sine_fit
            phi,dt,tilt = params[0],params[1],params[2]

            # sine_fit_2 = sci.curve_fit(SI_function,baz_and_incident_2,binned_truncated_results[:,1],sigma=binned_truncated_results[:,2],bounds=bounds,full_output=True,nan_policy='omit') #,p0=p0)
            # params_2,cov_2,fitdata_2,trash1_2,trash2_2 = sine_fit_2
            # phi_2,dt_2,tilt_2 = params_2[0],params_2[1],params_2[2]


        # # # raw measurements 
        # # sine_fit = sci.curve_fit(SI_function,eventdata['calcbaz'],eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
        # # params,cov,fitdata,trash1,trash2 = sine_fit
        # # phi,dt = params[0],params[1]

        sks_tick = 0 
        s_depth_tick = 0 
        s_tick = 0 

        y = (abs(eventdata['evbaz'] - eventdata['calcbaz']) + 180) % 180 #-eventdata['calcbaz']

        # # if abs(calcbaz - baz) > 90:
        # #             calcbaz = (calcbaz+180) % 360
                    
        # colors = plt.cm.plasma((y - y.min()) / (y.max() - y.min()))



        # plt.figure(figsize=(10,6))
        # plt.ylim(-2,2)
        # plt.title(station +' measurements')
        # plt.xlabel('Azimuth')
        # plt.ylabel("SI")
        plt.subplot(1,2,1)

        for i in range(0,len(eventdata)):
            if eventdata['evphase'].iloc[i] == 'SKS':
                label = 'SKS phase' if sks_tick == 0 else None
                sks_tick += 1
                plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='blue',alpha=0.5,label=label)
                plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='blue',barsabove=True,alpha=0.5)
                print('sks_tick:  ', sks_tick)
            
            if eventdata['evphase'].iloc[i] == 'S': #and abs(eventdata['evincident'].iloc[i] - 27) < 5:
                if eventdata['evdp'].iloc[i] >= s_depth:
                    s_depth_tick += 1
                    plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='purple',alpha=0.5)
                    plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='purple',barsabove=True,alpha=0.5)
                    print('s_depth_tick:  ', s_depth_tick)

                else:
                    label = 'S phase' if s_tick == 0 else None
                    s_tick += 1
                    plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='red',alpha=0.5,label=label)
                    plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='red',barsabove=True,alpha=0.5)
                    print('s_tick:  ', s_tick)

                    # plot=plt.scatter(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],marker='.',s=100,c=y.iloc[i], alpha=1, cmap='plasma',vmin=0, vmax=50)
        # plt.colorbar()
        if single: 
            plt.plot(degs , SI_function(degs, *params),color='black',linewidth=5)          #parameters calculated from binned back azimuths 
            # plt.plot(degs , SI_function(degs, *params_2),color='brown',linewidth=5)         # binned but for sks and s depth only 

        else:
            incs = np.full_like(degs,30)
            vars = np.vstack((degs,incs))
            plt.plot(degs , SI_function(vars, *params),color='black',linewidth=5)          #parameters calculated from binned back azimuths 
            # plt.plot(degs , SI_function(vars, *params_2),color='brown',linewidth=5)         # binned but for sks and s depth only 

        # if single:
        #     plt.title(f"{station_name} SI measurements for SKS and S, phi={round(phi,3)} dt={round(dt,3)} ")
        # else:
        #     plt.title(f"{station_name} SI measurements for SKS and S, phi={round(phi,3)} dt={round(dt,3)} tilt={round(tilt,3)} ")


        # plt.plot(degs , SI_function(degs, *params),color='green')          #parameters calculated from binned back azimuths 
        # plt.plot(degs , SI_function(degs, *params) + fitdata['fvec'])

        # plt.colorbar(plot)
        # plt.savefig('../data/sample/'+station+'/SI_measurements_apr7.png',dpi=200)
        # plt.plot(binned_results[:,3],binned_results[:,1],'o',markersize=20,color='orange',markeredgecolor='black')
        # plt.errorbar(binned_results[:,3],binned_results[:,1],yerr=binned_results[:,4],linestyle='',xerr=None,ecolor='orange',markeredgecolor='black',barsabove=True)

        plt.ylim(-2.0,2.0)

        # plt.subplot(1,2,2)

        # # for i in range(0,len(eventdata)):
        # #     if eventdata['evphase'].iloc[i] == 'SKS':
        # #         sks_tick += 1
        # #         plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='blue',alpha=0.5)
        # #         plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='blue',barsabove=True,alpha=0.5)
        # #         print('sks_tick:  ', sks_tick)
            
        # #     if eventdata['evphase'].iloc[i] == 'S': #and abs(eventdata['evincident'].iloc[i] - 27) < 5:
        # #         if eventdata['evdp'].iloc[i] >= s_depth:
        # #             s_depth_tick += 1
        # #             plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='purple',alpha=0.5)
        # #             plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='purple',barsabove=True,alpha=0.5)
        # #             print('s_depth_tick:  ', s_depth_tick)

        # #         else:
        # #             s_tick += 1
        # #             plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='red',alpha=0.5)
        # #             plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='red',barsabove=True,alpha=0.5)
        # #             print('s_tick:  ', s_tick)

        #             # plot=plt.scatter(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],marker='.',s=100,c=y.iloc[i], alpha=1, cmap='plasma',vmin=0, vmax=50)
        # # plt.colorbar()
        # if single: 
        #     plt.plot(degs , SI_function(degs, *params),color='black',linewidth=5)          #parameters calculated from binned back azimuths 
        #     # plt.plot(degs , SI_function(degs, *params_2),color='brown',linewidth=5)         # binned but for sks and s depth only 

        # else:
        #     incs = np.full_like(degs,30)
        #     vars = np.vstack((degs,incs))
        #     plt.plot(degs , SI_function(vars, *params),color='black',linewidth=5)          #parameters calculated from binned back azimuths 
        #     # plt.plot(degs , SI_function(vars, *params_2),color='brown',linewidth=5)         # binned but for sks and s depth only 


        # # plt.plot(degs , SI_function(degs, *params),color='green')          #parameters calculated from binned back azimuths 
        # # plt.plot(degs , SI_function(degs, *params) + fitdata['fvec'])

        # # plt.colorbar(plot)
        # # plt.savefig('../data/sample/'+station+'/SI_measurements_apr7.png',dpi=200)
        # # plt.plot(binned_s_depth_results[:,3],binned_s_depth_results[:,1],markersize=20,color='orange',markeredgecolor='black',marker='o',linestyle='')
        # # plt.plot(binned_results[:,3],binned_results[:,1],markersize=20,color='orange',markeredgecolor='black',marker='x',linestyle='')
        # plt.plot(binned_sks_results[:,3],binned_sks_results[:,1],markersize=10,color='blue',markeredgecolor='blue',marker='s',linestyle='')
        # plt.errorbar(binned_sks_results[:,3],binned_sks_results[:,1],yerr=binned_sks_results[:,4],linestyle='',xerr=None,ecolor='blue',markeredgecolor='blue',barsabove=True)

        # plt.plot(binned_s_results[:,3],binned_s_results[:,1],markersize=10,color='red',markeredgecolor='red',marker='s',linestyle='')
        # plt.errorbar(binned_s_results[:,3],binned_s_results[:,1],yerr=binned_s_results[:,4],linestyle='',xerr=None,ecolor='red',markeredgecolor='red',barsabove=True)

        # plt.plot(binned_s_depth_results[:,3],binned_s_depth_results[:,1],markersize=10,color='purple',markeredgecolor='purple',marker='s',linestyle='')
        # plt.errorbar(binned_s_depth_results[:,3],binned_s_depth_results[:,1],yerr=binned_s_depth_results[:,4],linestyle='',xerr=None,ecolor='purple',markeredgecolor='purple',barsabove=True)



        # # # plt.errorbar(binned_results[:,3],binned_results[:,1],yerr=binned_results[:,4],linestyle='',xerr=None,ecolor='grey',markeredgecolor='black',barsabove=True)

        # plt.ylim(-2.0,2.0)

        # if single:
        #     plt.title(f"{station_name} SI measurements for SKS and S, phi={round(phi_2,3)} dt={round(dt_2,3)} ")
        # else:
        #     plt.title(f"{station_name} SI measurements for SKS and S, phi={round(phi_2,3)} dt={round(dt_2,3)} tilt={round(tilt_2,3)} ")

        # plt.savefig('../data/complete/'+station+'/SI_'+station_name+'_binned_measurements_apr26_2.png',dpi=200)
        # plt.show()
        # plt.savefig('../data/complete/'+station+'/L_may_13_'+station_name+'_binned_SI_measurements.png',dpi=200,bbox_inches='tight')
        # plt.savefig('../data/complete/'+station+'/L_may_13_diagnostic_binned_SI_measurements.png',dpi=200,bbox_inches='tight')





        # Ask for approval via terminal input
        user_input = input(f"Do you want to save this measurement ({station_name} phi = {phi}, dt = {dt})? (y/n): ").strip().lower()

        if user_input == 'y':
            print('saving measurement...')
            # Save the measurement if user approves
            new_row = pd.DataFrame({'Station':[station_name],'Latitude':[eventdata['evstationlat'].iloc[0]],'Longitude':[eventdata['evstationlon'].iloc[0]],'phi': [phi], 'dt': [dt]})
            saved_df = pd.concat([saved_df, new_row], ignore_index=True)
        else:
            print(f"Measurement {i + 1} not saved.")



        #### plot for presentation
        plt.xlabel("Polarization Azimuth",fontsize=20)
        plt.ylabel("Splitting Intensity",fontsize=20)
        plt.title('Station '+station_name+' Splitting Intensity Results: SKS and S phases',fontsize=30)
        plt.legend(fontsize=20)
        print(phi,dt)

        plt.savefig('../data/complete/'+station+'/SI_'+station_name+'_presentation_measurements_may15.png',dpi=200,bbox_inches='tight')
        # #####

        # plt.close()

        # plt.figure(figsize=(6,6))
        # plt.title(f"{station_name} direct phi/dt measurements")
        # plt.xlim(-90,90)
        # plt.ylim(0,3)

        # for i in range(0,len(eventdata)):
        #     if eventdata['evphase'].iloc[i] == 'SKS':
        #         plt.plot(eventdata['phi'].iloc[i],eventdata['dt'].iloc[i],'.',color='blue',alpha=0.5)
            
        #     if eventdata['evphase'].iloc[i] == 'S': #and abs(eventdata['evincident'].iloc[i] - 27) < 5:
        #         if eventdata['evdp'].iloc[i] >= s_depth:
        #             plt.plot(eventdata['phi'].iloc[i],eventdata['dt'].iloc[i],'.',color='purple',alpha=0.5)

        #         else:
        #             plt.plot(eventdata['phi'].iloc[i],eventdata['dt'].iloc[i],'.',color='red',alpha=0.5)


        # plt.savefig('../data/complete/'+station+'/SI_'+station_name+'_direct_measurements_may1.png',dpi=200,bbox_inches='tight')






        print('successfully saved plot for station' +station_name)



        # plt.figure(figsize=(8,8))
        # plt.title(station_name + ' spread of quality')
        # plt.hist(eventdata['quality'],bins=20)
        # # plt.hist(fitdata['fvec'],bins=20)
        # plt.ylabel("Number")
        # plt.xlabel("quality")
        # plt.savefig('../data/complete/'+station+'/quality_spread_apr21.png',dpi=200)
    
    except Exception as exception: 
            print('didnt work lol')
            print(f"Error: {exception}")
            traceback.print_exc()
            # time.sleep(5)


saved_df.to_csv(output_dir+'_L_saved_SI_measurments_may13.txt', sep='|', index=False)


# # print(params)
# print(f"Total number of waveforms used: {len(eventdata)}")
# print(f"SKS waveforms used: {sks_tick}")
# print(f"S depth waveforms used: {s_depth_tick}")
# print(f"S waveforms used: {s_tick}")
# print(f"Fast Angle: {round(phi,4)} Delay Time {round(dt,4)}")
# if single != True:
#     print(f"Tilt Angle: {round(tilt,4)}")

# print("Covariance Matrix:")
# print(cov)

# print(len(eventdata))
# print(len(fitdata['fvec']))

# import pygmt

# catalog = eventdata
# km = 111 * 2  #111km per degree

# circle1 = 60 * km  # degree times deg of radius times 2 for diameter
# circle2 = 90 * km 
# circle3 = 40 * km
# circles = np.zeros(len(catalog))
# circles[:] = circle3
# # radius_km = 200

# # ### station: MEEK
# # center_lat = -26.638 
# # center_lon = 118.614998

# ### station: HYB
# center_lat = 17.41867
# center_lon = 78.55213

# # ### station: CAN
# # center_lat = -35.318715
# # center_lon = 148.996325

# # ### station: FITZ
# # center_lat = -18.0982
# # center_lon = 125.640297




# # Define the figure
# fig = pygmt.Figure()

# # Set region (xmin, xmax, ymin, ymax)
# region = [0, 300, -70, 70]

# # Define stereographic projection (S) with the center at (longitude=0, latitude=90, scale=15 cm)
# # projection = "S130/-13/15c"
# # projection = "N130/15c"
# projection = "M130/-20/15c"


# # Create a stereographic plot
# fig.coast(region=region, projection=projection, shorelines=True, frame="afg", land="gray", water="skyblue")

# # Optional: Add points to the plot
# # fig.plot(x=[0], y=[90], style="c0.5c", color="red", pen="1p,black")
# # fig.plot(x=[0],y=[45])
# fig.plot(x=catalog['evlon'],y=catalog['evlat'], style="c0.3c", fill="white", pen="black")
# # fig.plot(x=[center_lon],y=[center_lat], style="c0.3c", fill="red", pen="black")

# # fig.plot(x=[130], y=[-10], size=[90*km], style=f"SE-{radius_km1}k", pen="1.5p,purple2")
# # fig.plot(x=[130], y=[-10], size=[120*km], style=f"SE-{radius_km2}k", pen="1.5p,purple2")
# # fig.plot(x=[130], y=[-10], style=f"SE-{radius_km}k", pen="2p,red")
# # fig.plot(x=[center_lon], y=[center_lat], size=[circle1], style="E-", pen="1.5p,purple2")
# # fig.plot(x=[center_lon], y=[center_lat], size=[circle2], style="E-", pen="1.5p,purple2")

# # fig.plot(x=catalog['longitude'], y=catalog['latitude'], size=circles, style="E-", pen="1.5p,purple2",transparency=90)
# # fig.plot(x=catalog['longitude'], y=catalog['latitude'], size=circles, style="E-", pen="1.5p,purple2",transparency=98,fill='purple')





# # Save the figure
# # fig.savefig('./deep500km_2.png')
# fig.savefig('../data/sample/'+station+'/eq_map_apr7.png',dpi=200)




# really good loop for weighted CI in one loop one CI per azimuthal bin 

# # Loop over each bin (chunks)
# for chunk in range(0, len(array)):
#     # Extract SI values for the current bin
#     all_si_values_in_bin = eventdata[eventdata['chunk'] == chunk+1]
#     sks_phases = all_si_values_in_bin['evphase'] == 'SKS'
#     weights = np.where(sks_phases,1.0,0.5)
#     si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']
#     print(weights)
#     # si_values_in_bin = all_si_values_in_bin[all_si_values_in_bin['evphase'] == 'S']['wiggleSI_calc']
#     # si_values_in_bin = all_si_values_in_bin['wiggleSI_calc']

#     print(f"Azimuth Bin {chunk+1} ({array[chunk%len(array)]}° to {array[(chunk+1)%len(array)]}°):")
#     print(f"Numberof waveforms: {len(si_values_in_bin)}")

#     ### calculate sample mean/median and confidence interval assuming normal distribution ###


#     n = len(si_values_in_bin)

#     if n <= 1:
#         binned_sks_results[chunk,:] = chunk,0,0,0,0
#         continue

#     # mean = np.mean(si_values_in_bin)
#     mean=np.average(si_values_in_bin,weights=weights)
#     # std = np.std(si_values_in_bin, ddof=1)
#     variance = np.average((si_values_in_bin - mean) ** 2, weights=weights)
#     std = np.sqrt(variance)
#     conf = 0.95

#     # Get t-critical value
#     t_crit = stats.t.ppf((1 + conf) / 2, df=n-1)

#     # Margin of error
#     margin = t_crit * (std / np.sqrt(n))

#     # Confidence interval
#     ci_lower = mean - margin
#     ci_upper = mean + margin

#     # weight = 1 / ((2*margin)**2)
#     weight = margin**2

#     print(f"margin of error for SI: {margin}")
#     print(weight)

#     print(f"95% Confidence Interval for the Mean: ({ci_lower:.3f}, {ci_upper:.3f})")

#     # if n <= 1:
#     #     binned_results[chunk,:] = chunk,mean,0,array[chunk]+bin_width/2,margin
#     # else:
#     #     binned_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin

#     binned_sks_results[chunk,:] = chunk,mean,weight,array[chunk]+bin_width/2,margin 


#     print('-' * 50)