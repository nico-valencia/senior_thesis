# from support import plot_trace , lambda_filter_2
# import os 
# import glob 
# import numpy as np 
# import pandas as pd 
# from obspy import read
# import splitwavepy as sw 
# from obspy import geodetics
# from splitwavepy.core import data
# from support import SplittingIntensity
# import sys
# import traceback
# import time as ttime

# # import matplotlib
# # matplotlib.use('Agg')  # Use a non-interactive backend

# station = 'TAM'

data_tag = 'TAM'
station = 'TAM'
data_type = 'sample'

# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_apr16.txt',sep='|')
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_windowed_events_may30.txt',sep='|')

output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag
dir_station     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'

# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_apr1.txt',sep='|')

# output_dir      = '../data/sample/'+station+'/'+station

# dir_station     = '../data/sample/'+station+'/'+station+'_stations_info.txt'

station_df  = pd.read_csv(dir_station,sep='|',header=0)
station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]


# results_df = pd.DataFrame(columns=['evtime','evphase','phi','dt','snr','lambda','dfast','dlag'])

eventdata['phi'] =0
eventdata['dt']= 0
eventdata['snr']= 0
eventdata['lambda']= 0
eventdata['dfast']=0
eventdata['dlag']=0

eventdata['wiggleSI_act'] = 0 
eventdata['wiggleSIlow_act'] = 0 
eventdata['wiggleSIhigh_act'] = 0 

eventdata['wiggleSI_calc'] = 0 
eventdata['wiggleSIlow_calc'] = 0 
eventdata['wiggleSIhigh_calc'] = 0 



# print(eventdata)

phase = ['SKS','S']
# phase = ['S']

# station = 'STKA'
waveforms = []
for ph in phase:
    # dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'
    # dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'
    dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/'


    waveforms.extend(glob.glob(dir+'*'))

bruh = 0 
result_ind = 0 

# print('tst')

# for i in waveforms:
#     # print(i)
# # for i in range(30,40,1):
#     x = read(i)
#     # print(x)
#     # x = read(dir+i)
#     pieces  = i.split('-')
#     stat = pieces[1]
#     bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
#     time = bit.split('.')[0]+'.'+bit.split('.')[1]
#     # evind  = eventdata[eventdata.iloc[:, 0]==time].index[0] # this is [0] to turn from index to str 
#     evind = eventdata[(eventdata['evtime'] == time) & (eventdata['evstation'] == stat)].index

#     # FIX THIS BEFORE I RUN THIS CODE BY TOMORROW LMAO
#     evlat = eventdata['evlat'].iloc[evind]
#     evlon = eventdata['evlon'].iloc[evind]
#     evdp = eventdata['evdp'].iloc[evind]
#     evphase = eventdata['evphase'].iloc[evind]
#     evbaz = eventdata['evbaz'].iloc[evind]
#     evwindowstart = eventdata['windowstart'].iloc[evind]
#     evwindowend = eventdata['windowend'].iloc[evind]
#     evcalcbaz = eventdata['calcbaz'].iloc[evind]

#     # print(evphase.to_string())
#     # print(evwindowstart)
#     # print(evwindowend)
#     print(f"station: {stat} and phase: {str(evphase.iloc[0])}")
#     # print('printing evcalcbaz')
#     # print(evcalcbaz)
#     # # print(evind)
#     # print(eventdata[(eventdata['evtime'] == time) & (eventdata['evstation'] == stat)])
#     # sys.exit()
#     # print(evdp)

#     # print(evdp)
#     try:
#     # if 1 ==1 :
#             # filtering parameters
#             freqmin = 0.01        # 0.01
#             freqmax = 0.125  # 0.125   (8-100s)
            
#             # demean, detrend, and bandpass
#             x.detrend('demean')
#             x.detrend('linear')
#             x.filter("bandpass",freqmin=freqmin,freqmax=freqmax)


#             north           = x.select(component="N")[0].data
#             east            = x.select(component="E")[0].data
#             sample_interval = x[0].stats.delta


#             # makes number of samples odd for processing
#             if len(north) % 2 == 0:
#                     north   = north[:-1]
#             else:
#                     pass
            
#             if len(east) % 2 == 0: 
#                     east    = east[:-1]
#             else: 
#                     pass

#             pair    = sw.Pair(north, east, delta=sample_interval)
#             # pair.rotateto(evbaz)

#             # change this to see how window changes things
            
#             pair.set_window(float(evwindowstart.iloc[0]),float(evwindowend.iloc[0]))

#             # calculate beta 
#             print('beta start is here ----------------')
#             print(type(evcalcbaz))
#             print(type(evbaz))
#             print(evcalcbaz)
#             print(evbaz)
#             beta = ((((float(evcalcbaz.iloc[0]) - float(evbaz.iloc[0])) % 180) + 90) % 180) - 90 
#             # eventdata['deltaphi'].iloc[evind] = beta
#             #             eventdata.loc[evind, 'phi'] = fast

#             eventdata.loc[evind,'beta'] = beta

#             # print('evcalcbaz:')
#             # print(evcalcbaz)
#             # print('evbaz')
#             # print(evbaz)
#             # print(beta)
#             print('yo we are here')
#             wiggleSIlow_act , wiggleSI_act , wiggleSIhigh_act = SplittingIntensity(pair,float(evbaz.iloc[0]))
#             wiggleSIlow_calc , wiggleSI_calc , wiggleSIhigh_calc = SplittingIntensity(pair,float(evcalcbaz.iloc[0]))
#             print('yo SI ran')


#             # if str(evphase) == str('S'):
#             #         print('gonna do S phase measurement')
#             #         # measure = sw.EigenM(pair,lags=np.linspace(0,3,61),degs=180)
#             #         measure = sw.TransM(pair, pol=evcalcbaz, lags=np.linspace(0,3,61))
#             #         # obtains fast angle and delay time values (phi and dt)
#             #         fast    = measure.fast # - beta
#             #         lag     = measure.lag


#             # if str(evphase) == str('SKS'):
#             #     print('gonna do SKS phase measurement')
#             #     if abs(beta) <= 100:
#             #         measure = sw.TransM(pair, pol=evcalcbaz, lags=np.linspace(0,3,61))
                    
#             #         # obtains fast angle and delay time values (phi and dt)
#             #         fast    = measure.fast
#             #         lag     = measure.lag

#             # measure = sw.TransM(pair, pol=evcalcbaz, lags=np.linspace(0,3,61))
#             # fast = measure.fast 
#             # lag = measure.lag
#             # # Signal / Noise Filter (Restivo & Helffrich 1999)
#             # snr    = measure.snr() 
#             # snrstr = str(snr)  

#             # print('yo this ran too')
   

            

#             # # matrix and grid search parameters 
#             # filters = lambda_filter_2(measure)
#             # lambda_ratio = str(filters[0])
#             # dfast = str(filters[1])
#             # dlag = str(filters[2])

#             # # Assuming result_ind and evind are defined and valid indices
#             # if result_ind < len(results_df) and evind < len(eventdata):
#             #         results_df.at[result_ind, 'evtime'] = eventdata['evtime'].iloc[evind]
#             #         results_df.at[result_ind, 'evphase'] = evphase
#             #         results_df.at[result_ind, 'phi'] = float(fast)
#             #         results_df.at[result_ind, 'dt'] = float(lag)
#             #         results_df.at[result_ind, 'snr'] = float(snr)
#             #         results_df.at[result_ind, 'lambda'] = float(lambda_ratio)
#             #         results_df.at[result_ind, 'dfast'] = float(dfast)
#             #         results_df.at[result_ind, 'dlag'] = float(dlag)
                    
#             #         result_ind += 1
#             # else:
#             #         print("Index out of bounds")
#             #         print(result_ind)
#             #         print(evind)
#             # results_df['evtime'].loc[result_ind] = eventdata['evtime'].iloc[evind]
#             # results_df['evphase'].loc[result_ind] = evphase
#             # results_df['phi'].loc[result_ind] = float(fast)
#             # results_df['dt'].loc[result_ind] = float(lag)
#             # results_df['snr'].loc[result_ind] = float(snr)
#             # results_df['lambda'].loc[result_ind] = float(lambda_ratio)
#             # results_df['dfast'].loc[result_ind] = float(dfast)
#             # results_df['dlag'].loc[result_ind] = float(dlag)
#             # result_ind += 1

#             # eventdata['phi'].iloc[evind] = fast
#             # eventdata['dt'].iloc[evind] = lag
#             # eventdata['snr'].iloc[evind] = snr

#             # eventdata['wiggleSIlow_act'].iloc[evind] = wiggleSIlow_act
#             # eventdata['wiggleSI_act'].iloc[evind] = wiggleSI_act
#             # eventdata['wiggleSIhigh_act'].iloc[evind] = wiggleSIhigh_act

#             # eventdata['wiggleSIlow_calc'].iloc[evind] = wiggleSIlow_calc
#             # eventdata['wiggleSI_calc'].iloc[evind] = wiggleSI_calc
#             # eventdata['wiggleSIhigh_calc'].iloc[evind] = wiggleSIhigh_calc

#             # eventdata.loc[evind, 'phi'] = fast
#             # eventdata.loc[evind, 'dt'] = lag
#             # eventdata.loc[evind, 'snr'] = snr

#             eventdata.loc[evind, 'wiggleSIlow_act'] = wiggleSIlow_act
#             eventdata.loc[evind, 'wiggleSI_act'] = wiggleSI_act
#             eventdata.loc[evind, 'wiggleSIhigh_act'] = wiggleSIhigh_act

#             eventdata.loc[evind, 'wiggleSIlow_calc'] = wiggleSIlow_calc
#             eventdata.loc[evind, 'wiggleSI_calc'] = wiggleSI_calc
#             eventdata.loc[evind, 'wiggleSIhigh_calc'] = wiggleSIhigh_calc



#             # eventdata['lambda'].iloc[evind] = lambda_ratio
#             # eventdata['dfast'].iloc[evind] = float(dfast)
#             # eventdata['dlag'].iloc[evind] = float(dlag)

            
#             # print('check')
#             # print(fast,lag)

#             # plot_trace(pair,baz,figname='../data/sample/CAN/plots/'+i+'_fr_example_plot.png',title=str(tit))
#             # plt.savefig
#             print('code worked, SI value is:', wiggleSI_calc)
#     except Exception as exception: 
#             print('didnt work lol')
#             print(f"Error: {exception}")
#             traceback.print_exc()
#             # ttime.sleep(5)
#             bruh += 1

# # print(len(waveforms))
# # print(len(eventdata))
# # print(bruh)
# # # print(len(results_df))

# # print(eventdata)

# # for i in range(0,len(eventdata)):
# #         for j in range(0,len(results_df)):
# #                 if eventdata['evtime'].iloc[i] == results_df['evtime'].iloc[j]:
# #                         print('yes')
# eventdata = eventdata.drop_duplicates(subset=['evtime'])
# eventdata.to_csv(output_dir+'_events_measurements_may30.txt', sep='|', index=False)


import os
import glob
import numpy as np
import pandas as pd
from obspy import read
import splitwavepy as sw
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from support import SplittingIntensity

def process_waveform(i, eventdata):
    try:
        x = read(i)
        pieces  = i.split('-')
        stat = pieces[1]
        bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
        time = bit.split('.')[0]+'.'+bit.split('.')[1]

        # Match event
        evrow = eventdata[(eventdata['evtime'] == time) & (eventdata['evstation'] == stat)]
        if evrow.empty:
            return None  # skip if not found
        evind = evrow.index[0]

        evbaz = evrow['evbaz'].iloc[0]
        evcalcbaz = evrow['calcbaz'].iloc[0]
        evwindowstart = evrow['windowstart'].iloc[0]
        evwindowend = evrow['windowend'].iloc[0]

        # Filter
        x.detrend('demean')
        x.detrend('linear')
        x.filter("bandpass", freqmin=0.01, freqmax=0.125)

        north = x.select(component="N")[0].data
        east  = x.select(component="E")[0].data
        sample_interval = x[0].stats.delta

        if len(north) % 2 == 0: north = north[:-1]
        if len(east) % 2 == 0:  east = east[:-1]

        pair = sw.Pair(north, east, delta=sample_interval)
        pair.set_window(float(evwindowstart), float(evwindowend))

        beta = ((((float(evcalcbaz) - float(evbaz)) % 180) + 90) % 180) - 90

        # Splitting intensity
        wiggleSIlow_act , wiggleSI_act , wiggleSIhigh_act = SplittingIntensity(pair,float(evbaz))
        wiggleSIlow_calc, wiggleSI_calc, wiggleSIhigh_calc = SplittingIntensity(pair,float(evcalcbaz))

        return {
            'evtime': time,
            'evstation': stat,
            'beta': beta,
            'wiggleSIlow_act': wiggleSIlow_act,
            'wiggleSI_act': wiggleSI_act,
            'wiggleSIhigh_act': wiggleSIhigh_act,
            'wiggleSIlow_calc': wiggleSIlow_calc,
            'wiggleSI_calc': wiggleSI_calc,
            'wiggleSIhigh_calc': wiggleSIhigh_calc
        }

    except Exception as e:
        print(f"Error processing {i}: {e}")
        traceback.print_exc()
        return None
    

results = []
with ProcessPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(process_waveform, wf, eventdata) for wf in waveforms]
    for fut in as_completed(futures):
        res = fut.result()
        if res:
            results.append(res)

# Merge results back into eventdata
results_df = pd.DataFrame(results)
eventdata = eventdata.merge(results_df, on=['evtime','evstation'], how='left')

# eventdata = eventdata.drop_duplicates(subset=['evtime'])
eventdata.to_csv(output_dir+'_events_measurements_may30.txt', sep='|', index=False)