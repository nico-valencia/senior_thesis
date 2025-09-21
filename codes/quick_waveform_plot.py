import pandas as pd 
import numpy as np 
import glob 
import os 
from obspy import read 
from obspy import geodetics
import splitwavepy as sw 
from support import plotdata,covar_particle_motion, plot_trace_final, envelopes
import sys
import copy 
import matplotlib.pyplot as plt 
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.filter import envelope

from support import SplittingIntensity
import traceback

'''
SELECT STATION, EVENTDATA FILE, OUTPUT DIRECTORY, STATION INFO
'''
station = 'TAM'
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_mar5.txt',sep='|')
output_dir      = '../data/sample/'+station+'/'+station
dir_station     = '../data/sample/'+station+'/'+station+'_stations_info.txt'
station_df  = pd.read_csv(dir_station,sep='|',header=0)
station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]



'''
CHOOSE PHASES TYPE AND CREATE LIST OF WAVEFORMS, CONFIRM CODE IS IN AUTOMATIC MODE
'''
auto = True
# phase = ['S']
phase=['SKS','S']
waveforms = []
for ph in phase:
    # dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'
    dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'



    waveforms.extend(glob.glob(dir+'*'))
bruh = 0 

if auto:
    pass
else:
    sys.exit()

'''
MAIN ANALYSIS LOOP PER WAVEFORM
'''
for i in waveforms:
# for i in range(100,101,1):
    x = read(i)
    

    # DEFINE EVENT DATA NEEDED FOR ANY LATER PROCESS 
    pieces  = i.split('-')
    # pieces  = waveforms[i].split('-')
    bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
    time = bit.split('.')[0]+'.'+bit.split('.')[1]
    evind  = eventdata[eventdata.iloc[:, 0]==time].index
    evtime = eventdata['evtime'].iloc[evind]
    evlat = eventdata['evlat'].iloc[evind]
    evlon = eventdata['evlon'].iloc[evind]
    evdp = eventdata['evdp'].iloc[evind]
    evphase = eventdata['evphase'].iloc[evind]
    evbaz = eventdata['evbaz'].iloc[evind]
    evmag = eventdata['evmg'].iloc[evind]
    evepic = eventdata['evepic'].iloc[evind]
    # print(evdp)

    
   
    try:

        # copied trace to check particle motion later 
        motiontrace = copy.deepcopy(x)

        ### usual code for measuring anisotropy and window creation 

        # filtering parameters
        freqmin = 0.01        # 0.01
        freqmax = 0.125  # 0.125   (8-100s)
        
        # demean, detrend, and bandpass
        x.detrend('demean')
        x.detrend('linear')
        x.filter("bandpass",freqmin=freqmin,freqmax=freqmax)

        north           = x.select(component="N")[0].data
        east            = x.select(component="E")[0].data
        sample_interval = x[0].stats.delta

        print('north and east array lengths: ',len(north),len(east))


        # makes number of samples odd for processing
        if len(north) % 2 == 0:
                north   = north[:-1]
        else:
                pass
        
        if len(east) % 2 == 0: 
                east    = east[:-1]
        else: 
                pass


        ### new window calculated using sta/lta algorithm from obspy 
        n_sta_lta = classic_sta_lta(north,200,800)
        e_sta_lta = classic_sta_lta(east,200,800)

        time = list(range(len(north)))

        sum_sta_lta = n_sta_lta+e_sta_lta
        phase_ind = np.argmax(sum_sta_lta)
        sta_lta_max = sum_sta_lta[phase_ind]

        print(f"STA LTA maximum: {sta_lta_max} at {time[phase_ind]*sample_interval} seconds.")

        windowstart = phase_ind - 200   # play around with this window half width value...
        windowend = phase_ind + 200

        #############################################

        # determine particle motion axes
        particle_freqmin = 0.02 
        particle_freqmax = 0.066   # 15 - 50s
        motiontrace.detrend('demean')
        motiontrace.detrend('linear') 
        motiontrace.filter('bandpass',freqmin=particle_freqmin,freqmax=particle_freqmax)

        motion_north           = motiontrace.select(component="N")[0].data
        motion_east            = motiontrace.select(component="E")[0].data
        motion_sample_interval = motiontrace[0].stats.delta


        # makes number of samples odd for processing
        if len(motion_north) % 2 == 0:
                motion_north   = motion_north[:-1]
        else:
                pass
        
        if len(motion_east) % 2 == 0: 
                motion_east    = motion_east[:-1]
        else: 
                pass
        
        calcbaz,xlam1,xlam2 = covar_particle_motion(motion_north[windowstart:windowend],motion_east[windowstart:windowend])


        ### make sure that calcbaz is within same quadrant as evbaz 

        if abs(calcbaz - evbaz.iloc[0]) > 90:
            calcbaz = (calcbaz+180) % 360

        print('Calculated Polarization Direction: ',calcbaz)

        #############################################

        # envelope filtering for SNR (based on technique used previously

        pair    = sw.Pair(north, east, delta=sample_interval)

        mean_snr , max_snr, ne_env, mean_snr2, max_snr2  = envelopes(pair,float(evbaz.iloc[0]),float(calcbaz),windowstart,windowend)
        print("envelope data:")
        print(mean_snr , mean_snr2)
        print(max_snr , max_snr2)
        print("-----")

        print(evbaz,calcbaz)

        
        eventdata.loc[evind,'windowstart'] = windowstart*sample_interval
        eventdata.loc[evind,'windowend'] = windowend*sample_interval
        eventdata.loc[evind,'calcbaz'] = calcbaz
        eventdata.loc[evind,'sta_lta'] = sta_lta_max
        eventdata.loc[evind,'snrbefore'] = mean_snr
        eventdata.loc[evind,'snrafter'] = mean_snr2
        eventdata.loc[evind,'ellipticity'] = xlam1/xlam2
        
        print('worked')
        bruh+=1

        # fig_string = str('/home/gcl/BR/nicolasv03/deep_eq_plots/'+time+'.png')

        kwargs = {'title': f"Event: {evtime.iloc[0]} Phase: {evphase.iloc[0]} Depth: {evdp.iloc[0]} Mag: {evmag.iloc[0]} Epic: {round(evepic.iloc[0],1)} Ratio: {round(xlam1/xlam2,2)}"}
 
        # if (sta_lta_max > 6) & (mean_snr > 2) & (mean_snr2 > 1.5): # & (abs(calcbaz - evbaz.iloc[0]) < 20):
        #     print('let us test this')
        #     # print(i)
        #     trouble_wave = '1990-05-30T02:34:02.560000Z' + '.mseed'
        #     print(trouble_wave)
        #     print(bit)
        #     # print('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_'+evphase+'_waveforms/G-TAM-1990-03-25T13:22:55.930000Z.mseed')
        #     plot_trace_final(pair,float(evbaz.iloc[0]),float(calcbaz),windowstart,windowend,ne_env,**kwargs)

        # trouble_wave = '1990-05-30T02:34:02.560000Z' + '.mseed'
        # trouble_wave = '2001-07-07T09:38:43.700000Z' + '.mseed'
        # trouble_wave = '1995-05-18T00:06:27.680000Z' + '.mseed'
        # trouble_wave = '2000-05-12T18:43:15.160000Z' + '.mseed'
        # trouble_wave = '2005-06-02T10:56:00.540000Z' + '.mseed'
        trouble_wave = '2005-01-02T15:35:55.920000Z' + '.mseed'

        if trouble_wave == bit:
            print('yo')
            print(trouble_wave)
            print(bit)
            plot_trace_final(pair,float(evbaz.iloc[0]),float(calcbaz),windowstart,windowend,ne_env,**kwargs)



        
    except Exception as exception: 
        print('didnt work lol')
        print(f"Error: {exception}")
        traceback.print_exc()
        bruh += 1