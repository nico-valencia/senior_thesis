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
MANUAL IDENTIFICATION OF WAVE QUALITY 
MANUAL IDENTIFICATION OF NOISE AND EVENT WINDOWS
IDENTIFICATION OF PARTICLE MOTION AXES
SAVES VALUES OF WINDOWS AND WAVEFORM QUALITY
APPLIED TO ALL WAVES

*testing automatic version of this right now* 
'''


'''
SELECT STATION, EVENTDATA FILE, OUTPUT DIRECTORY, STATION INFO
'''
data_tag = 'AU'
station = 'australia'
data_type = 'complete'
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_apr1.txt',sep='|')
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_apr16.txt',sep='|')
output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag
dir_station     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
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
    dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/'

    waveforms.extend(glob.glob(dir+'*'))
bruh = 0 

if auto:
    pass
else:
    sys.exit()

print(station_df.head())
print(eventdata.head())
print(len(eventdata))
print(len(station_df))
# sys.exit()

'''
MAIN ANALYSIS LOOP PER WAVEFORM
'''
for i in waveforms:
# for i in range(100,101,1):
    x = read(i)
    # x = read(waveforms[i])

    # DEFINE EVENT DATA NEEDED FOR ANY LATER PROCESS 
    pieces  = i.split('-')
    # pieces  = waveforms[i].split('-')
    bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
    time = bit.split('.')[0]+'.'+bit.split('.')[1]
    evind  = eventdata[eventdata.iloc[:,'evtime']==time].index
    # evind = eventdata[eventdata['evtime'] == time].index
    evtime = eventdata['evtime'].iloc[evind]
    evlat = eventdata['evlat'].iloc[evind]
    evlon = eventdata['evlon'].iloc[evind]
    evdp = eventdata['evdp'].iloc[evind]
    evphase = eventdata['evphase'].iloc[evind]
    evbaz = eventdata['evbaz'].iloc[evind]
    evmag = eventdata['evmg'].iloc[evind]
    evepic = eventdata['evepic'].iloc[evind]
    # print(evdp)
    # print(x)
    # print(pieces)
    # print(bit)
    # print(time)
    # print(eventdata[eventdata['evtime'] == time])
    # print(evtime,evlat,evlon,evbaz,evmag,evepic)
    sys.exit()

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

        # mar 11 2026 is this correction necessary when unraveling SI plots 0-360? 
        # yes, this is important since the roatation correlation only looks at degrees between 0-180

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
                    
                    
                    
                #   'figname': fig_string
                #   }
        
        # if sta_lta_max > 0:
        # if (max_snr2 or max_snr) < 1:

        # print("HERE IT ISSSSSSSS:")
        # print(SplittingIntensity(pair,float(calcbaz)))
        # print("AND THAT IS ALL")

        # if (sta_lta_max > 6) & (mean_snr > 2) & (mean_snr2 > 1.5): # & (abs(calcbaz - evbaz.iloc[0]) < 20):
        #     plot_trace_final(pair,float(evbaz.iloc[0]),float(calcbaz),windowstart,windowend,ne_env,**kwargs)

        #     print("HERE IT ISSSSSSSS:")
        #     print(SplittingIntensity(pair,float(calcbaz)))
        #     print("AND THAT IS ALL")
        #     # plt.savefig('/home/gcl/BR/nicolasv03/deep_eq_plots/'+time+'.png',dpi=800)
        #     # fig.savefig('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/FITZ/pm_plots'+time+'.png',dpi=800)
        

        # plot_trace(pair,baz,figname='../data/sample/CAN/plots/'+i+'_fr_example_plot.png',title=str(tit))
            # plt.savefig
        sys.exit()
    except Exception as exception: 
        print('didnt work lol')
        print(f"Error: {exception}")
        traceback.print_exc()
        bruh += 1
        sys.exit()

print(len(waveforms))
print(len(eventdata))
print(bruh)


eventdata.to_csv(output_dir+'_events_apr17.txt', sep='|', index=False)


# %----- construct parent window -----% # stealing this from RF synthetic code to identify the peak ends !!
# Par_rf  = P_conv; 
# Dau_rf  = SV_conv;

# Par_rf_win  = timewin(time,Par_rf,-3,3);
# [a,b]       = findpeaks(Par_rf_win);

# for j=1:length(b)                           % identifies peak of parent wave
#    if max(Par_rf_win) == a(j);n=b(j);end
# end

# for k=n:length(Par_rf)                      % identifies end of parent wave
#     if Par_rf(k)<0; npte=k;break;end
#     if Par_rf(k+1)>Par_rf(k); npte=k;break;end
# end

# for k=n:-1:1                                % identifies beginning of parent wave
#     if Par_rf(k)<0; nptb=k;break;end    
#     if Par_rf(k-1)>Par_rf(k); nptb=k;break;end
# end

# Par_rf(1:nptb)      = 0.0;                  % all trace values outside parent wave = 0 
# Par_rf(npte:end)    = 0.0;  
