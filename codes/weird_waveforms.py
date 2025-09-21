import pandas as pd 
import numpy as np 
import glob 
import os
import shutil
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
identify weird waveforms that deviate from expected SI results, look at them here!! 
uses plot trace final to look at waveform in different rotations
(instrument, back azimuth, polarization azimuth)
nico valencia may 2025 post grad 
'''

# data_tag = 'swp'
# station = 'swp'
# data_type = 'complete'
# # eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_apr1.txt',sep='|')
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_may9.txt',sep='|')
# output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag
# dir_station     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
# station_df  = pd.read_csv(dir_station,sep='|',header=0)
# station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]

data_type = 'sample'
data_tag = 'TAM'
station = 'TAM'
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_measurements_may30.txt',sep='|')

# data_tag = 'swp'
# station = 'swp'
# eventdataa = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_measurements_may11.txt',sep='|')

# eventdataa = pd.concat([eventdataa,eventdataa2])
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_events_apr16.txt',sep='|')
dir_stations  = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
stations_df  = pd.read_csv(dir_stations,sep='|',header=0)
output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag



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


snr_threshold = 4
snrbefore_threshold = 1.0
snrafter_threshold = 1.0
sta_lta_threshold = 6
ellipticity_threshold = 20
mag_threshold = 7.0
sierrmax = 0.5

s_depth = 200 

# ### truncate sta, meansnr before and after, ellip

# # # actually truncate the dataset lmao 
# eventdata = eventdata[eventdata['snr'] > snr_threshold]
eventdata = eventdata[eventdata['snrbefore'] > snrbefore_threshold]
eventdata = eventdata[eventdata['snrafter'] > snrafter_threshold]
# eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]
eventdata = eventdata[eventdata['evmg'] < mag_threshold]

station_name = 'TAM'
eventdata = eventdata[eventdata['evstation']==station_name].copy()

eventdata.loc[:, 'SI_calc_error'] = abs(eventdata['wiggleSI_calc'] - eventdata['wiggleSIlow_calc'])

# eventdata = eventdata[abs(eventdata['wiggleSI_calc'])<2.0]

# eventdata = eventdata[eventdata['SI_act_error'] < sierrmax]
eventdata = eventdata[eventdata['SI_calc_error'] < sierrmax]


# eventdata = eventdata.dropna()[
# eventdata = eventdata[eventdata['wiggleSI_calc'] != 0]
eventdata = eventdata[eventdata['SI_calc_error'] != 0]
eventdata = eventdata[eventdata['evphase']=='S']

# eventdata = eventdata[eventdata['evdp']>=s_depth]
eventdata = eventdata[eventdata['wiggleSI_calc']<=0]
eventdata = eventdata[(eventdata['calcbaz']>=0)&(eventdata['calcbaz']<=50)]

# 70 - 100 polarization azimuth for CMSA deep events

print(stations_df.head())
print(eventdata.head())
print(len(eventdata))
print(len(stations_df))


filtered_waveforms = [wf for wf in waveforms if "G-TAM" in wf]

# print(filtered_waveforms)

# plt.plot(eventdata['calcbaz'],eventdata['wiggleSI_calc'],'o')
# plt.show()
# sys.exit()

'''
MAIN ANALYSIS LOOP PER EVENT (NOT WAVEFORM)
'''


# for i in range(len(eventdata)):











for i in range(len(filtered_waveforms)):
# for i in range(100,101,1):
    # x = read(i)
    x = read(filtered_waveforms[i])
    # print(x)
    # print(filtered_waveforms[i])

    # # DEFINE EVENT DATA NEEDED FOR ANY LATER PROCESS 
    # pieces  = filtered_waveforms[i].split('-')
    # stat = pieces[1]
    pieces  = filtered_waveforms[i].split('-')
    stat = pieces[1]
    bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
    time = bit.split('.')[0]+'.'+bit.split('.')[1]
    # # evind  = eventdata[eventdata.iloc[:,'evtime']==time].index
    # print('this it the time and station')
    # print(time,stat)
    evind = eventdata[(eventdata['evtime'] == time) & (eventdata['evstation'] == stat)].index
    print(bit)
    # continue
    if len(evind) == 0:
        continue
    # evtime = eventdata['evtime'].iloc[evind]
    # evlat = eventdata['evlat'].iloc[evind]
    # evlon = eventdata['evlon'].iloc[evind]
    # evdp = eventdata['evdp'].iloc[evind]
    # evphase = eventdata['evphase'].iloc[evind]
    # evbaz = eventdata['evbaz'].iloc[evind]
    # evmag = eventdata['evmg'].iloc[evind]
    # evepic = eventdata['evepic'].iloc[evind]
    # # print(evdp)
    # print(x)
    # print(i)
    # print(stat)
    # print(pieces)
    # print(bit)
    # print(time)
    # print(eventdata[(eventdata['evtime'] == time) & (eventdata['evstation'] == stat)])
    # print(evtime,evlat,evlon,evbaz,evmag,evepic)
    # continue

    try:


        # if stat == 'CAN':
        #       continue
        # if stat == 'NOUC':
        #       continue
        # copied trace to check particle motion later 
        # motiontrace = copy.deepcopy(x)

        ### usual code for measuring anisotropy and window creation 

        # filtering parameters
        freqmin = 0.01        # 0.01
        freqmax = 0.125  # 0.125   (8-100s)
        
        # demean, detrend, and bandpass
        x.detrend('demean')
        x.detrend('linear')
        x.filter("bandpass",freqmin=freqmin,freqmax=freqmax)

        north           = x.select(component="N",channel='BH*')[0].data
        east            = x.select(component="E",channel='BH*')[0].data
        sample_interval = x[0].stats.delta

        # print('north and east array lengths: ',len(north),len(east))
        # print(stat)


        # makes number of samples odd for processing
        if len(north) % 2 == 0:
                north   = north[:-1]
        else:
                pass
        
        if len(east) % 2 == 0: 
                east    = east[:-1]
        else: 
                pass


        # ### new window calculated using sta/lta algorithm from obspy 
        # n_sta_lta = classic_sta_lta(north,200,800)
        # e_sta_lta = classic_sta_lta(east,200,800)

        # time = list(range(len(north)))

        # sum_sta_lta = n_sta_lta+e_sta_lta
        # phase_ind = np.argmax(sum_sta_lta)
        # sta_lta_max = sum_sta_lta[phase_ind]

        # print(f"STA LTA maximum: {sta_lta_max} at {time[phase_ind]*sample_interval} seconds.")

        # windowstart = phase_ind - 250   # play around with this window half width value...
        # windowend = phase_ind + 250

        # #############################################

        # # determine particle motion axes
        # particle_freqmin = 0.02 
        # particle_freqmax = 0.066   # 15 - 50s
        # motiontrace.detrend('demean')
        # motiontrace.detrend('linear') 
        # motiontrace.filter('bandpass',freqmin=particle_freqmin,freqmax=particle_freqmax)

        # motion_north           = motiontrace.select(component="N")[0].data
        # motion_east            = motiontrace.select(component="E")[0].data
        # motion_sample_interval = motiontrace[0].stats.delta


        # # makes number of samples odd for processing
        # if len(motion_north) % 2 == 0:
        #         motion_north   = motion_north[:-1]
        # else:
        #         pass
        
        # if len(motion_east) % 2 == 0: 
        #         motion_east    = motion_east[:-1]
        # else: 
        #         pass
        
        # calcbaz,xlam1,xlam2 = covar_particle_motion(motion_north[windowstart:windowend],motion_east[windowstart:windowend])


        # ### make sure that calcbaz is within same quadrant as evbaz 

        # # mar 11 2026 is this correction necessary when unraveling SI plots 0-360? 
        # # yes, this is important since the roatation correlation only looks at degrees between 0-180

        # if abs(calcbaz - evbaz.iloc[0]) > 90:
        #     calcbaz = (calcbaz+180) % 360

        # print('Calculated Polarization Direction: ',calcbaz)

        # #############################################

        # # envelope filtering for SNR (based on technique used previously
        pair    = sw.Pair(north, east, delta=sample_interval)
        # print(pair.x,pair.y)

        print('let us see why envelopes is spitting out some bs')
        # print(float(eventdata.loc[evind, 'evbaz'].iloc[0]),float(eventdata.loc[evind, 'calcbaz'].iloc[0]))
        # print(eventdata.loc[evind, 'windowstart'].iloc[0],eventdata.loc[evind, 'windowend'].iloc[0])
        # print(type(float(eventdata.loc[evind, 'evbaz'].iloc[0]),float(eventdata.loc[evind, 'calcbaz'].iloc[0]),eventdata.loc[evind, 'windowstart'].iloc[0],eventdata.loc[evind, 'windowend'].iloc[0]))
        mean_snr , max_snr, ne_env, mean_snr2, max_snr2  = envelopes(pair,float(eventdata.loc[evind, 'evbaz'].iloc[0]),float(eventdata.loc[evind, 'calcbaz'].iloc[0]),eventdata.loc[evind, 'windowstart'].iloc[0]/sample_interval - 200,(eventdata.loc[evind, 'windowend'].iloc[0]/sample_interval))
        print("envelope data:")
        print(mean_snr , mean_snr2)
        print(max_snr , max_snr2)
        print("-----")

        # # print(evbaz,calcbaz)

        print(ne_env)
        print('---------------')

        # come back to this tomorrow morning 
        
        
        # # eventdata.loc[evind,'windowstart'] = windowstart*sample_interval
        # # eventdata.loc[evind,'windowend'] = windowend*sample_interval
        # # eventdata.loc[evind,'calcbaz'] = calcbaz
        # # eventdata.loc[evind,'sta_lta'] = sta_lta_max
        # # eventdata.loc[evind,'snrbefore'] = mean_snr
        # # eventdata.loc[evind,'snrafter'] = mean_snr2
        # # eventdata.loc[evind,'ellipticity'] = xlam1/xlam2
        
        # # print('worked')
        # # bruh+=1

        # # fig_string = str('/home/gcl/BR/nicolasv03/deep_eq_plots/'+time+'.png')
        # kwargs = {'title': f"Event: {time} Phase: {eventdata['evphase'].iloc[evind]} Depth: {eventdata['evdp'].iloc[evind]} Mag: {eventdata['evmag'].iloc[evind]} Epic: {round(eventdata['evepic'].iloc[evind])} Ratio: {round(eventdata['ellipticity'].iloc[evind])}"}

        # kwargs = {'title': f"Event: {evtime.iloc[0]} Phase: {evphase.iloc[0]} Depth: {evdp.iloc[0]} Mag: {evmag.iloc[0]} Epic: {round(evepic.iloc[0],1)} Ratio: {round(xlam1/xlam2,2)}"}
        print("DEBUGGING INFO HERE")
        # print(evind)
        # print(len(eventdata))
        # print(eventdata['evbaz'].iloc[evind])
        # print(eventdata.loc[evind, 'evbaz'])

        kwargs = {'title': f"Event: {eventdata.loc[evind, 'evtime'].iloc[0]} Phase: {eventdata.loc[evind, 'evphase'].iloc[0]} Depth: {float(eventdata.loc[evind, 'evdp'].iloc[0])} Mag: {float(eventdata.loc[evind, 'evmg'].iloc[0])} Epic: {float(eventdata.loc[evind, 'evepic'])} Ratio: {float(eventdata.loc[evind, 'ellipticity'])}"}

        print(float(eventdata.loc[evind, 'evbaz']),float(eventdata.loc[evind, 'calcbaz']),int(eventdata.loc[evind, 'windowstart'].iloc[0]),int(eventdata.loc[evind, 'windowend'].iloc[0]),ne_env)
        fig = plot_trace_final(pair,float(eventdata.loc[evind, 'evbaz']),float(eventdata.loc[evind, 'calcbaz']),int(eventdata.loc[evind, 'windowstart'].iloc[0]/sample_interval)-200,(int(eventdata.loc[evind, 'windowend'].iloc[0]/sample_interval)),ne_env,**kwargs)
        fig.savefig('/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/waveform_plots_TAM/chopped'+time+'.png',dpi=800)
                    
        # if filename.endswith('.txt') and 'target_keyword' in filename:
        # src_path = os.path.join('/work/gcl3/BR/nicolasv03/senior_thesis/data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/', filtered_waveforms[i])
        # dst_path = os.path.join('/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/miniseeds/', filtered_waveforms[i])

        # src_path = os.path.join(filtered_waveforms[i])
        # dst_path = os.path.join('/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/miniseeds/'+bit)

        # # # Copy or move the file
        # shutil.copy(src_path, dst_path)     # use shutil.move() to move instead
        # print(f"Saved {filtered_waveforms[i]} to {'/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/miniseeds/'}")
                    
        #         #   'figname': fig_string
        #         #   }
        
        # # if sta_lta_max > 0:
        # # if (max_snr2 or max_snr) < 1:

        # # print("HERE IT ISSSSSSSS:")
        # # print(SplittingIntensity(pair,float(calcbaz)))
        # # print("AND THAT IS ALL")

        # if (sta_lta_max > 6) & (mean_snr > 2) & (mean_snr2 > 1.5): # & (abs(calcbaz - evbaz.iloc[0]) < 20):
        #     plot_trace_final(pair,float(evbaz.iloc[0]),float(calcbaz),windowstart,windowend,ne_env,**kwargs)

        #     print("HERE IT ISSSSSSSS:")
        #     print(SplittingIntensity(pair,float(calcbaz)))
        #     print("AND THAT IS ALL")
            # plt.savefig('/home/gcl/BR/nicolasv03/deep_eq_plots/'+time+'.png',dpi=800)
            # fig.savefig('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/FITZ/pm_plots'+time+'.png',dpi=800)
        

        # plot_trace(pair,baz,figname='../data/sample/CAN/plots/'+i+'_fr_example_plot.png',title=str(tit))
            # plt.savefig
        # sys.exit()
    except Exception as exception: 
        print('didnt work lol')
        print(f"Error: {exception}")
        print(x)
        traceback.print_exc()
        # ttime.sleep(5)
        bruh += 1
        # sys.exit()

# eventdata.to_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/eventdata.txt', sep='|', index=False)

