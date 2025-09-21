from support import plot_trace
import os 
import glob 
import numpy as np 
import pandas as pd 
from obspy import read
import splitwavepy as sw 
from obspy import geodetics

'''
VIEW TRACE AND PARTICLE MOTION PLOT 
- INCLUDES EVENT DATA (MAG,BAZ,EPI,DEPTH,PHASE)
- should incorporate way to save tags to df if wanted
'''

eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/STKA/STKA_events.txt',sep='|')
station = 'STKA'

output_dir      = '../data/sample/'+station+'/'+station

dir_station     = '../data/sample/'+station+'/'+station+'_stations_info.txt'

station_df  = pd.read_csv(dir_station,sep='|',header=0)
station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]

phase = ['S']

waveforms = []
for ph in phase:
    # dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'
    dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/STKA/STKA_'+ph+'_waveforms/'

    waveforms.extend(glob.glob(dir+'*'))

for i in waveforms:
        x = read(i)

        # x = read(dir+i)
        pieces  = i.split('-')
        bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
        time = bit.split('.')[0]+'.'+bit.split('.')[1]
        name = i.split('/')[-1]
        evind  = eventdata[eventdata.iloc[:, 0]==time].index
        # evlat = eventdata['evlat'].iloc[evind]
        # evlon = eventdata['evlon'].iloc[evind]
        evdp = eventdata['evdp'].iloc[evind]
        evphase = eventdata['evphase'].iloc[evind]
        evbaz = eventdata['evbaz'].iloc[evind]
        evepic = eventdata['evepic'].iloc[evind]

# evbaz , epepic
        try:

                #filtering parameters
                freqmin = 0.01        # 0.01
                freqmax = 0.125  # 0.125   (8-100s)

                # demean, detrend, and bandpass
                x.detrend('demean')
                x.detrend('linear')
                x.filter("bandpass",freqmin=freqmin,freqmax=freqmax)


                north           = x.select(component="N")[0].data
                east            = x.select(component="E")[0].data
                sample_interval = x[0].stats.delta


                # makes number of samples odd for processing
                if len(north) % 2 == 0:
                        north   = north[:-1]
                else:
                        pass

                if len(east) % 2 == 0: 
                        east    = east[:-1]
                else: 
                        pass

                # remove indices from string representation of df values
                phase = evphase.astype(str).str.cat()
                epic = evepic.astype(str).str.cat()
                depth = evdp.astype(str).str.cat()
                baz = evbaz.astype(str).str.cat()
                

                # title = f"Phase: {phase} EpicDist: {round(float(evepic),3)} Depth: {round(float(evdp),3)} BackAz: {round(float(evbaz),3)}"
                title = f"Phase: {phase} EpicDist: {epic} Depth: {depth} BackAz: {baz}"
                print(title)

                pair    = sw.Pair(north, east, delta=sample_interval)
                figname = '../data/sample/STKA/plots_detailed/'+name+'_traceandinfo_plot.png'
                print(figname)
                plot_trace(pair,baz,figname='../data/sample/STKA/plots_detailed/'+name+'_traceandinfo_plot.png',title=title)
                # plot_trace(pair,float(evbaz),figname=output_dir+'/plots_detailed/'+i+'__traceandinfo_plot.png',)


        except Exception as exception: 
                print('did not work')
    






### OLD CODE ###

# ##### EXTREMELY PRELIMINARY CODE TO VIEW WAVEFORMS BEFORE ANY PROCESSING IS DONE WHATSOEVER 
# ##### WAVEFORMS ONLY DEMEANED, DETRENDED, AND BANDPASS FILTERED

# # dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/CAN/CAN_S_waveforms/'

# # eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/CAN/_events.txt',sep='|')
# # eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/CAN/CAN_event_info/G-CAN-1987-2599-events-info-S.txt',sep=',')
# # eventdata.columns = eventdata.columns.str.replace(' ', '')

# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/CAN/CAN_events.txt',sep='|')


# # print(eventdata)

# station = 'CAN'
# phase = ['S']

# waveforms = []
# for ph in phase:
#     # dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+ph+'_waveforms/'
#     dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/CAN/CAN_'+ph+'_waveforms/'

#     waveforms.extend(glob.glob(dir+'*'))

# # print(waveforms)

# # UPDATE THIS, ISNT READING ALL WAVEFORMS BUT JUST THE ONES FROM THE MOST RECENTLY CREATED DIRECTORY

# for i in os.listdir(dir):
#     x = read(dir+i)
#     # # print(i)
#     pieces  = i.split('-')
#     bit    = pieces[2] + "-" + pieces[3] + "-" + pieces[4]
#     time = bit.split('.')[0]+'.'+bit.split('.')[1]
#     # print(time)
#     evind  = eventdata[eventdata.iloc[:, 0]==time].index
#     # print(event_index)
#     # print(float(eventdata['evdp'].iloc[evind]))
#     # print(event['evdp'])
#     # event = eventdata.iloc[evind]
#     evlat = eventdata['evlat'].iloc[evind]
#     evlon = eventdata['evlon'].iloc[evind]
#     evdp = eventdata['evdp'].iloc[evind]
#     evphase = eventdata['evphase'].iloc[evind]
#     # if not isinstance(event['evdp'], pd.Series):
#     #     print('true')

#     # if event['evdp'] >= 0:
#     # print(event['evdp'])

#     # if event['evdp'] == False:
#     print(evdp)
#     if float(evdp) > 100:
#         # print(evdp)

#         try:
            
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

#             print(len(north))
#             print(len(east))

#             lat = -35.318715
#             long = 148.996325

#             dist , az , baz = geodetics.base.gps2dist_azimuth(float(evlat),float(evlon),float(lat),float(long))


#             pair    = sw.Pair(north, east, delta=sample_interval)
#             pair.rotateto(baz)

#             tit = [evdp,evphase]

#         #     plot_trace(pair,baz,figname='../data/sample/CAN/plots/'+i+'_fr_example_plot.png',title=str(tit))
#                 # plt.savefig
#         except Exception as exception: 
#             print('didnt work lol ')

# # print(waveforms)