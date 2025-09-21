import pandas as pd 
import glob 
import os 
from obspy import geodetics
import numpy as np 
import matplotlib.pyplot as plt
import math

# '''
# COMPILE ALL DATA ABOUT WAVEFORMS MEASURED AT SPECIFIED STATION, NO REPETITIONS:
# 1. Compile all event data for station into one dataframe and remove duplicates.
# 2. Construct dataframe with:
#     - event time, latitude, longitude, depth, magnitude, phase
#     - event distance, azimuth, back azimuth, epicentral distance
#     - prepares dataframe for window selection 
#     - saves dataframe to txt file for measurements
# '''

# ### meant to work within senior thesis directory, all data within data directory

# data_tag = 'swpmag7'
# station = 'swp_mag7'
# data_type = 'complete'
# phase = ['S','SKS']
# # phase = ['SKS']


# waveforms = []
# for ph in phase:
#     dir_waveforms   = '../data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/'
#     waveforms.extend(glob.glob(dir_waveforms+'*'))



# # dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+phase+'_waveforms/'
# dir_events      = '../data/'+data_type+'/'+station+'/'+data_tag+'_event_info/'
# dir_station     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
# output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag

# # waveforms = glob.glob(dir_waveforms+'*')

# station_df  = pd.read_csv(dir_station,sep='|',header=0)
# station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]

# events_df = pd.DataFrame(columns=['evtime','evlat','evlon','evdp','evmg',
#                                   'evphase','evdist','evaz','evbaz','evepic',
#                                   'windowstart','windowend','calcbaz','deltaphi',
#                                   'sta_lta','snrbefore','snrafter','ellipticity'])

# '''
# PARAMETERS:
# 'evtime' : event time 
# 'evlat' : event latitude 
# 'evlon' : event longitude
# 'evdp' : event depth 
# 'evmg' : event magnitude 
# 'evphase' : trace phase 
# 'evdist' : event distance 
# 'evaz' : event azimuth 
# 'evbaz' : event back azimuth 
# 'evepic' : event epicentral distance 
# 'windowstart' : start of measured window for trace 
# 'windowend' : end of measured window for trace 
# 'calcbaz' : back azimuth calculated from covariance method 
# 'deltaphi' : angular distance between R/T and R'/T' (particle motion coordinate system)
# 'sta_lta' : summed sta/lta phase maximum detection value for phase
# 'snrbefore' : signal noise ratio before window 
# 'snrafter' : signal noise ratio after window 
# 'ellipticity' : xlam1/xlam2 ratio from covariance function


# '''

# ### combine all event files into one dataframe
# for files in glob.glob(dir_events+'*'):
#     events      = pd.read_csv(files,sep=',',header=0)
#     events_df   = events_df._append(events,ignore_index=True)


# # events_df.drop_duplicates()


# # print(events_df)


# ### get event baz and epicentral distance from receiver
# for i in range(len(events_df)):
#     ev      = events_df.iloc[i]
#     evlat   = ev['evlat']
#     evlon   = ev['evlon']
#     dist,az,baz = geodetics.base.gps2dist_azimuth(float(ev['evlat']),float(ev['evlon']),
#                                                   float(station_lat),float(station_lon))

#     epic = geodetics.base.locations2degrees(float(ev['evlat']),float(ev['evlon']),
#                                             float(station_lat),float(station_lon))

#     events_df['evdist'].iloc[i] = dist
#     events_df['evaz'].iloc[i]   = az
#     events_df['evbaz'].iloc[i]  = baz
#     events_df['evepic'].iloc[i] = epic

#     # events_df['evphase'] = 'SKS' # only for old data

# print(len(events_df))
# # events_df = events_df.drop_duplicates(inplace=False)
# print(len(events_df))


# ### remove events that do not have an mseed file 

# # cleans waveforms list to just have the date of each event
# waveforms = [s.split('/')[5] for s in waveforms]
# waveforms = [s.split('-')[2]+'-'+s.split('-')[3]+'-'+s.split('-')[4] for s in waveforms]
# waveforms = [s.replace('.mseed','') for s in waveforms]

# # print(waveforms)

# ### gets indices of events that did not have mseed file 
# dropped_event_indices = [] 
# for i in range(len(events_df)):
#     if events_df['evtime'].iloc[i] in waveforms:
#         pass
#     else:
#         dropped_event_indices.append(i)

# # drop specified events and save to file 
# # events_df.drop(dropped_event_indices, inplace=True)

# print(len(events_df))

# events_df = events_df.drop_duplicates()
#     # events_df['evphase'] = 'SKS' # only for old data
# # events_df['evphase'].fillna('S', inplace=True) # do this for fitz for some reason

# print(events_df)

# print(len(events_df))
# events_df.to_csv(output_dir+'_events_dec4.txt', sep='|', index=False)


### rose diagram of back azimuths for individual stations
print('test')
station = 'TOO'
data_type = 'sample'
# file_dir      = '../data/'+data_type+'/'+station+'/'+station
file_dir = '/work/gcl3/BR/nicolasv03/senior_thesis/data/complete/australia/AU_events_apr16.txt'

data = pd.read_csv(file_dir,sep='|')
# data = pd.read_csv(file_dir+'_events_dec6.txt',sep='|')


data = data[data['evdp']>=0]

data = data[data['evmg']>=6.5]


s_phase = data[data['evphase']=='S']
sks_phase = data[data['evphase']=='SKS']

# print(data[''])

s_bazs = np.deg2rad(s_phase['evbaz'])
sks_bazs = np.deg2rad(sks_phase['evbaz'])

# Create the rose diagram
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))

# Histogram of directional data
counts_s, bin_edges_s = np.histogram(s_bazs, bins=24)
width = np.diff(bin_edges_s)  # Width of each bin

# Histogram of directional data
counts_sks, bin_edges_sks = np.histogram(sks_bazs, bins=24)
width = np.diff(bin_edges_s)  # Width of each bin

base = 500
log_counts_s = np.log(counts_s + 1) / np.log(base)
log_counts_sks = np.log(counts_sks + 1) / np.log(base)

# log_counts_s = math.log(5,counts_s + 1)
# log_counts_sks = math.log(5,counts_sks + 1)

# 
normalized_counts_s = counts_s / counts_s.sum() * 100


counts_s = np.clip(counts_s, 0, 5000)
counts_sks = np.clip(counts_sks, 0, 7000)

# Plot on polar coordinates
ax.bar(bin_edges_s[:-1], counts_s, width=width, align='edge', edgecolor='black', alpha=0.5,label='S Phases',color='red')
ax.bar(bin_edges_s[:-1], counts_sks, width=width, align='edge', edgecolor='black', alpha=0.5,label='SKS Phases',color='blue')


# Add labels and title
ax.set_title('Event Distribution for '+station, va='bottom')
ax.legend(loc="best", bbox_to_anchor=(1.2, 1.1))
ax.set_theta_direction(-1)  # Set direction to clockwise
ax.set_theta_zero_location('N')  # Set 0 degrees to north

# Show the diagram
# plt.show()
# ax.set_xticks(np.deg2rad(np.arange(0, 360, 45)))
# ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
plt.savefig('./'+station+'_rose_diagram.png')
