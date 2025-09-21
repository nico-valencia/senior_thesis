# import pandas as pd 
# import glob 
# import os 
# from obspy import geodetics
# from obspy.taup import TauPyModel
# import sys


# '''
# COMPILE ALL DATA ABOUT WAVEFORMS MEASURED AT NETWORK, NO REPETITIONS:
# 1. Compile all event data for network into one dataframe and remove duplicates.
# 1a. Make sure each event has corresponding station name and lat/lon
# 2. Construct dataframe with:
#     - event time, latitude, longitude, depth, magnitude, phase
#     - event distance, azimuth, back azimuth, epicentral distance
#     - prepares dataframe for window selection 
#     - saves dataframe to txt file for measurements
# '''

# ### meant to work within senior thesis directory, all data within data directory
# model = TauPyModel(model="iasp91")

# data_tag = 'swp'
# station = 'swp'
# data_type = 'complete'
# phase = ['S','SKS']
# # phase = ['SKS']

# ### also can calculate incidence angle of the waveform right here for future purposes...
# waveforms = []
# for ph in phase:
#     dir_waveforms   = '../data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/'
#     waveforms.extend(glob.glob(dir_waveforms+'*'))



# dir_events      = '../data/'+data_type+'/'+station+'/'+data_tag+'_event_info/'
# dir_stations     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
# output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag

# # waveforms = glob.glob(dir_waveforms+'*')

# stations_df  = pd.read_csv(dir_stations,sep='|',header=0)
# # station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]

# # print(stations_df)
# # print(waveforms)


# events_df = pd.DataFrame(columns=['evnetwork','evstation','evstationlat','evstationlon',
#                                   'evtime','evlat','evlon','evdp','evmg',
#                                   'evphase','evdist','evaz','evbaz','evepic','evincident',
#                                   'windowstart','windowend','calcbaz','deltaphi',
#                                   'sta_lta','snrbefore','snrafter','ellipticity'])

# '''
# PARAMETERS:
# 'evnetwork' : station network
# 'evstation' : station name 
# 'evstationlat' : station latitude
# 'evstationlon' : station longitude
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
#     print(files)
#     network = files.split('/')[-1].split('-')[0]
#     station = files.split('/')[-1].split('-')[1]
#     print(network)
#     print(station)
#     events      = pd.read_csv(files,sep=',',header=0)
#     events['evnetwork'] = network
#     events['evstation'] = station
#     events_df   = events_df._append(events,ignore_index=True)
#     # print('done with one of the files!')

# print('done now')
# print(events_df.head())
# print(stations_df.head())
# print('length of events df:',len(events_df))
# print('length of stations df:',len(stations_df))

# # sys.exit()

# stations_df['Station'] = stations_df['Station'].str.strip().str.upper()
# events_df['evstation'] = events_df['evstation'].str.strip().str.upper()


# ### get event baz and epicentral distance from receiver
# for i in range(len(events_df)):
#     ev      = events_df.iloc[i]
#     evlat   = ev['evlat']
#     evlon   = ev['evlon']
#     print(i)

#     station_row = stations_df[stations_df['Station'] == ev['evstation']]
#     # print(station_row)
#     if not station_row.empty:
#         station_lat = station_row['Latitude'].values[0]
#         station_lon = station_row['Longitude'].values[0]
#         # print(station_lat,station_lon)
#     else:
#         print(f"[Warning] No station info found for {ev['evstation']}")

#     # evstationlat = stations_df[stations_df['Station']==]
#     # correct_stat = stations_df[stations_df['Station']==ev['evstation']]
#     # print(ev['evstation'])
#     # print(correct_stat)
#     # station_lat = stations_df[stations_df['Station']==ev['evstation']]['Latitude'][0]
#     # station_lon = stations_df[stations_df['Station']==ev['evstation']]['Longitude'][0]
#     # print(station_lat,station_lon)
#     # station_lon = stations_df['Longitude']
#     # print(evlat,evlon)
#     # print(type(station_lat))
# # SOMEDTHING IS UP HERE FIX IT TOMORROW 
#     # sys.exit()


#     dist,az,baz = geodetics.base.gps2dist_azimuth(float(ev['evlat']),float(ev['evlon']),
#                                                   float(station_lat),float(station_lon))

#     epic = geodetics.base.locations2degrees(float(ev['evlat']),float(ev['evlon']),
#                                             float(station_lat),float(station_lon))
    
#     try:
#         incident = model.get_travel_times(source_depth_in_km=events_df['evdp'].iloc[i],
#                                       distance_in_degree=epic,
#                                       phase_list=[events_df['evphase'].iloc[i]]
#                                       )[0].incident_angle
#     except Exception as exception:
#         incident = None
    
#     # print(incident)


#     # events_df['evdist'].iloc[i] = dist
#     # events_df['evaz'].iloc[i]   = az
#     # events_df['evbaz'].iloc[i]  = baz
#     # events_df['evepic'].iloc[i] = epic
#     # events_df['evincident'].iloc[i] = incident

#     events_df.loc[i, ['evstationlat', 'evstationlon', 'evdist', 'evaz', 'evbaz', 'evepic', 'evincident']] = [station_lat, station_lon, dist, az, baz, epic, incident]

#     # events_df['evphase'] = 'SKS' # only for old data

# # print(len(events_df))
# # events_df = events_df.drop_duplicates(inplace=False)
# print(events_df.head())
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
# # events_df['evphase'].fillna('SKS', inplace=True) # do this for fitz for some reason

# print(events_df.head())

# print(len(events_df))

# # events_df = events_df.sort_values(by='Station', ascending=True)

# events_df.to_csv(output_dir+'_events_may6.txt', sep='|', index=False)


# import pandas as pd
# import glob
# import os
# from obspy import geodetics
# from obspy.taup import TauPyModel

# # Parameters
# model = TauPyModel(model="iasp91")
# data_tag = 'swp'
# station = 'swp'
# data_type = 'complete'
# phase = ['S', 'SKS']

# # Directories
# waveforms = []
# for ph in phase:
#     dir_waveforms = f'../data/{data_type}/{station}/{data_tag}_{ph}_waveforms/'
#     waveforms.extend(glob.glob(dir_waveforms + '*'))

# dir_events = f'../data/{data_type}/{station}/{data_tag}_event_info/'
# station_file = f'../data/{data_type}/{station}/{data_tag}_stations_info.txt'
# output_path = f'../data/{data_type}/{station}/{data_tag}_events_may6.txt'

# # Load station data
# stations_df = pd.read_csv(station_file, sep='|', header=0)
# stations_df['evstation'] = stations_df['Station'].str.strip().str.upper()
# stations_df = stations_df.rename(columns={'Latitude': 'evstationlat', 'Longitude': 'evstationlon'})

# # Compile event files
# events_df = pd.DataFrame()
# for file in glob.glob(dir_events + '*'):
#     network = os.path.basename(file).split('-')[0]
#     sta = os.path.basename(file).split('-')[1]
#     ev = pd.read_csv(file, sep=',', header=0)
#     ev['evnetwork'] = network
#     ev['evstation'] = sta.strip().upper()
#     events_df = pd.concat([events_df, ev], ignore_index=True)

# print(f"Loaded {len(events_df)} event rows")

# # Merge station metadata
# events_df['evstation'] = events_df['evstation'].str.strip().str.upper()
# events_df = events_df.merge(stations_df[['evstation', 'evstationlat', 'evstationlon']], on='evstation', how='left')

# # Compute geodetic info
# def compute_geodetics(row):
#     try:
#         dist, az, baz = geodetics.base.gps2dist_azimuth(row['evlat'], row['evlon'],
#                                                         row['evstationlat'], row['evstationlon'])
#         epic = geodetics.base.locations2degrees(row['evlat'], row['evlon'],
#                                                 row['evstationlat'], row['evstationlon'])
#     except Exception:
#         return pd.Series([None, None, None, None])
#     return pd.Series([dist, az, baz, epic])

# events_df[['evdist', 'evaz', 'evbaz', 'evepic']] = events_df.apply(compute_geodetics, axis=1)

# # Compute incidence angle
# def get_incident_angle(row):
#     try:
#         arrivals = model.get_travel_times(source_depth_in_km=row['evdp'],
#                                           distance_in_degree=row['evepic'],
#                                           phase_list=[row['evphase']])
#         return arrivals[0].incident_angle if arrivals else None
#     except:
#         return None

# events_df['evincident'] = events_df.apply(get_incident_angle, axis=1)

# # Prepare waveform list for filtering
# waveform_set = set()
# for wf in waveforms:
#     base = os.path.basename(wf).replace('.mseed', '')
#     parts = base.split('-')
#     if len(parts) >= 5:
#         event_time = f"{parts[2]}-{parts[3]}-{parts[4]}"
#         waveform_set.add(event_time)

# # Filter events by available waveforms
# events_df = events_df[events_df['evtime'].isin(waveform_set)]

# # Drop duplicates
# events_df = events_df.drop_duplicates()

# # Save result
# events_df.to_csv(output_path, sep='|', index=False)

# print(f"Final event count: {len(events_df)}")
# print(f"Saved to: {output_path}")


######## parallel version 

import pandas as pd
import glob
import os
from obspy import geodetics
from obspy.taup import TauPyModel
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ----------------- PARAMETERS -----------------
# data_tag = 'swp'
# station = 'swp'
# data_type = 'complete'
data_tag = 'AU'
station = 'australia'
data_type = 'complete'
phase = ['S', 'SKS']
max_workers = 64  # Adjust this for parallelism
model_name = "iasp91"

# ----------------- FILE PATHS -----------------
waveforms = []
for ph in phase:
    dir_waveforms = f'../data/{data_type}/{station}/{data_tag}_{ph}_waveforms/'
    waveforms.extend(glob.glob(dir_waveforms + '*'))

dir_events = f'../data/{data_type}/{station}/{data_tag}_event_info/'
station_file = f'../data/{data_type}/{station}/{data_tag}_stations_info.txt'
output_path = f'../data/{data_type}/{station}/{data_tag}_events_sep16.txt'

# ----------------- LOAD STATION DATA -----------------
stations_df = pd.read_csv(station_file, sep='|', header=0)
stations_df['evstation'] = stations_df['Station'].str.strip().str.upper()
stations_df = stations_df.rename(columns={'Latitude': 'evstationlat', 'Longitude': 'evstationlon'})

# ----------------- COMPILE EVENTS -----------------
events_df = pd.DataFrame()
for file in glob.glob(dir_events + '*'):
    network = os.path.basename(file).split('-')[0]
    sta = os.path.basename(file).split('-')[1]
    ev = pd.read_csv(file, sep=',', header=0)
    ev['evnetwork'] = network
    ev['evstation'] = sta.strip().upper()
    events_df = pd.concat([events_df, ev], ignore_index=True)

print(f"Loaded {len(events_df)} event rows")

# ----------------- MERGE STATION COORDINATES -----------------
events_df['evstation'] = events_df['evstation'].str.strip().str.upper()
events_df = events_df.merge(stations_df[['evstation', 'evstationlat', 'evstationlon']], on='evstation', how='left')

print('ok now we are here')
# ----------------- CALCULATE GEODETIC METRICS -----------------
def compute_geodetics(row):
    try:
        dist, az, baz = geodetics.base.gps2dist_azimuth(row['evlat'], row['evlon'],
                                                        row['evstationlat'], row['evstationlon'])
        epic = geodetics.base.locations2degrees(row['evlat'], row['evlon'],
                                                row['evstationlat'], row['evstationlon'])
        return pd.Series([dist, az, baz, epic])
    except Exception:
        return pd.Series([None, None, None, None])

events_df[['evdist', 'evaz', 'evbaz', 'evepic']] = events_df.apply(compute_geodetics, axis=1)

print('just finished computing geodetic info')

# ----------------- PARALLELIZE INCIDENT ANGLE -----------------
def compute_incident_angle(input_tuple):
    evdp, epic, phase = input_tuple
    try:
        model_local = TauPyModel(model=model_name)
        arrivals = model_local.get_travel_times(source_depth_in_km=evdp,
                                                distance_in_degree=epic,
                                                phase_list=[phase])
        return arrivals[0].incident_angle if arrivals else None
    except Exception:
        return None

input_data = list(zip(events_df['evdp'], events_df['evepic'], events_df['evphase']))
results = []

print("Computing incident angles with TauPyModel (parallel)...")
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = {executor.submit(compute_incident_angle, item): i for i, item in enumerate(input_data)}
    for future in tqdm(as_completed(futures), total=len(futures)):
        results.append(future.result())

# Apply incident angle results (preserving order)
events_df['evincident'] = results

# ----------------- MATCH EVENT TIMES TO WAVEFORMS -----------------
waveform_set = set()
for wf in waveforms:
    base = os.path.basename(wf).replace('.mseed', '')
    parts = base.split('-')
    if len(parts) >= 5:
        event_time = f"{parts[2]}-{parts[3]}-{parts[4]}"
        waveform_set.add(event_time)

events_df = events_df[events_df['evtime'].isin(waveform_set)]

# Ensure all expected columns exist
required_cols = [
    'evnetwork','evstation','evstationlat','evstationlon',
    'evtime','evlat','evlon','evdp','evmg',
    'evphase','evdist','evaz','evbaz','evepic','evincident',
    'windowstart','windowend','calcbaz','deltaphi',
    'sta_lta','snrbefore','snrafter','ellipticity'
]

for col in required_cols:
    if col not in events_df.columns:
        events_df[col] = None


# ----------------- FINAL CLEANUP & SAVE -----------------
events_df = events_df.drop_duplicates()
events_df.to_csv(output_path, sep='|', index=False)
print(f"Saved {len(events_df)} events to {output_path}")