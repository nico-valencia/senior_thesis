import pandas as pd 
import glob 
import os 
from obspy import geodetics
from obspy.taup import TauPyModel


'''
COMPILE ALL DATA ABOUT WAVEFORMS MEASURED AT SPECIFIED STATION, NO REPETITIONS:
1. Compile all event data for station into one dataframe and remove duplicates.
2. Construct dataframe with:
    - event time, latitude, longitude, depth, magnitude, phase
    - event distance, azimuth, back azimuth, epicentral distance
    - prepares dataframe for window selection 
    - saves dataframe to txt file for measurements
'''

### meant to work within senior thesis directory, all data within data directory
model = TauPyModel(model="iasp91")



data_tag = 'TAM'
station = 'TAM'
data_type = 'sample'
phase = ['S','SKS']
# phase = ['SKS']

### also can calculate incidence angle of the waveform right here for future purposes...



waveforms = []
for ph in phase:
    dir_waveforms   = '../data/'+data_type+'/'+station+'/'+data_tag+'_'+ph+'_waveforms/'
    waveforms.extend(glob.glob(dir_waveforms+'*'))



# dir_waveforms   = '../data/sample/'+station+'/'+station+'_'+phase+'_waveforms/'
dir_events      = '../data/'+data_type+'/'+station+'/'+data_tag+'_event_info/'
dir_station     = '../data/'+data_type+'/'+station+'/'+data_tag+'_stations_info.txt'
output_dir      = '../data/'+data_type+'/'+station+'/'+data_tag

# waveforms = glob.glob(dir_waveforms+'*')

station_df  = pd.read_csv(dir_station,sep='|',header=0)
station_lat = station_df['Latitude'].iloc[0] ; station_lon = station_df['Longitude'].iloc[0]

events_df = pd.DataFrame(columns=['evtime','evlat','evlon','evdp','evmg',
                                  'evphase','evdist','evaz','evbaz','evepic','evincident',
                                  'windowstart','windowend','calcbaz','deltaphi',
                                  'sta_lta','snrbefore','snrafter','ellipticity'])

'''
PARAMETERS:
'evtime' : event time 
'evlat' : event latitude 
'evlon' : event longitude
'evdp' : event depth 
'evmg' : event magnitude 
'evphase' : trace phase 
'evdist' : event distance 
'evaz' : event azimuth 
'evbaz' : event back azimuth 
'evepic' : event epicentral distance 
'windowstart' : start of measured window for trace 
'windowend' : end of measured window for trace 
'calcbaz' : back azimuth calculated from covariance method 
'deltaphi' : angular distance between R/T and R'/T' (particle motion coordinate system)
'sta_lta' : summed sta/lta phase maximum detection value for phase
'snrbefore' : signal noise ratio before window 
'snrafter' : signal noise ratio after window 
'ellipticity' : xlam1/xlam2 ratio from covariance function


'''

### combine all event files into one dataframe
for files in glob.glob(dir_events+'*'):
    events      = pd.read_csv(files,sep=',',header=0)
    events_df   = events_df._append(events,ignore_index=True)


# events_df.drop_duplicates()


# print(events_df)


### get event baz and epicentral distance from receiver
for i in range(len(events_df)):
    ev      = events_df.iloc[i]
    evlat   = ev['evlat']
    evlon   = ev['evlon']
    dist,az,baz = geodetics.base.gps2dist_azimuth(float(ev['evlat']),float(ev['evlon']),
                                                  float(station_lat),float(station_lon))

    epic = geodetics.base.locations2degrees(float(ev['evlat']),float(ev['evlon']),
                                            float(station_lat),float(station_lon))
    
    try:
        incident = model.get_travel_times(source_depth_in_km=events_df['evdp'].iloc[i],
                                      distance_in_degree=epic,
                                      phase_list=[events_df['evphase'].iloc[i]]
                                      )[0].incident_angle
    except Exception as exception:
        incident = None
    
    print(incident)


    events_df['evdist'].iloc[i] = dist
    events_df['evaz'].iloc[i]   = az
    events_df['evbaz'].iloc[i]  = baz
    events_df['evepic'].iloc[i] = epic
    events_df['evincident'].iloc[i] = incident

    # events_df['evphase'] = 'SKS' # only for old data

print(len(events_df))
# events_df = events_df.drop_duplicates(inplace=False)
print(len(events_df))


### remove events that do not have an mseed file 

# cleans waveforms list to just have the date of each event
waveforms = [s.split('/')[5] for s in waveforms]
waveforms = [s.split('-')[2]+'-'+s.split('-')[3]+'-'+s.split('-')[4] for s in waveforms]
waveforms = [s.replace('.mseed','') for s in waveforms]

# print(waveforms)

### gets indices of events that did not have mseed file 
dropped_event_indices = [] 
for i in range(len(events_df)):
    if events_df['evtime'].iloc[i] in waveforms:
        pass
    else:
        dropped_event_indices.append(i)

# drop specified events and save to file 
# events_df.drop(dropped_event_indices, inplace=True)

print(len(events_df))

events_df = events_df.drop_duplicates()
    # events_df['evphase'] = 'SKS' # only for old data
# events_df['evphase'].fillna('SKS', inplace=True) # do this for fitz for some reason

print(events_df)

print(len(events_df))
events_df.to_csv(output_dir+'_events_apr8.txt', sep='|', index=False)
