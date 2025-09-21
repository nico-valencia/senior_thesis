import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy.optimize as sci 
# from astropy.stats import sigma_clip



station = 'STKA'
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_apr1.txt',sep='|')

ravel = False

single = False



snr_threshold = 6
snrbefore_threshold = 1.0
snrafter_threshold = 1.0
sta_lta_threshold = 7
ellipticity_threshold = 40
mag_threshold = 7.0
sierrmax = 0.5

s_depth = 400 

# ### truncate sta, meansnr before and after, ellip

# # # actually truncate the dataset lmao 
eventdata = eventdata[eventdata['snr'] > snr_threshold]
# eventdata = eventdata[eventdata['snrbefore'] > snrbefore_threshold]
# eventdata = eventdata[eventdata['snrafter'] > snrafter_threshold]
eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]
# eventdata = eventdata[eventdata['evmg'] < mag_threshold]

# # eventdata = eventdata[eventdata['ellipticity'] > ellipticity_threshold]

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
        eventdata['phase_weight'].iloc[i] = 0.46
    if eventdata['evphase'].iloc[i] == 'S':
        if eventdata['evdp'].iloc[i] >= s_depth:
            eventdata['phase_weight'].iloc[i] = 1.0
        else:
            eventdata['phase_weight'].iloc[i] = 0.1


eventdata['quality'] = np.sqrt((inv_error**2)+(snr**2)+(stalta**2)) * eventdata['phase_weight']

quality_threshold = np.percentile(eventdata['quality'],99)
eventdata = eventdata[eventdata['quality'] <= quality_threshold]

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

# else:
#     def SI_function(baz,fast1,lag1,fast2,lag2):
#         SI =    (((lag1 * 2 * np.sin(np.deg2rad(baz-fast1)) * np.cos(np.deg2rad(baz-fast1)))) 
#                 + ((lag2 * 2 * np.sin(np.deg2rad(baz-fast2)) * np.cos(np.deg2rad(baz-fast2)))))
#         return SI
#     bounds = ([-360, 0, -360, 0], [360, 3, 360, 3])
#     # p0 = [90.0 , 1.0 , 90,0 , 1.0]

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
    bounds = ([-360, 0, 0], [360, 3, 90])

    # raw measurements 
    # baz_and_incident = np.array(np.meshgrid(eventdata['calcbaz'], eventdata['evincident']))
    baz_and_incident = np.vstack((eventdata['calcbaz'], eventdata['evincident']))


    print("baz_and_incident shape:", np.shape(baz_and_incident))  # Should be (2, N)
    print("eventdata['wiggleSI_calc'] shape:", np.shape(eventdata['wiggleSI_calc']))  # Should be (N,)
    print("inv_weight shape:", np.shape(inv_weight))  # Should match `eventdata['wiggleSI_calc']`

    sine_fit = sci.curve_fit(SI_function,baz_and_incident,eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
    params,cov,fitdata,trash1,trash2 = sine_fit
    phi,dt,tilt = params[0],params[1],params[2]


# # raw measurements 
# sine_fit = sci.curve_fit(SI_function,eventdata['calcbaz'],eventdata['wiggleSI_calc'],sigma=inv_weight,bounds=bounds,full_output=True) #,p0=p0)
# params,cov,fitdata,trash1,trash2 = sine_fit
# phi,dt = params[0],params[1]

sks_tick = 0 
s_depth_tick = 0 
s_tick = 0 

y = (abs(eventdata['evbaz'] - eventdata['calcbaz']) + 180) % 180 #-eventdata['calcbaz']

# if abs(calcbaz - baz) > 90:
#             calcbaz = (calcbaz+180) % 360
            
colors = plt.cm.plasma((y - y.min()) / (y.max() - y.min()))




plt.figure(figsize=(10,6))
plt.ylim(-2,2)
plt.title(station +' measurements')
plt.xlabel('Azimuth')
plt.ylabel("SI")

for i in range(0,len(eventdata)):
    # if eventdata['evphase'].iloc[i] == 'SKS':
    #     sks_tick += 1
    #     plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='blue',alpha=0.5)
    #     plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='blue',barsabove=True)
    
    if eventdata['evphase'].iloc[i] == 'S': #and abs(eventdata['evincident'].iloc[i] - 27) < 5:
        if eventdata['evdp'].iloc[i] >= s_depth:
            s_depth_tick += 1
            # plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='purple',alpha=0.5)
            # plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='purple',barsabove=True)
        else:
            s_tick += 1
            # plt.plot(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],'.',color='red',alpha=0.5)
            # plt.errorbar(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],yerr=eventdata['SI_calc_error'].iloc[i],xerr=None,ecolor='red',barsabove=True)
            plot=plt.scatter(eventdata['calcbaz'].iloc[i],eventdata['wiggleSI_calc'].iloc[i],marker='.',s=100,c=y.iloc[i], alpha=1, cmap='plasma',vmin=0, vmax=50)
plt.colorbar()
if single: 
    plt.plot(degs , SI_function(degs, *params),color='green')          #parameters calculated from binned back azimuths 

else:
    incs = np.full_like(degs,27)
    vars = np.vstack((degs,incs))
    plt.plot(degs , SI_function(vars, *params),color='green')          #parameters calculated from binned back azimuths 


# plt.plot(degs , SI_function(degs, *params),color='green')          #parameters calculated from binned back azimuths 
# plt.plot(degs , SI_function(degs, *params) + fitdata['fvec'])

# plt.colorbar(plot)
plt.savefig('../data/sample/'+station+'/SI_measurements_apr7.png',dpi=200)

plt.figure(figsize=(8,8))
plt.title(station + ' spread of quality')
plt.hist(eventdata['quality'],bins=20)
# plt.hist(fitdata['fvec'],bins=20)
plt.ylabel("Number")
plt.xlabel("quality")
plt.savefig('../data/sample/'+station+'/quality_spread_apr7.png',dpi=200)


# print(params)
print(f"Total number of waveforms used: {len(eventdata)}")
print(f"SKS waveforms used: {sks_tick}")
print(f"S depth waveforms used: {s_depth_tick}")
print(f"S waveforms used: {s_tick}")
print(f"Fast Angle: {round(phi,4)} Delay Time {round(dt,4)}")
if single != True:
    print(f"Tilt Angle: {round(tilt,4)}")

print("Covariance Matrix:")
print(cov)

print(len(eventdata))
print(len(fitdata['fvec']))

import pygmt

catalog = eventdata
km = 111 * 2  #111km per degree

circle1 = 60 * km  # degree times deg of radius times 2 for diameter
circle2 = 90 * km 
circle3 = 40 * km
circles = np.zeros(len(catalog))
circles[:] = circle3
# radius_km = 200

# ### station: MEEK
# center_lat = -26.638 
# center_lon = 118.614998

### station: HYB
center_lat = 17.41867
center_lon = 78.55213

# ### station: CAN
# center_lat = -35.318715
# center_lon = 148.996325

# ### station: FITZ
# center_lat = -18.0982
# center_lon = 125.640297




# Define the figure
fig = pygmt.Figure()

# Set region (xmin, xmax, ymin, ymax)
region = [0, 300, -70, 70]

# Define stereographic projection (S) with the center at (longitude=0, latitude=90, scale=15 cm)
# projection = "S130/-13/15c"
# projection = "N130/15c"
projection = "M130/-20/15c"


# Create a stereographic plot
fig.coast(region=region, projection=projection, shorelines=True, frame="afg", land="gray", water="skyblue")

# Optional: Add points to the plot
# fig.plot(x=[0], y=[90], style="c0.5c", color="red", pen="1p,black")
# fig.plot(x=[0],y=[45])
fig.plot(x=catalog['evlon'],y=catalog['evlat'], style="c0.3c", fill="white", pen="black")
# fig.plot(x=[center_lon],y=[center_lat], style="c0.3c", fill="red", pen="black")

# fig.plot(x=[130], y=[-10], size=[90*km], style=f"SE-{radius_km1}k", pen="1.5p,purple2")
# fig.plot(x=[130], y=[-10], size=[120*km], style=f"SE-{radius_km2}k", pen="1.5p,purple2")
# fig.plot(x=[130], y=[-10], style=f"SE-{radius_km}k", pen="2p,red")
# fig.plot(x=[center_lon], y=[center_lat], size=[circle1], style="E-", pen="1.5p,purple2")
# fig.plot(x=[center_lon], y=[center_lat], size=[circle2], style="E-", pen="1.5p,purple2")

# fig.plot(x=catalog['longitude'], y=catalog['latitude'], size=circles, style="E-", pen="1.5p,purple2",transparency=90)
# fig.plot(x=catalog['longitude'], y=catalog['latitude'], size=circles, style="E-", pen="1.5p,purple2",transparency=98,fill='purple')





# Save the figure
# fig.savefig('./deep500km_2.png')
fig.savefig('../data/sample/'+station+'/eq_map_apr7.png',dpi=200)
