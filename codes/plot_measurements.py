import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy.optimize as sci 
from astropy.stats import sigma_clip



station = 'TAM'
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements.txt',sep='|')
# eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_SKS_measurements_nov4.txt',sep='|')
eventdata = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/'+station+'/'+station+'_events_measurements_mar11.txt',sep='|')
eventdata['si'] = 0

# print(eventdata)


snr_threshold = 4
snrbefore_threshold = 2.0
snrafter_threshold = 2.0
sta_lta_threshold = 6
ellipticity_threshold = 20
dt_threshold = 2.8
mag_threshold = 8.0

### truncate sta, meansnr before and after, ellip

# # actually truncate the dataset lmao 
# eventdata = eventdata[eventdata['snr'] > snr_threshold]
# eventdata = eventdata[eventdata['snrbefore'] > snrbefore_threshold]
# eventdata = eventdata[eventdata['snrafter'] > snrafter_threshold]
# eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]


# eventdata = eventdata[eventdata['ellipticity'] > ellipticity_threshold]


# eventdata = eventdata[eventdata['dt'] < dt_threshold]
print(eventdata)


# eventdata['si'].fillna(10, inplace=True)

# eventdata['calcbaz'] = eventdata['calcbaz'] % 180
# eventdata['evbaz'] = eventdata['evbaz'] % 180
eventdata = eventdata.sort_values(by='calcbaz', ascending=True)

# Define bins and categorize data
array = np.arange(0, 360, 40)
eventdata['chunk'] = np.digitize(eventdata['calcbaz'], array, right=True)


eventdata['si'] = eventdata['dt'] * np.sin(2*np.deg2rad(eventdata['calcbaz']-eventdata['phi']))


# event si errors so they can be filtered out

eventdata['SI_act_error_low'] = abs(eventdata['wiggleSI_act']-eventdata['wiggleSIlow_act'])
eventdata['SI_act_error_high'] = abs(eventdata['wiggleSI_act']-eventdata['wiggleSIhigh_act'])

eventdata['SI_calc_error_low'] = abs(eventdata['wiggleSI_calc']-eventdata['wiggleSIlow_calc'])
eventdata['SI_calc_error_high'] = abs(eventdata['wiggleSI_calc']-eventdata['wiggleSIhigh_calc'])

sierrmax = 2

eventdata = eventdata[eventdata['SI_act_error_low'] < sierrmax]
eventdata = eventdata[eventdata['SI_act_error_high'] < sierrmax]

eventdata = eventdata[eventdata['SI_calc_error_low'] < sierrmax]
eventdata = eventdata[eventdata['SI_calc_error_high'] < sierrmax]

eventdata = eventdata[eventdata['SI_calc_error_high'] != 0]


# error_SI_act = [abs(eventdata['wiggleSI_act']-eventdata['wiggleSIlow_act']),abs(eventdata['wiggleSI_act']-eventdata['wiggleSIhigh_act'])]

error_SI_act = [eventdata['SI_act_error_low'],eventdata['SI_act_error_high']]

error_SI_calc = [eventdata['SI_calc_error_low'],eventdata['SI_calc_error_high']]


# print(error_SI_calc)

# print(error_SI_calc.iloc[1])


print(eventdata)
# print(array)
# print(len(array))
total = 0 
for i in range(len(array)):
    # print(i)
    eventdatachunk = eventdata[eventdata['chunk']==i+1]
    x = sigma_clip(eventdatachunk['si'],sigma=10,maxiters=5)
    test = ~x.mask
    total += len(x)
    # print(eventdatachunk)
    # print(x)
    # print(total)
    # print(test)

    # eventdata[eventdata['chunk']==i+1]['cliptest'] = test
    eventdata.loc[eventdata['chunk'] == i + 1, 'cliptest'] = test


# print(eventdata[eventdata['cliptest']==False])

### CHANGING EVBAZ TO CALCBAZ TO SEE WHAT IS UP 

tick = 0 
# for i in range(0,len(eventdata)):

#     # if eventdata['dt'].iloc[i] >= 2.9:
#     #     eventdata['dt'].iloc[i] = 0 
#     si = eventdata['si'].iloc[i]
#     if eventdata['cliptest'].iloc[i] == True:
#         tick += 1
        
#         # SI plot with particle motion baz
#         # eventdata['si'].iloc[i] = eventdata['dt'].iloc[i] * np.sin(2*np.deg2rad(eventdata['calcbaz'].iloc[i]-eventdata['phi'].iloc[i]))
        
#         # sii = eventdata['dt'].iloc[i] * np.sin(2*np.deg2rad(eventdata['calcbaz'].iloc[i]-eventdata['phi'].iloc[i] + 2*eventdata['deltaphi'].iloc[i] ))
        
#         if eventdata['evphase'].iloc[i] == 'S':
#             plt.plot(eventdata['calcbaz'].iloc[i],si,'o',c='red',alpha=0.5)
#             # plt.plot(eventdata['evbaz'].iloc[i]%180,sii,'o',c='green',alpha=0.5)

#         if eventdata['evphase'].iloc[i] == 'SKS':
#             plt.plot(eventdata['calcbaz'].iloc[i],si,'o',c='blue',alpha=0.5)
#             # plt.plot(eventdata['evbaz'].iloc[i]%180,sii,'o',c='purple',alpha=0.5)

#     else:
#         plt.plot(eventdata['calcbaz'].iloc[i],si,'o',c='grey',alpha=0.5)



eventdatapass = eventdata[eventdata['cliptest']==True]
# # binning the results for graph
dfbaz = eventdatapass.groupby('chunk')['calcbaz'].mean().to_frame()
# dfSI = eventdatapass.groupby('chunk')['si'].mean().to_frame()                
# dfstd = eventdatapass.groupby('chunk')['si'].std(ddof=0).to_frame(name='std')

dfSI = eventdatapass.groupby('chunk')['wiggleSI_calc'].mean().to_frame()                
dfstd = eventdatapass.groupby('chunk')['wiggleSI_calc'].std(ddof=0).to_frame(name='std')
result = pd.concat([dfbaz,dfSI,dfstd],axis=1)

eventdata = eventdata.dropna()

# # SI function and best fit according to measurements or binned measurements
# def SI_function(baz,fast,lag):
#     SI = lag * 2 * np.sin(np.deg2rad(baz-fast)) * np.cos(np.deg2rad(baz-fast))
#     return SI

# SI function for double layer anisotropy according to measurements or binned measurements 
def SI_function(baz,fast1,lag1,fast2,lag2):
    SI =    (((lag1 * 2 * np.sin(np.deg2rad(baz-fast1)) * np.cos(np.deg2rad(baz-fast1)))) 
            + ((lag2 * 2 * np.sin(np.deg2rad(baz-fast2)) * np.cos(np.deg2rad(baz-fast2)))))
    return SI

# bounds = ([-360, 0], [360,3])
# p0 = [90.0 , 1.0]


# # binned measurements
# sine_fit = sci.curve_fit(SI_function,result['calcbaz'],result['wiggleSI_calc'],sigma=result['std'],bounds=bounds)
# params,cov = sine_fit

bounds = ([-360, 0, -360, 0], [360,3,360,3])

# raw measurements 
sine_fit = sci.curve_fit(SI_function,eventdata['calcbaz'],eventdata['wiggleSI_calc'],sigma=eventdata['SI_calc_error_high'],bounds=bounds)
params,cov = sine_fit

print(len(eventdata))


degs = np.arange(0,360,1)
rads = np.radians(degs)
plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.ylim(-2,2)

plt.title(station +' measurements')
plt.xlabel('Azimuth')
plt.ylabel("SI")

# plt.plot(result['calcbaz'],result['si'],'.',color='purple',markersize='10',label="Bin SI Mean")
# plt.errorbar(result['calcbaz'],result['si'],yerr=2*result['std'], xerr=None,fmt='none',ecolor='purple',barsabove=True,label="Bin Error Bar (1\u03C3)")

# for i in range(0,len(eventdata)):
#     if eventdata['evphase'].iloc[i] == 'S':

#         plt.plot(eventdata['evbaz'],eventdata['wiggleSI_act'],'.',color='red',alpha=0.5,label='SI from event back azimuth')
#         plt.errorbar(eventdata['evbaz'],eventdata['wiggleSI_act'],yerr=error_SI_act,xerr=None,fmt='none',ecolor='red',barsabove=True)
#     else:
#         plt.plot(eventdata['evbaz'],eventdata['wiggleSI_act'],'.',color='blue',alpha=0.5,label='SI from event back azimuth')
#         plt.errorbar(eventdata['evbaz'],eventdata['wiggleSI_act'],yerr=error_SI_act,xerr=None,fmt='none',ecolor='blue',barsabove=True)

# plt.plot(degs , SI_function(degs, *params))          #parameters calculated from binned back azimuths 
# eventdata['SI_calc_error_low'] = abs(eventdata['wiggleSI_calc']-eventdata['wiggleSIlow_calc'])


plt.plot(eventdata['evbaz'],eventdata['wiggleSI_act'],'.',color='blue',alpha=0.5,label='SI from event back azimuth')
plt.errorbar(eventdata['evbaz'],eventdata['wiggleSI_act'],yerr=eventdata['SI_calc_error_low'],xerr=None,fmt='none',ecolor='blue',barsabove=True)

# plt.legend()
plt.subplot(1,2,2)
plt.ylim(-2,2)

# for i in range(0,len(eventdata)):
#     if eventdata['evphase'].iloc[i] == 'S':

#         plt.plot(eventdata['calcbaz'],eventdata['wiggleSI_calc'],'.',color='red',alpha=0.5,label='SI from calculated azimuth')
#         plt.errorbar(eventdata['calcbaz'],eventdata['wiggleSI_calc'],yerr=error_SI_act,xerr=None,fmt='none',ecolor='red',barsabove=True)
#     else:
#         plt.plot(eventdata['calcbaz'],eventdata['wiggleSI_calc'],'.',color='blue',alpha=0.5,label='SI from calculated azimuth')
#         plt.errorbar(eventdata['calcbaz'],eventdata['wiggleSI_calc'],yerr=error_SI_act,xerr=None,fmt='none',ecolor='blue',barsabove=True)

plt.plot(eventdata['evbaz'],eventdata['wiggleSI_act'],'.',color='blue',alpha=0.5,label='SI from event back azimuth')
plt.errorbar(eventdata['evbaz'],eventdata['wiggleSI_act'],yerr=eventdata['SI_calc_error_low'],xerr=None,fmt='none',ecolor='blue',barsabove=True)

plt.plot(degs , SI_function(degs, *params))          #parameters calculated from binned back azimuths 

# plt.plot(eventdata['calcbaz'],eventdata['wiggleSI_calc'],'.',color='green',alpha=0.5,label='SI from calculated azimuth')
# plt.errorbar(eventdata['calcbaz'],eventdata['wiggleSI_calc'],yerr=error_SI_calc,xerr=None,fmt='none',ecolor='green',barsabove=True)

# plt.show()
plt.xlabel('Azimuth')
plt.ylabel("SI")
# plt.legend()
plt.savefig('../data/sample/'+station+'/measurements_mar11.png',dpi=200)
# print(tick)
print(params)
print(cov)

output_dir      = '../data/sample/'+station+'/'+station
# eventdata.to_csv(output_dir+'_final_plotted_events_mar10.txt', sep='|', index=False)

