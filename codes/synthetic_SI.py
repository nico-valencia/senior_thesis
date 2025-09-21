# import pandas as pd 
# import numpy as np 
# import matplotlib.pyplot as plt 
# # import scipy.stats as scistats
# # import scipy.optimize as sci 
# # import cartopy.crs as ccrs
import splitwavepy as sw
import numpy as np
# import pandas as pd 
import matplotlib.pyplot as plt
# # import pyvista as pv 
# # import gc 
# # import math 
# import copy 
# from obspy.signal.filter import envelope
# from obspy import read 
# from scipy.interpolate import griddata
# import glob
# from obspy import geodetics
# import itertools 
# from multiprocessing import Pool 
# import multiprocessing
# import os 
from support import plot_trace  #plot_eigen_measure, plot_trans_measure , lambda_filter_2 , plotdata
import random 
import itertools 
from support import SplittingIntensity, covar_particle_motion




bazs = np.linspace(0,360,21)

param_grid = {
                'fast'      : np.around(np.linspace(0,360,20),3), 
                'dt'        : np.around(np.linspace(0,1,5),3)
                
                    }
param_combination = list(itertools.product(*param_grid.values()))

data = np.zeros((len(bazs)*len(param_combination), 5))

# print(param_combination)
ind = 0 

for i in range(0,len(bazs)):
    for j in range(0,len(param_combination)):
        baz = bazs[i] 
        fast1,dt1 = param_combination[j]

        wavelet = sw.Pair(split=[(fast1,dt1),(30,0.5)], noise=0.1, pol=baz, delta=0.02,width=150,nsamps=8001,noisewidth=80)


        calcbaz,xlam1,xlam2 = covar_particle_motion(wavelet.x,wavelet.y)

        wiggleSIlow_act , wiggleSI_act , wiggleSIhigh_act = SplittingIntensity(wavelet,float(baz))
        wiggleSIlow_calc , wiggleSI_calc , wiggleSIhigh_calc = SplittingIntensity(wavelet,float(calcbaz))

        if abs(calcbaz - baz) > 90:
            calcbaz = (calcbaz+180) % 360
        # print(baz,fast1,dt1)
        print(baz,calcbaz)

        data[ind,:] = baz,calcbaz,wiggleSI_act,wiggleSI_calc,abs(wiggleSIhigh_calc-wiggleSI_calc)
        ind += 1

selected_indices = np.random.choice(len(data), size=len(data) // 4 , replace=False)
trimdata = data[selected_indices]


plt.figure(figsize=(20,10))
plt.subplot(1,3,1)

# plt.plot(data[:,0],data[:,2],'o',c='blue')
plt.plot(data[:,1],data[:,3],'o',c='blue')
plt.errorbar(data[:,1],data[:,3],yerr=data[:,4],xerr=None,fmt='none',ecolor='blue',barsabove=True)
plt.title('Synthetic S Phase Data: Varied Source-Side Anisotropy',fontsize=30)
plt.xlabel('Polarization Azimuth',fontsize=20)
plt.ylabel('Splitting Intensity',fontsize=20)

plt.subplot(1,3,3)
truncated = data[(data[:,1]>40)&(data[:,1]<60)]
plt.hist(truncated[:,3],bins=7)

# plt.plot(trimdata[:,1],trimdata[:,3],'o',c='red')
# plt.errorbar(trimdata[:,1],trimdata[:,3],yerr=trimdata[:,4],xerr=None,fmt='none',ecolor='red',barsabove=True)



# plt.show()
plt.savefig('../data/sample/synthetic.png',dpi=200)










# ang1 = -90 
# ang2 = 90 
# dt1 = 0 
# dt2 = 2

# baz = random.randint(0,360)
# # baz = 90
# print(baz)
# # pair    = sw.Pair(split=[(60,0.5),(60.0,0.8)], noise=0.005, pol=baztopolar(baz), delta=0.05)

# # wavelet = sw.Pair(split=[(30,0.5),(0.0,0.0)], noise=0.005, pol=baztopolar(baz), delta=0.05)
# wavelet = sw.Pair(split=[(30,0.6),(70.0,1.0)], noise=0.03, pol=baz, delta=0.02,width=150,nsamps=8001,noisewidth=80)
# wavelet.rotateto(baz)
# # wavelet.unsplit(40,0.9)

# kwargs = {'figname': 'synth_trace.png','cmplabels':['Radial','Transverse']}
# # plot_trace(wavelet,baz, **kwargs)
# plt.close()
# n = len(wavelet.x)
# d = 0.01
# t = np.linspace(0.0, n*d, n, endpoint=False)
# plt.plot(t,wavelet.x,label='x',color='blue')
# plt.plot(t,wavelet.y,label='y',color='red')
# plt.legend()
# plt.show()
# plt.close()
# # wavelet.unsplit(30,0.5)
# wavelet.unsplit()
# kwargs = {'figname': 'synth_trace.png'}
# plot_trace(wavelet,baz, **kwargs)



# filename = './synthetic_data_fr.txt'

# file = open(filename, "w")
# file.write("name"           +"|"+ 
#            "#Network"       +"|"+ 
#            "Station"        +"|"+ 
#            "Latitude"       +"|"+ 
#            "Longitude"      +"|"+ 
#            "fast1"          +"|"+ 
#            "lag1"           +"|"+ 
#            "fast2"          +"|"+ 
#            "lag2"           +"|"+ 
#            "dist"           +"|"+ 
#            "az"             +"|"+ 
#            "baz"            +"|"+
#            "event"          +"|"+
#            "evlat"          +"|"+
#            "evlon"          +"|"+
#            "evdepth"        +"|"+
#            "evmag"          +"|"+
#            "SNR"            +"|"+
#            "Emin"           +"|"+
#         #    "dfast"          +"|"+
#         #    "dlag"           +"|"+
#            "radsnr"         +"|"+
#            "transsnr"       +"|"+
#            "radenv"         +"|"+   
#            "transenv"       +"|"+ 
#            "\n")

# bazrange = np.arange(0,360,10)
# # for i in range(0,len(path_list),1):
# for i in range(0,len(bazrange),1):
#     print(bazrange[i])
#     # try:
#     if 1 ==1:  

#         wavelet = sw.Pair(split=[(30,0.5),(0.0,0.0)], noise=0.001, pol=bazrange[i], delta=0.05)
#         x = wavelet.splitting_intensity()
#         # print(len(wavelet.x))
#         wave = copy.deepcopy(wavelet)
#         ogwave = copy.deepcopy(wavelet)
#         R = ogwave.x
#         T = ogwave.y
#         Rmod = wave.x
#         Tmod = wave.y
#         print(x)



