import pandas as pd
import os, sys, glob
import numpy as np
import enum # Enum for size units
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime as dt
import pygmt
from math import * #added for error fix

nameoffile = "/Users/BarronNguyen/ShearWaveSplitting/splittingDB.txt"
nameoffilenewdata = '/work/gcl3/BR/nicolasv03/senior_thesis/measurements/splittingDB3.txt'
nameoffile = '/work/gcl3/BR/nicolasv03/senior_thesis/measurements/splittingDB3.txt'

df = pd.read_csv(nameoffile, sep='|', encoding='cp1252')
dfnewdata = pd.read_csv(nameoffilenewdata, sep='|', encoding='cp1252') #NEW DATA

print(df)
print(df['Latitude'])


# #Hash this out with splittingDB2
# df.dropna(inplace=True)
df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
print("TEST \n")
print(df['Latitude'])
df['Longitude'] = pd.to_numeric(df['Longitude'], errors='coerce')
df['phi'] = pd.to_numeric(df['phi'], errors='coerce')
df['dt'] = pd.to_numeric(df['dt'], errors='coerce')

#New Data
dfnewdata['Latitude'] = pd.to_numeric(dfnewdata['Latitude'], errors='coerce')
dfnewdata['Latitude'] = pd.to_numeric(dfnewdata['Latitude'], errors='coerce')
dfnewdata['Longitude'] = pd.to_numeric(dfnewdata['Longitude'], errors='coerce')
dfnewdata['phi'] = pd.to_numeric(dfnewdata['phi'], errors='coerce')
dfnewdata['dt'] = pd.to_numeric(dfnewdata['dt'], errors='coerce')

# #WestPacificSplittingMapCoordinates:
# minlon, maxlon = 70.00, 170.00
# minlat, maxlat = -60.00, 70.00

# #SEAsiaSplittingMapCoordinates:
# minlon, maxlon = 60.00, 160.00
# minlat, maxlat = -20.00, 20.00

# #USASplittingMapCoordinates:
# minlon, maxlon = -132.89, -66.09
# minlat, maxlat = 22.62, 51.64

#PacificSplittingMapCoordinates
minlon, maxlon = 90.00, 190.00
minlat, maxlat = -50.00, 50.00


dftmp1 = df[(minlon < df['Longitude']) & (df['Longitude'] < maxlon)]
sks_meas_all = dftmp1[(minlat < df['Latitude']) & (df['Latitude'] < maxlat)]
sks_meas_all.reset_index(inplace=True)

#New Data
dftmp1newdata = dfnewdata[(minlon < dfnewdata['Longitude']) & (dfnewdata['Longitude'] < maxlon)]
sks_meas_all_newdata = dftmp1newdata[(minlat < dfnewdata['Latitude']) & (dfnewdata['Latitude'] < maxlat)]
sks_meas_all_newdata.reset_index(inplace=True)

# print(sks_meas_all.head())
# print(sks_meas_all)

# #Added Decimal Degrees to DMS Conversion for ErrorThatPopsUp
# def decdeg2dms(dd):
#    is_positive = dd >= 0
#    flooreddegrees = floor(dd)
#    dd = abs(dd)
#    minutes,seconds = divmod(dd*3600,60)
#    degrees,minutes = divmod(minutes,60)
#    degrees = degrees if is_positive else -degrees
#    return(str(flooreddegrees)+":"+str(minutes))

# #INPUT COORDINATES HERE:
# minlonDD, maxlonDD = -132.89, -66.09
# minlatDD, maxlatDD = 22.62, 51.64

# #OriginalInputCode
# minlon, maxlon = minlonDD, maxlonDD
# minlat, maxlat = minlatDD, maxlatDD
# dftmp1 = df[(minlon < df['Longitude']) & (df['Longitude'] < maxlon)]
# sks_meas_all = dftmp1[(minlat < df['Latitude']) & (df['Latitude'] < maxlat)]
# sks_meas_all.reset_index(inplace=True)
# print(sks_meas_all.head())

# #ModifiedInput
# modifiedminlon, modifiedmaxlon = decdeg2dms(minlonDD), decdeg2dms(maxlonDD)
# modifiedminlat, modifiedmaxlat = decdeg2dms(minlatDD), decdeg2dms(maxlatDD)
# regionDMS = [modifiedminlon, modifiedmaxlon, modifiedminlat, modifiedmaxlat]

# #variables
# boxcoordinates=[minlon, maxlon, minlat, maxlat]
# dcoord=0.1 #original 0.5
# figname='./splitting_map.png'
# frame=["WSen","x10+lLatitude(°)", "y10+lLongutide(°)"] #for xNUM, sets one tick in x/y direction every num degrees
# #frame=["a1f0.25", "WSen"], for splitting_map9.png frame=["a"], for splitting_map10.png frame=["a0.5", "y+ls"]
# #for splitting_map11 frame=["WSen","a1f1"], #alfnum where num = amount tick colored-in ex. "a1f2"
# #for splitting_map12 frame=["WSen","a1f2"]
# topo_data="@earth_relief_01m" #For default
# # topo_data = "@earth_relief_30s" #For topographic layer
# colormap=None #For default
# # colormap="etopo1" #For topographic layer
# proj="M14c" #Was M5c
# markerstyle="cc"
# markersizescale=0.06 #was 0.05(original), 0.04(For reduced), 0.1 (For enhanced color gradient)
# markerstyle_nocmap="c0.25c"
# markercolor="blue"
# pencolor="black" #must change paramters manually in pen="penwidth,pencolor" functions
# penwidth="1p"  #must change paramters manually in pen="penwidth,pencolor" functions
# dtscale=0.1*2 #*2 added
# measurement_cmap=True
# measurement_cmapnewdata=True
# markercolormap="grayC"#Original: "grayC" #FIXED: "cool" #GMT COLORMAP CODES NOT MATHPLOTLIB
# markercolormapnewdata="cool" #FIXED: "cool" #GMT COLORMAP CODES NOT MATHPLOTLIB
# colorbar=False 
# legend=True
                
# def plot_splitting_map(sks_meas_all, boxcoordinates, dcoord, figname, frame, topo_data, colormap, proj, markerstyle, markersizescale, markerstyle_nocmap, markercolor, markercolormapnewdata, pencolor, penwidth, dtscale, measurement_cmap, measurement_cmapnewdata, markercolormap, colorbar, legend):

# #Original Code
# # def plot_splitting_map(sks_meas_all, boxcoordinates,
# #                        dcoord=0.5, figname='splitting_map.png',
# #                        frame=["a1f0.25", "WSen"],
# #                        topo_data="@earth_relief_01m",
# #                        colormap=None,
# #                        proj="M5c",
# #                        markerstyle="cc",
# #                        markersizescale=0.05*2,
# #                        markerstyle_nocmap="c0.25c",
# #                        markercolor="blue",
# #                        pencolor='black',
# #                        penwidth="1p",
# #                        dtscale=0.5,
# #                        measurement_cmap=True,
# #                        markercolormap="viridis",
# #                        colorbar=False,
# #                       legend=True):


#     # Plot the shear wave splitting measurements using pygmt
#     # param sks_meas_all: pandas dataframe with `Longitude`, `Latitude`, `phi`, and `dt` columns.
#     # param boxcoordinates: list with minimum longitude, maximum longitude, minimum latitude, maximum latitude
#     # param dcoord: offset of the map from the given coordinates
#     # param figname: output figure name
#     # param topo_data: topographic data str
#     # param colormap: colormap for the topographic data. Defaults to None
#     # param proj: projection of the map
#     # param measurement_cmap: boolean. plot markers using the colormap
#     # param markerstyle: marker style for colormapped markers
#     # param markerstyle_nocmap: marker style without colormap
#     # param markercolor: marker color without colormap
#     # param pencolor: pen color for the splitting lines
#     # param dtscale: scale for delay time on the map



#     minlon, maxlon = boxcoordinates[0], boxcoordinates[1]
#     minlat, maxlat = boxcoordinates[2], boxcoordinates[3]
#     minlon, maxlon, minlat, maxlat = (
#         minlon - dcoord,
#         maxlon + dcoord,
#         minlat - dcoord,
#         maxlat + dcoord,
# )

# res = "f"

# #Original Code
# fig = pygmt.Figure()
# fig.basemap(region=[minlon, maxlon, minlat, maxlat],
#             projection=proj, frame=frame)

# #Modified DMS Code
# fig = pygmt.Figure()
# fig.basemap(region=regionDMS, projection=proj, frame=frame)


# if colormap is not None:
#     pygmt.makecpt(cmap=colormap, series="-12000/9000/1000",
#                     continuous=True)
#     fig.grdimage(
#         grid=topo_data,
#         shading=True,
#         cmap=colormap,
#     )

#     fig.coast(
#         frame=frame,
#         resolution=res,
#         shorelines=["1/0.2p,black", "2/0.05p,gray"],
#         borders=1,
#     )
# else:
#     fig.coast(land="lightgray", water="skyblue", resolution=res,
#                 shorelines=["1/0.2p,black", "2/0.05p,gray"],
#                 borders=1,)

# phivals = 90-sks_meas_all['phi'].values
# phivalsnewdata =  90-sks_meas_all_newdata['phi'].values
# dtvals = dtscale*sks_meas_all['dt'].values
# dtvalsnewdata = dtscale*sks_meas_all_newdata['dt'].values

# #Old Data
# if measurement_cmap:
#     pygmt.makecpt(cmap=markercolormap, series=[
#         round(sks_meas_all['dt'].min()-0.05, 1), round(sks_meas_all['dt'].max()+0.05, 1)])
#     fig.plot(
#         x=sks_meas_all['Longitude'],
#         y=sks_meas_all['Latitude'],
#         size=markersizescale*2**dtvals,
#         color=sks_meas_all['dt'],
#         cmap=True,
#         style=markerstyle,
#         pen="0.1p,black", #added 0.1p size, default is "black" for markers
#     )

# else:
#     fig.plot(
#         x=sks_meas_all['Longitude'].values,
#         y=sks_meas_all['Latitude'].values,
#         style=markerstyle_nocmap,
#         color=markercolor,
#         pen="0.1p, black", #added 0.1p size, default is "black" for markers
#         label="Station",
#     )

# fig.plot(
#     x=sks_meas_all['Longitude'],
#     y=sks_meas_all['Latitude'],
#     style="v0i+e",
#     pen="0.1p,black", #pen=[pencolor, penwidth], #pen="penwidth, pencolor",
#     direction=[
#         phivals,
#         dtvals/2,
#     ],  # angle (from xaxis) and magnitude
# )

# fig.plot(
#     x=sks_meas_all['Longitude'],
#     y=sks_meas_all['Latitude'],
#     style="v0i+e",
#     pen="0.1p,black", #pen=[pencolor, penwidth], #pen="penwidth, pencolor",
#     direction=[
#         180+phivals,
#         dtvals/2,
#     ],  # angle (from xaxis) and magnitude
# )

# #New Data
# if measurement_cmapnewdata:
#     pygmt.makecpt(cmap=markercolormapnewdata, series=[
#             round(sks_meas_all_newdata['dt'].min()-0.05, 1), round(sks_meas_all_newdata['dt'].max()+0.05, 1)])
#     fig.plot(
#         x=sks_meas_all_newdata['Longitude'],
#         y=sks_meas_all_newdata['Latitude'],
#         size=markersizescale*2**dtvalsnewdata,
#         color=sks_meas_all_newdata['dt'],
#         cmap=True,
#         style=markerstyle,
#         pen="0.1p,black", #added 0.1p size, default is "black" for markers
#     )

# fig.plot(
#     x=sks_meas_all_newdata['Longitude'],
#     y=sks_meas_all_newdata['Latitude'],
#     style="v0i+e",
#     pen="0.1p,black", #pen=[pencolor, penwidth], #pen="penwidth, pencolor",
#     direction=[
#         phivalsnewdata,
#         dtvalsnewdata/2,
#     ],  # angle (from xaxis) and magnitude
# )

# fig.plot(
#     x=sks_meas_all_newdata['Longitude'],
#     y=sks_meas_all_newdata['Latitude'],
#     style="v0i+e",
#     pen="0.1p,black", #pen=[pencolor, penwidth], #pen="penwidth, pencolor",
#     direction=[
#         180+phivalsnewdata,
#         dtvalsnewdata/2,
#     ],  # angle (from xaxis) and magnitude
# )

# #for the legend (box)
# if legend:
#     delaytimearrays = [1.0, 2.0, 3.0]
#     ydiff = 0.2 #original=0.3 #affects
#     # ydifflinespace = 1.2 #affects spacing vertical for splitting lines and numbers--> ydifflinespace #was 1.2
#     ydifflinespace = (maxlat-minlat)/150 #affects spacing vertical for splitting lines and nukbers--> ydifflinespace
#         #USA constant = 24.183
#         #SEPacific constant = 16
#         #WestAsia constant = 120
#         #Pacfic constant = 60
#         #Pacific M14c constant = 150

#     # xdiff = 8 #original=3 #affects spacing between splitting lines and numbers
#     xdiff = (maxlon-minlon)/15.5 #original=3 #affects spacing between splitting lines and numbers #was 8.35, smaller num = larger spacing
#     yloc = (maxlat-dcoord+2*ydiff) 
#     xloc = minlon+dcoord+0.2 #affects object spacing from box, added +1 to shift right #was 1

#     ##plots actual white box and frame + box position
#     # #adjustable variables:
#     # xscale = 0.2 #scales box laterally #original=0.2
#     # yscale = 0.2 #scales box vertically #original=0.2
#     # scalingconstant = 2 #original=2

#     #coordinates (x,y) represent center of legend
#     fig.plot(x=minlon+2, y=maxlat+0, style=f"r{len(delaytimearrays)/1.75}/{len(delaytimearrays)/3}", #where len(delaytimearrays)/NUM, NUM scales the box
#                 color="white", pen="0.5p,black") #legend border and (+2, +0) for new projection 

#     # fig.plot(x=xloc+xdiff/scalingconstant+xscale, y=yloc-len(delaytimearrays)/scalingconstant*ydiff-yscale, style=f"r{len(delaytimearrays)}/{len(delaytimearrays)/2}",
#     #             color="white", pen="0.5p,black") #legend border
    
#     # fig.plot(x=xloc+xdiff/2+xscale, y=yloc-len(delaytimearrays)/2*ydiff-yscale, style=f"r{len(delaytimearrays)}/{len(delaytimearrays)/2}",
#     #             color="white", pen="0.5p,black") #legend border


# ##for legend (objects inside box)
#     xloc += xdiff

#     for dltime in delaytimearrays:
#         yloc -= ydifflinespace #original line: yloc -= ydiff
#         penwidth = penwidth

#         ##does the splitting lines key
#         fig.plot(
#             x=xloc-xdiff,
#             y=yloc,
#             style="v0i+e",
#             pen="0.1p,black", #pen=[pencolor, penwidth], #pen="penwidth, pencolor",
#             direction=[ ##orientation+magnitude of splitting lines
#                 [0],
#                 [dtscale*dltime],
#             ],  # angle (from xaxis) and magnitude
#         )

#         #does the text
#         fig.text(x=xloc, y=yloc,
#                     text=f"{dltime} s", font="0.1c,Helvetica") #textsize changed from 0.3c to 0.1c, legend text size





# if colorbar:
#     fig.colorbar(frame=["a0.5", "y+ls"])

fig.savefig(figname, crop=True, dpi=720) #original dpi=300

# import pandas as pd
# import numpy as np
# import pygmt

# # File paths
# nameoffile = '/work/gcl3/BR/nicolasv03/senior_thesis/measurements/splittingDB3.txt'

# # Map region limits
# minlon, maxlon = 110, 160
# minlat, maxlat = -45, -10

# # Read and clean old SKS splitting data
# df = pd.read_csv(nameoffile, sep='|', encoding='cp1252')
# df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
# df['Longitude'] = pd.to_numeric(df['Longitude'], errors='coerce')
# df['phi'] = pd.to_numeric(df['phi'], errors='coerce')
# df['dt'] = pd.to_numeric(df['dt'], errors='coerce')

# # Filter data within the map region
# dftmp1 = df[(minlon < df['Longitude']) & (df['Longitude'] < maxlon)]
# sks_meas_all = dftmp1[(minlat < dftmp1['Latitude']) & (dftmp1['Latitude'] < maxlat)]
# sks_meas_all.reset_index(inplace=True)

# # Initialize the map
# fig = pygmt.Figure()

# fig.basemap(
#     region=[minlon, maxlon, minlat, maxlat],
#     projection="M15c",
#     frame=["af", '+t"SKS splitting measurements"']
# )

# # fig.coast(shorelines="1/0.1p,black", resolution="10m", water="skyblue", land="lightgray")

# # fig.gridline(
# #     frame=["af"],
# #     pen="0.5p,gray,--",
# #     interval="auto"
# # )

# fig.plot(
#     x=sks_meas_all["Longitude"],
#     y=sks_meas_all["Latitude"],
#     style="c0.15c",
#     pen="black",
#     color="black"
# )

# fig.plot(
#     x=sks_meas_all["Longitude"],
#     y=sks_meas_all["Latitude"],
#     direction=[sks_meas_all["phi"], sks_meas_all["dt"]],
#     style="v0.15c+e",
#     pen="0.5p,black"
# )

# fig.savefig('./testmap.png')

# # savefig(figname, crop=True, dpi=720) #original dpi=300