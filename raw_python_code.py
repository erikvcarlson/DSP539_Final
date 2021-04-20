import subprocess
import sys
import pandas as pd
import datetime as dt
import time
import numpy as np
import multiprocessing as mp
import final_proj_functions as fpf
import navpy 
import matplotlib.pyplot as plt
import csv
import convert_data_for_BRI as cdfB

startTime = time.time()


gps_data_csv = 'apr7_data.csv'

gps_data_raw = pd.read_csv(gps_data_csv)
gps_data_sorted = gps_data_raw[['Y', 'X', 'ele','time']].copy()
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 



antenna_data_csv = 'apr7_tag.csv'

antenna_1_data = fpf.antenna_puller(antenna_data_csv,1)
antenna_2_data = fpf.antenna_puller(antenna_data_csv,2)
antenna_3_data = fpf.antenna_puller(antenna_data_csv,3)
antenna_4_data = fpf.antenna_puller(antenna_data_csv,4)
antenna_5_data = fpf.antenna_puller(antenna_data_csv,5)

lat_long_of_detections_1 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_1_data)
#lat_long_of_detections_2 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_2_data)
#lat_long_of_detections_3 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_3_data)
#lat_long_of_detections_4 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_4_data)
#lat_long_of_detections_5 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_5_data)


station_location = [41.520752667,-71.599254667,40.0];
station_offset = navpy.lla2ecef(station_location[0],station_location[1],station_location[2])

ant1_xyz = fpf.xyz_converter(station_location,station_offset,lat_long_of_detections_1)
#ant2_xyz = fpf.xyz_converter(station_location,station_offset,lat_long_of_detections_2)
#ant3_xyz = fpf.xyz_converter(station_location,station_offset,lat_long_of_detections_3)
#ant4_xyz = fpf.xyz_converter(station_location,station_offset,lat_long_of_detections_4)
#ant5_xyz = fpf.xyz_converter(station_location,station_offset,lat_long_of_detections_5)    

#print(ant1_xyz)
x_1 = []
for x in range(0,len(ant1_xyz)):
    x_1.append(ant1_xyz[x][0])

y_1 = []
for y in range(0,len(ant1_xyz)):
    y_1.append(ant1_xyz[y][1])

z_1 = []
for z in range(0,len(ant1_xyz)):
    z_1.append(ant1_xyz[z][2])

rssi_1 = []
for dbm in range(0,len(ant1_xyz)):
    rssi_1.append(ant1_xyz[dbm][3])



spherical_coordinates_radius = []
spherical_coordinates_az = []

for xyz in range(0,len(ant1_xyz)):
    [radius, ele, az] = fpf.cart2sph(x_1[xyz],y_1[xyz],z_1[xyz])
    spherical_coordinates_radius.append(radius)    
    spherical_coordinates_az.append(az)



df = pd.DataFrame([[i[0],i[1],i[2],i[3]] for i in ant1_xyz],columns = ['x','y','z','rssi'])

df.to_csv('calibration_data.csv',index = False)

#plt.plot(spherical_coordinates_radius,rssi_1,'o')
#plt.savefig('foo.png', bbox_inches='tight')

fig = plt.figure()
#ax = fig.add_subplot(111, polar=True)
#ax.plot(spherical_coordinates_az,spherical_coordinates_radius,'o')
#plt.savefig('bar.png', bbox_inches='tight')


curve_set_1 = cdfB.dbm_to_curve(-100, ant1_xyz)
#print(len(ans))

ax = fig.add_subplot(projection='3d')
ax.plot(curve_set_1[:, 0], curve_set_1[:, 1], curve_set_1[:, 2])
#ax.plot(points_set_1[:, 0], points_set_1[:, 1], points_set_1[:, 2], 'o:')
plt.show()

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))

#print(lat_long_of_detections_1)
#for time in antenna_1_data['time']
