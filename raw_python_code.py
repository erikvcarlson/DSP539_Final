import subprocess
import sys
import pandas as pd
import datetime as dt
import time
import numpy as np
import multiprocessing as mp

def antenna_puller(antenna_data_csv,antenna_number):
    ant_data_raw = pd.read_csv(antenna_data_csv)
    ant_data_raw = pd.read_csv(antenna_data_csv)
    ant_data_raw['time']= pd.to_datetime(ant_data_raw['time'])
    antenna_data = ant_data_raw[ant_data_raw['antenna'] == antenna_number].copy()
    return antenna_data

def comparison_loop(antenna_data,gps_sub_data):
    for m  in range(len(antenna_data['time'])):
        for n in range(len(gps_sub_data['time'])):
            t1 =  antenna_data['time'][m].to_pydatetime().timestamp()
            t2 = gps_sub_data['time'][n].to_pydatetime().timestamp()
            time_difference = t1 - t2
            if time_difference == 0: 
                lat_long_of_detects.append(antenna_data[n])  

def latlong_of_detections_creator(gps_data, antenna_data):
    lat_long_of_detects = pd.DataFrame([])
    #print(len(antenna_data['time']))
    #print(len(gps_data['time']))

    #break the data into chunks and then run multithreaded
    print(mp.cpu_count())
    gps_split = np.array_split(gps_data,mp.cpu_count()) 
    print('split complete')
    a_pool = mp.Pool(mp.cpu_count())
    print('pool complete')
    #print(gps_split)
    a_pool.map(comparison_loop,gps_split)

    return lat_long_of_detects

gps_data_csv = 'apr7_data.csv'

gps_data_raw = pd.read_csv(gps_data_csv)
gps_data_sorted = gps_data_raw[['Y', 'X', 'ele','time']].copy()
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 



antenna_data_csv = 'apr7_tag.csv'

antenna_1_data = antenna_puller(antenna_data_csv,1)
antenna_2_data = antenna_puller(antenna_data_csv,2)
antenna_3_data = antenna_puller(antenna_data_csv,3)
antenna_4_data = antenna_puller(antenna_data_csv,4)
antenna_5_data = antenna_puller(antenna_data_csv,5)

#print(antenna_1_data)

lat_long_of_detections_1 = latlong_of_detections_creator(gps_data_sorted, antenna_1_data)
lat_long_of_detections_2 = latlong_of_detections_creator(gps_data_sorted, antenna_2_data)
lat_long_of_detections_3 = latlong_of_detections_creator(gps_data_sorted, antenna_3_data)
lat_long_of_detections_4 = latlong_of_detections_creator(gps_data_sorted, antenna_4_data)
lat_long_of_detections_5 = latlong_of_detections_creator(gps_data_sorted, antenna_5_data)
print(lat_long_of_detections_1)
#for time in antenna_1_data['time']
