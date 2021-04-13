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

def comparison_loop(gps_sub_data):
    for m  in range(len(ant_matrix['time'])):
        for n in range(len(gps_sub_data['time'])):
            t1 =  ant_matrix['time'][m].to_pydatetime().timestamp()
            t2 = gps_sub_data['time'][n].to_pydatetime().timestamp()
            time_difference = t1 - t2
            if time_difference == 0: 
                print(time_difference)

def split_ident(x):
    gps_data_compare = gps_data_break[x]
    print(x)
    a = comparison_loop(gps_data_compare)
    return a 

def latlong_of_detections_creator(gps_data, antenna_data):
    lat_long_of_detects = pd.DataFrame([])
    global ant_matrix 
    global gps_data_break
    ant_matrix = antenna_data
    #break the data into chunks and then run multithreaded
    print(mp.cpu_count())
    gps_split = np.array_split(gps_data,mp.cpu_count())
    gps_data_break = gps_split
    print('split complete')
    a_pool = mp.Pool(mp.cpu_count())
    print('pool complete')
    cores = list(range(0,mp.cpu_count()))
    print(cores)
    a_pool.map(split_ident,cores)
    
    #print(gps_split)
    #pool_res = [mp.Process]a_pool.map(comparison_loop,gps_split)
    #for ii, chunk in enumerate(gps_split)
    #presult = astropy_example.parallel_map(comparison_loop, gps_data)
    
    a_pool.close()
    a_pool.join()
   # with ProcessPoolExecutor() as executor:
   #     res = executor.map(comparison_loop,gps_split) 
        
    return lat_long_of_detects
