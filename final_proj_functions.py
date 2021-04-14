import subprocess
import sys
import pandas as pd
import datetime as dt
import time
import numpy as np
import multiprocessing as mp
import navpy 
import math as m


def antenna_puller(antenna_data_csv,antenna_number):
    ant_data_raw = pd.read_csv(antenna_data_csv)
    ant_data_raw = pd.read_csv(antenna_data_csv)
    ant_data_raw['time']= pd.to_datetime(ant_data_raw['time'])
    antenna_data = ant_data_raw[ant_data_raw['antenna'] == antenna_number].copy()
    return antenna_data

def comparison_loop(gps_sub_data):
    coordinated_data = [];
    gps_time = gps_sub_data['time'].to_numpy()
    gps_lat_long_data = gps_sub_data[['Y','X','ele']].to_numpy().tolist()
    for m  in range(len(ant_matrix)):
        for n in range(len(gps_time)):
            t1 =  ant_matrix[m].astype('uint64')/(10**9)
            t2 = gps_time[n].to_pydatetime().timestamp()
            time_difference = t1 - t2
            if time_difference == 0: 
                holding_matrix = gps_lat_long_data[n] + [ant_matrix_rssi[m]]
                coordinated_data.append(holding_matrix)
    return coordinated_data

def split_ident(x):
    gps_data_compare = gps_data_break[x]
    a = comparison_loop(gps_data_compare)
    return a 

def latlong_of_detections_creator(gps_data, antenna_data):
    lat_long_of_detects = pd.DataFrame([])
    global ant_matrix 
    global ant_matrix_rssi
    global gps_data_break
    ant_matrix = antenna_data['time'].to_numpy()
    ant_matrix_rssi =  antenna_data['rssi'].to_numpy().tolist()
    
    #break the data into chunks and then run multithreaded
    gps_split = np.array_split(gps_data,mp.cpu_count())
    gps_data_break = gps_split
    a_pool = mp.Pool(mp.cpu_count())
    cores = list(range(0,mp.cpu_count()))
    lat_long_of_detects = a_pool.map(split_ident,cores)
    #combine the chunks into a master list
    lat_long_of_detects = [ent for sublist in lat_long_of_detects for ent in sublist]       
    return lat_long_of_detects


def xyz_converter(station_location,station_offset,detections):
    xyz_location_1 = []
    for i in range(len(detections)):
        source_location_array = detections[i][0:3] 
        source_location_xyz_offset = navpy.lla2ecef(source_location_array[0],source_location_array[1],source_location_array[2])
        source_location_xyz = source_location_xyz_offset - station_offset
        source_location_xyz[2] = source_location_array[2] - station_location[2] 
        source_location_xyz = np.append(source_location_xyz,detections[i][3])
        xyz_location_1.append(source_location_xyz.tolist())
    return xyz_location_1


def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az
