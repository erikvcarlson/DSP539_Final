import subprocess
import sys
import pandas as pd
import datetime as dt
import time
import numpy as np
import multiprocessing as mp
import navpy 
import math as m
from Bezier import Bezier

def antenna_puller(antenna_data_csv,antenna_number):
    #################
    #Objective: Seperate data from individual antennas from an input csv
    #Inputs: antenna_data_csv must be a csv as described in tag data section of the read me file. The antenna number is an integer from 1 to 5
    #Outputs: A pandas dataframe with time and rssi detection information from a singular antenna 
    #################
    #intake the data using the pandas read_csv function 
    ant_data_raw = pd.read_csv(antenna_data_csv)
    #ensure that the time column is the appropriate datatype for use later
    ant_data_raw['time']= pd.to_datetime(ant_data_raw['time'])
    #take a copy of the data and append it to the antenna_data dataframe and retrun this dataframe to the user
    antenna_data = ant_data_raw[ant_data_raw['antenna'] == antenna_number].copy()
    return antenna_data

def comparison_loop(gps_sub_data):
    #################
    #Objective: Identify GPS detections that occured at the same time as a tag was detected 
    #Inputs: A numpy array from the split_ident function. The antenna_matrix is a global variable passed from the latlong_of_detections_creator
    #Outputs: A list which includes the latitude, longitude, elevation and RSSI of a putative detection
    #################
    #generate a blank list to append our data to 
    coordinated_data = [];
    #convert our gps time dataframe to a numpy array. This aids in the iteration over the data in the for loop
    gps_time = gps_sub_data['time'].to_numpy()
    #convert our gps_sub_data to first a numpy array and then a python list
    gps_lat_long_data = gps_sub_data[['Y','X','ele']].to_numpy().tolist()
    #iterate over every detection in the antenna matrix and the gps_time matrix to compare detection time
    for m  in range(len(ant_matrix)):
        for n in range(len(gps_time)):
            #convert detection times to epoch time
            t1 =  ant_matrix[m].astype('uint64')/(10**9)
            t2 = gps_time[n].to_pydatetime().timestamp()
            #subtract the two epoch times, if they are equal their result will be 0 
            time_difference = t1 - t2
            if time_difference == 0: 
                #if the data is the same time, add the lat, long, elevation and rssi to a holding_matrix which is then appended to the matrix passed back to the split ident matrix
                holding_matrix = gps_lat_long_data[n] + [ant_matrix_rssi[m]]
                coordinated_data.append(holding_matrix)
    return coordinated_data

def split_ident(x):
    #################
    #Objective: Aid in the mapping of the multithreaded comparision of the latlong_of_detections_creator using the comparison loop function
    #Inputs: A core number passed from the multiprocessing mapping
    #Outputs: A list which includes the latitude, longitude, elevation and RSSI of a putative detection
    #################
    gps_data_compare = gps_data_break[x]
    a = comparison_loop(gps_data_compare)
    return a 

def latlong_of_detections_creator(gps_data, antenna_data):
    #################
    #Objective: To identify the lat,long, and elevation of the tag by comparing GPS and tag detections and pairing identical transmissions
    #Inputs: The raw cleaned gps_data in the form of a numpy array, and the antenna data in the form of a pandas data frame. 
    #Outputs: A list which includes the latitude, longitude, elevation and RSSI of a putative detection 
    #################
    
    #generate a dataframe which will be passed back to the user.
    lat_long_of_detects = pd.DataFrame([])
    
    #define global arrays to pass back and forth between the various helper functions called by this function
    global ant_matrix 
    global ant_matrix_rssi
    global gps_data_break
    ant_matrix = antenna_data['time'].to_numpy()
    ant_matrix_rssi =  antenna_data['rssi'].to_numpy().tolist()
    
    #break the data into chunks and then run multithreaded
    gps_split = np.array_split(gps_data,mp.cpu_count())
    gps_data_break = gps_split
    a_pool = mp.Pool(mp.cpu_count())
    #generate a list with the number of cores on the users current device
    cores = list(range(0,mp.cpu_count()))
    lat_long_of_detects = a_pool.map(split_ident,cores)
    #combine the chunks into a master list
    lat_long_of_detects = [ent for sublist in lat_long_of_detects for ent in sublist]       
    return lat_long_of_detects


def xyz_converter(station_location,station_offset,detections):
    #################
    #Objective: Convert latitude,longitude and elevation to XYZ coordinates centered on the antenna
    #Inputs: The station location, station offset (as generated from navpy in the Jupyter notebook) and a list of detections from a specific antenna are required 
    #Outputs: A list of XYZ detections as centered on the antenna 
    #################
    
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
    #################
    #Objective: A basic mathematical function for converting X,Y,Z coordinates to spherical coordinates
    #Inputs: X,Y,Z values passed as floats
    #Outputs: A radius, elevation angle (radians) and an azimuth angle (radians)
    #################
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az

def dbm_to_curve(rssi, points):
    #################
    # THIS FUNCTION IS NOT USED FOR THIS ASSIGNMENT AND IS A WORK IN PROGRESS 
    # Objective: Generate points on a 3d Bezier Surface from a set of control points with equal RSSI 
    #################
    control_points = []
    for n in range(0,len(points)):
        if points[n][3] == rssi:
            control_points.append(points[n])
    if len(control_points) < 3:
        control_points = [[1,0,0],[0,1,0],[0,0,1]]
    control_points.append(control_points[0])
    #potentially need to sort data so close points are adjacent data points are next to each other+
    points_set_1 = np.array(control_points)
    t_points = np.arange(0, 1, 0.0005)
    power_curve = Bezier.Curve(t_points, points_set_1)
    

    return power_curve

def lat_long_converter(xyz_positions,new_station_location):
    #################
    #Objective: Convert XYZ positions from an arbitrary station to lat, long and elevation points on a new station. Used to translate calibration data onto deployed stations
    #Inputs: A list containing the lat, long and elevation of a new station. A nested list containing the XYZ detection points from the calibration data
    #Outputs: A list containing the lat,long and elevation points on the new station
    #################
    station_offset = navpy.lla2ecef(new_station_location[0],new_station_location[1],new_station_location[2])
    lat_long_location = []
    for i in range(len(xyz_positions)):
        lat_long_location.append(navpy.ecef2lla(xyz_positions[i][0:3] + station_offset))
    return lat_long_location


