import subprocess
import sys
import pandas as pd
import datetime as dt
import time
import numpy as np
import multiprocessing as mp
import final_proj_functions as fpf

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

#print(antenna_1_data)

lat_long_of_detections_1 = fpf.latlong_of_detections_creator(gps_data_sorted, antenna_1_data)
#lat_long_of_detections_2 = latlong_of_detections_creator(gps_data_sorted, antenna_2_data)
#lat_long_of_detections_3 = latlong_of_detections_creator(gps_data_sorted, antenna_3_data)
#lat_long_of_detections_4 = latlong_of_detections_creator(gps_data_sorted, antenna_4_data)
#lat_long_of_detections_5 = latlong_of_detections_creator(gps_data_sorted, antenna_5_data)
print(lat_long_of_detections_1)
#for time in antenna_1_data['time']
