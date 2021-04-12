import subprocess
import sys
import pandas as pd
import datetime as dt

def gps_data_prepartion(gps_data_csv):
    gps_data_raw = pd.read_csv(gps_data_csv)
    gps_data_sorted = gps_data_raw[['Y', 'X', 'ele','time']].copy()
    gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time'])
    return

#gps_data_cleaned = gps_data_prepartion('apr7_data.csv')


gps_data_csv = 'apr7_data.csv'

gps_data_raw = pd.read_csv(gps_data_csv)
gps_data_sorted = gps_data_raw[['Y', 'X', 'ele','time']].copy()
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 
gps_data_sorted['time']= pd.to_datetime(gps_data_sorted['time']) 


print(gps_data_sorted)
