def dbm_to_curve(rssi, calibration_csv):
    #the dbm to curve function takes a rssi, expressed as an integer as well as a calibration csv and generates a list with the points on a bezier curve. 
    #The calibration csv must be structured as x position, y position, z postion, rssi. 
    import subprocess
    import sys
    #Implenting a function to insall more complex packages
    def install(package):
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    try:
        #ensure that all requesite packages are installed
        install('pandas')
        install('datetime')
        install('numpy')
        install('navpy')
        install('matplotlib')
    except:
        pass
    #import the required packages
    import pandas as pd
    import datetime as dt
    import time
    import numpy as np
    import multiprocessing as mp
    import navpy 
    import math as m
    import csv
    from Bezier import Bezier
    
    #take the calibration csv and convert it to the required points
    with open(calibration_csv, newline='') as f:
        reader = csv.reader(f)
        points = list(reader)
    #initalize a list to place our control points
    control_points = []
    #iterate over the points in the calibration data
    for n in range(0,len(points)):
        #check to see if the rssi is equal to the user's input
        if points[n][3] == rssi:
            #append the calibration point to the list of control points
            control_points.append(points[n])
    #if there are less than 3 control points then we set up a series of points within one meter of the antenna
    if len(control_points) < 3:
        control_points = [[1,0,0],[0,1,0],[0,0,1]]
    #append the first point to the end of the list so the curve comes back on itself  
    control_points.append(control_points[0])
    #convert the control points into a numpy array 
    points_set_1 = np.array(control_points)
    #how mnay points do you want to generate the curve, in this case, I select 2000 points
    t_points = np.arange(0, 1, 1/2000)
    #calculate the points on the curve for the constant value of rssi 
    power_curve = Bezier.Curve(t_points, points_set_1)
    
    return power_curve

