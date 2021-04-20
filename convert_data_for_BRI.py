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
        install('geomdl')
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
    from geomdl import BSpline
    from geomdl.visualization import VisMPL
    from geomdl import utilities

    #take the calibration csv and convert it to the required points
    with open(calibration_csv, newline='') as f:
        reader = csv.reader(f)
        points = my_list = [[float(x) for x in rec] for rec in csv.reader(f, delimiter=',')]
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
    control_points = [control_points,control_points]
    # Create a BSpline surface
    surf = BSpline.Surface()

    # Set degrees
    surf.degree_u = 1
    surf.degree_v = 2
    print(control_points)
    # Set control points
    surf.ctrlpts2d = control_points

    # Set knot vectors
    surf.knotvector_u = [0, 0, 1]
    surf.knotvector_v = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]

    # Set evaluation delta
    surf.delta = 0.025

    # Evaluate surface points
    power_curve = surf.evaluate()

    # Import and use Matplotlib's colormaps
    from matplotlib import cm

    # Plot the control points grid and the evaluated surface
    surf.vis = VisMPL.VisSurface()
    surf.render()

    #power_curve = surf.evalpts
    #surf.sample_size = 100

    #surf.vis = VisMPL.VisSurface(ctrlpts=True, legend=False)
    #surf.render(colormap=cm.terrain)

    #convert the control points into a numpy array 
    #points_set_1 = np.array(control_points)
    #how mnay points do you want to generate the curve, in this case, I select 2000 points
    #t_points = np.arange(0, 1, 1/2000)
    #calculate the points on the curve for the constant value of rssi 
    #power_curve = Bezier.Curve(t_points, points_set_1)
    #power_curve = 1
    return power_curve

