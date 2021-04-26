# Analyzing Yagi Antenna Radiation Patterns for Tracking Endangered Wildlife around the Block Island Wind Farms

This repository represents the final project of Erik Carlson for his DSP 539 class at the University of Rhode Island. This repository was last edited on 26 April 2021. 

# Installation 

This repository is meant to be run on the Quahog Jupyter Server at the University of Rhode Island. There are no operating specific requirement beyond package installation, and thus should theoretically be able to be run using any Jupyter Server running Python 3. 

# Required Packages 

All required packages that are require are installed as part of the Jupyter Script, however, they are also listed below: 

Pandas
Datetime
Numpy
Navpy 
Plotly

# Usage

This package can be used with any data that is generated as part of a calibration survey of radio antenna using the Cell Track Technologies Sensorstation. The data must be of the type CSV and be organized in the following manner: 

## Tag Data
The Tag Data CSV must have four columns and can handle N rows. It must be formatted using date in the form: YYYY-MM-DD HH:mm:SS, the antenna in form of an integer from 1:5 and then the RSSI as an integer with the units of dBm. 

## GPS Data 
The GPS Data CSV must be converted from a Garmin 130 FIT to GPX file using the https://gotoes.org/strava/Combine_GPX_TCX_FIT_Files.php website. The gpx file must then be converted to CSV using the following website: https://mygeodata.cloud/converter/gpx-to-csv

# Structure of Directories
It is imperative that no files are removed the repository directory as all directories in the main folder are required for proper running of the program. This is especially important for the Bezier.py and final_proj_functions.py which are called as part of the main Jupyter notebook. 

# Contributions
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

# License
[MIT](https://choosealicense.com/licenses/mit/)
