library(reticulate)
use_python("/usr/bin/python3.6")

source_python("convert_data_for_BRI.py")

dbm_to_curve(-100,'calibration_data.csv' )

