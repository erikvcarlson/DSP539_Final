import convert_data_for_BRI as c
import pandas as pd
var = c.dbm_to_curve(-100, 'calibration_data.csv')

df = pd.DataFrame([[i[0],i[1],i[2],] for i in var],columns = ['x','y','z'])

df.to_csv('-100dbm_curve.csv',index = False)


print(var)
