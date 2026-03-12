# -*- coding: utf-8 -*-

from mychron_csv import MychronCSV
from session import Session
import numpy as np
import pathlib
import matplotlib.pyplot as plt



#session_filename = pathlib.Path('../data/Clara_Kosmic_CentennialT_a_0062.csv')
session_filename = pathlib.Path('../data/Clara_Centennial_CCW_18092024_6.csv')

# Load data from file
chron = MychronCSV(session_filename)
session = Session(chron)

# Extract Laps
lap_data = []
for lap in range(1, session.num_laps()-1):
    lap_data.append(session.get_lap(lap))

# Combine laps into total dataset
time = lap_data[0]["gps_time"]
xgps = lap_data[0]["x"]
ygps = lap_data[0]["y"]
for i in range(1, len(lap_data)):
    time = np.append(time, lap_data[i]["gps_time"])
    xgps = np.append(xgps, lap_data[i]["x"])
    ygps = np.append(ygps, lap_data[i]["y"])

# Plot Acceleration
SMALL  = 10
MEDIUM = 12
LARGE  = 20
 
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(xgps, ygps, "k-")

ax.set_title("Track Positions", fontsize=LARGE)

plt.xticks(fontsize=SMALL)
ax.set_xlabel('X (m)', fontsize=MEDIUM)

plt.yticks(fontsize=SMALL)
ax.set_ylabel('Y (m)', fontsize=MEDIUM)

plt.grid()
plt.show()
