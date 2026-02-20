#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mychron_csv import MychronCSV
from session import Session
import numpy as np
import pathlib
import matplotlib.pyplot as plt

# Constants
g = 9.80665

# Input Parameters
weight = 355.0 * 0.4535924 # kg
tire_area = 1.0 # m^2


# Constants
g = 9.80665

session_filename = pathlib.Path('../data/Clara_Kosmic_CentennialT_a_0062.csv')
#session_filename = pathlib.Path('../data/Clara_Centennial_CCW_18092024_6.csv')

# Load data from file
chron = MychronCSV(session_filename)
session = Session(chron)

# Extract Laps
lap_data = []
for lap in range(1, session.num_laps()-1):
    lap_data.append(session.get_lap(lap))

# Combine laps into total dataset
velo = lap_data[0]["gps_velocity"]
accl = lap_data[0]["total_acc"]
for i in range(1, len(lap_data)):
    velo = np.append(velo, lap_data[i]["gps_velocity"])
    accl = np.append(accl, lap_data[i]["total_acc"])

# Select out percentage
perc = 98
pvel = np.percentile(velo, perc)
pacc = np.percentile(accl, perc)

# Plot Acceleration
SMALL = 10
MEDIUM = 12
LARGE = 20
 
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(velo, accl, "ko")
ax.hlines( y=pacc, xmin=0, xmax=30, colors='r', linestyles='--', label=f'{perc}th Percentile')

plt.axis([8, 26, 0, 20])

#ax.set_title("Total Acceleration vs. Velocity", fontsize=LARGE)

plt.xticks(fontsize=SMALL)
ax.set_xlabel('Speed (m/s)', fontsize=MEDIUM)

plt.yticks(fontsize=SMALL)
ax.set_ylabel('Total Acceleration (m/s^2)', fontsize=MEDIUM)

ax.legend(fontsize=SMALL)

plt.grid()
plt.show()

# Coeficient of Friction
cof = pacc / g

print(f'Total Acceleration = {pacc}')
print(f'Coefficient of Friction = {cof}')
