# -*- coding: utf-8 -*-

from mychron_csv import MychronCSV
from session import Session
import numpy as np
import pathlib
import matplotlib.pyplot as plt



session_filename = pathlib.Path('../data/Clara_Kosmic_CentennialT_a_0062.csv')
#session_filename = pathlib.Path('../data/Clara_Centennial_CCW_18092024_6.csv')

# Load data from file
chron = MychronCSV(session_filename)
session = Session(chron)

# Extract Laps
lap_data = []
for lap in range(1, session.num_laps()-1):
    lap_data.append(session.get_lap(lap))
    
time = lap_data[0]["gps_time"]
xgps = lap_data[0]["x"]
ygps = lap_data[0]["y"]
for i in range(1, len(lap_data)):
    time = np.append(time, lap_data[i]["gps_time"])
    xgps = np.append(xgps, lap_data[i]["x"])
    ygps = np.append(ygps, lap_data[i]["y"])
 
# Get Intersections
# -- Fast Turn
xi = -75
yi = -45.66
ft_entry_times = session.get_intersections(xi, yi, system='xyz', parameter='gps_time')

xi = -95.5
yi = --0.18
ft_exit_times = session.get_intersections(xi, yi, system='xyz', parameter='gps_time')

# -- Slow Turn
xi = 64.13
yi = 31.08
st_entry_times = session.get_intersections(xi, yi, system='xyz', parameter='gps_time')

xi = 46.915
yi = 54.30
st_exit_times = session.get_intersections(xi, yi, system='xyz', parameter='gps_time')

fast_turn_times = np.column_stack((ft_entry_times, ft_exit_times[1:])) # Remove track join
slow_turn_times = np.column_stack((st_entry_times, st_exit_times))

xpos_ft = session.get_data('gps_time', fast_turn_times[2], 'x')
ypos_ft = session.get_data('gps_time', fast_turn_times[2], 'y')
velo_ft = session.get_data('gps_time', fast_turn_times[2], 'gps_velocity')
accl_ft = session.get_data('gps_time', fast_turn_times[2], 'total_acc')  
curv_ft = session.get_data('gps_time', fast_turn_times[2], 'curvature')  
arad_ft = velo_ft**2 * curv_ft

xpos_st = session.get_data('gps_time', slow_turn_times[2], 'x')
ypos_st = session.get_data('gps_time', slow_turn_times[2], 'y')
velo_st = session.get_data('gps_time', slow_turn_times[2], 'gps_velocity')
accl_st = session.get_data('gps_time', slow_turn_times[2], 'total_acc')  
curv_st = session.get_data('gps_time', slow_turn_times[2], 'curvature')  
arad_st = velo_st**2 * curv_st

# Select out percentage
perc = 98
pvelo_ft = np.percentile(velo_ft, perc)
paccl_ft = np.percentile(accl_ft, perc)
parad_ft = np.percentile(arad_ft, perc)
pvelo_st = np.percentile(velo_st, perc)
paccl_st = np.percentile(accl_st, perc)
parad_st = np.percentile(arad_st, perc)

#
# ====================================================================
#

# Plot Acceleration
SMALL  = 10
MEDIUM = 12
LARGE  = 20

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(xgps, ygps, "k-")
ax.plot(xpos_ft, ypos_ft, "ro", label=f'Fast Turn')
ax.plot(xpos_st, ypos_st, "bo", label=f'Slow Turn')
plt.axis([-100, 120, -110, 110])
ax.legend(fontsize=SMALL)
plt.grid()
plt.show()

#
# Fast Turn Total Acceleration
#
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(velo_ft, accl_ft, "ko")
ax.hlines( y=paccl_ft, xmin=0, xmax=30, colors='r', linestyles='--', label=f'{perc}th Percentile')

plt.axis([19.5, 23, 0, 20])

plt.xticks(fontsize=SMALL)
ax.set_xlabel('Speed (m/s)', fontsize=MEDIUM)
plt.yticks(fontsize=SMALL)
ax.set_ylabel('Total Acceleration (m/s^2)', fontsize=MEDIUM)
ax.legend(fontsize=SMALL)
plt.grid()
plt.show()

print(f'Fast Turn Acceletation (GPS) = {paccl_ft}')

#
# Fast Turn Total Acceleration (Radius of Curvature)
#
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(velo_ft, arad_ft, "ko")
ax.hlines( y=parad_ft, xmin=0, xmax=30, colors='r', linestyles='--', label=f'{perc}th Percentile')

plt.axis([19.5, 23, 0, 35])

plt.xticks(fontsize=SMALL)
ax.set_xlabel('Speed (m/s)', fontsize=MEDIUM)
plt.yticks(fontsize=SMALL)
ax.set_ylabel('Total Acceleration (m/s^2)', fontsize=MEDIUM)
ax.legend(fontsize=SMALL)
plt.grid()
plt.show()

print(f'Fast Turn Acceletation (RadCurv) = {parad_ft}')

#
# Slow Turn Total Acceleration
#
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(velo_st, accl_st, "ko")
ax.hlines( y=paccl_st, xmin=0, xmax=30, colors='r', linestyles='--', label=f'{perc}th Percentile')

plt.axis([12.0, 20, 0, 20])

plt.xticks(fontsize=SMALL)
ax.set_xlabel('Speed (m/s)', fontsize=MEDIUM)
plt.yticks(fontsize=SMALL)
ax.set_ylabel('Total Acceleration (m/s^2)', fontsize=MEDIUM)
ax.legend(fontsize=SMALL)
plt.grid()
plt.show()

print(f'Slow Turn Acceletation (GPS) = {paccl_st}')

#
# Slow Turn Total Acceleration (Radius of Curvature)
#
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(velo_st, arad_st, "ko")
ax.hlines( y=parad_st, xmin=0, xmax=30, colors='r', linestyles='--', label=f'{perc}th Percentile')

plt.axis([12.0, 20, 0, 35])

plt.xticks(fontsize=SMALL)
ax.set_xlabel('Speed (m/s)', fontsize=MEDIUM)
plt.yticks(fontsize=SMALL)
ax.set_ylabel('Total Acceleration (m/s^2)', fontsize=MEDIUM)
ax.legend(fontsize=SMALL)
plt.grid()
plt.show()

print(f'Slow Turn Acceletation (RadCurv) = {parad_st}')

