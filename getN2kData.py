#!/usr/bin/python3

import argparse
import matplotlib.pyplot as plt
import sys
#import matplotlib as mpl
import json
import numpy as np
import statistics as stat
import scipy.optimize as opt;
import datetime
import math
import numpy as np

from modules.Boat import Boat
from modules.Utils import * 

import csv

def func(x, a, b, c,d):
     return a * np.exp(-b * x) + c

boat = Boat(
  Name = "Outremer51",
  Lwl = 15.41,      #Length of waterline in meters
  Bwl = 1.24,       #Beam of each hull
  Tc = 0.90,        #Maximum draft of haull body
  i = 5.44,         #half entry angle of the waterline at bow (neglecting the local rounded shape at the stem)
  Sw = 29.05,       #wetted area of one hull
  Dc = 6.8779,      #Displacment on one hull in m3
  Cp = 0.573,       #hull prismatic coefficient
  At = 0.00,        #immersed part of the transom area at zero speed
  LCB = 46.4,       # Longitudinal center of buoyancy, in % of Lwl from aft perpendicular
  S = 5.72,         #Space between hulls, transversal space between the 2 hulls axis
  S_aero = 11.35,   #Frontal Secion 
  Cx_aero = 0.40,   #The Cx to take into account for the wind drag comptutation (to estimate due to the overall shape of the superstructure)
  Prop_eff = 0.54,  #to estimate with the one for the propeller and the one of the mechanical transmission. Example : 0.54 = 0.6 (propeller) x 0.9 (mechaniccal transmission)
  Cms = 0.70,     # midship section coefficient (section / Tc Bw)
  Cwp = 0.70      #waterplane area coefficient (Sf / Lw Bw)
)

mps2knots= 0.000539957 * 60 * 60  #1m=0.000539957 nautical mile

parser = argparse.ArgumentParser(description='Process N2K Data.')
parser.add_argument('--from-time', type=int, default=0, required=False,
                    help='Process from this time')
parser.add_argument('--to-time', type=int,  default=-1, required=False,
                    help='Process to this time')

parser.add_argument('--y1-data', type=str, required=False, 
                    help="The data to plot 'STW,SOG'")

parser.add_argument('--y2-data', type=str, required=False, 
                    help="The data to plot on the other y axis 'STW,SOG'")

parser.add_argument('--fit-data', type=str, required=False, 
                    help="Fit a curve to the following data, example 'STW,SOG'")

parser.add_argument('--output-data', type=str, required=False, 
                    help="Output data to a csv file")

parser.add_argument('--title', type=str, default="", required=False, 
                    help="Title for the graph")

args = parser.parse_args()

data1_keys = []
data2_keys = []
fit_data_keys = []

if (args.y1_data):
  data1_keys = [item for item in args.y1_data.split(',')]
if (args.y2_data):
  data2_keys = [item for item in args.y2_data.split(',')]
if (args.fit_data):
  fit_data_keys = [item for item in args.fit_data.split(',')]

all_data_keys = data1_keys + data2_keys

line = sys.stdin.readline()
date = None
boat_data = {}
boat_data['t'] = 0

from_datetime = 0
to_datetime = 0

csv_writer = None
if args.output_data:
  data_file = open(args.output_data, 'w')
  csv_writer = csv.writer(data_file)
  csv_keys = ['t', 'date', 'time'] + all_data_keys
  csv_writer.writerow(csv_keys)

plot_data = {}
while line:
  data = json.loads(line) 

  pgn = -1
  if 'pgn' in data:
    pgn = data['pgn']
    src = data['src']
    if 'fields' in data:
      fields = data['fields']

  try:
    #Get the time
    boat_data['time'] = data['timestamp']
    boat_data['t'] = boat_data['t'] + 1
   
    if (boat_data['t'] == args.from_time):
      from_datetime = boat_data['date'] + ' T' + boat_data['time']
    if (boat_data['t'] == args.to_time):
      to_datetime = boat_data['date'] + ' T' + boat_data['time']

    if (args.to_time > 0 and boat_data['t'] > args.to_time):
      break

    if (pgn == 129029):
      if 'Date' in fields:
        boat_data['date'] = fields['Date']
        if (from_datetime == 0):
          from_datetime = boat_data['date'] + ' T' + boat_data['time']
    
    if (pgn == 127257 and src == 12):
      boat_data['pitch'] = fields['Pitch']
      boat_data['roll'] = fields['Roll']

    if (pgn == 129025 and src == 12):
      boat_data['lat'] = fields['Latitude']
      boat_data['long'] = fields['Longitude']

    if (pgn == 127251 and src == 12):
      boat_data['turn_rate'] = fields['Rate']

    if (pgn == 127250 and src == 12):
      boat_data['HDG'] = fields['Heading']

    if (pgn == 127489 and src == 57):
      boat_data['port_temp'] = fields['Temperature']
      boat_data['port_alt'] = fields['Alternator Potential']
      boat_data['port_hours'] = fields['Total Engine hours']

    if (pgn == 127488 and src == 57):
      boat_data['port_rpm'] = fields['Speed']

    if (pgn == 127489 and src == 56):
      boat_data['stb_temp'] = fields['Temperature']
      boat_data['stb_alt'] = fields['Alternator Potential']
      boat_data['stb_hours'] = fields['Total Engine hours']

    if (pgn == 127488 and src == 56):
      boat_data['stb_rpm'] = fields['Speed']

    if (pgn == 127245 and src == 17):
      boat_data['rudder'] = fields['Position']

    if (pgn == 129026 and src == 13):
      boat_data['SOG'] = fields['SOG']
      boat_data['COG'] = fields['COG']

    if (pgn == 128259 and src == 15):
      boat_data['STW'] = fields['Speed Water Referenced']
      boat_data['STW_corr'] = 0.8498807 * fields['Speed Water Referenced'] + 0.11

    
    if (pgn == 130306 and src == 15):
      if (fields['Reference'] == 'Apparent'):
        boat_data['AWA'] = fields['Wind Angle']
        boat_data['AWS'] = fields['Wind Speed']
      if (fields['Reference'] == 'True (boat referenced)'):
        boat_data['TWA'] = fields['Wind Angle']
        boat_data['TWS'] = fields['Wind Speed']
      if (fields['Reference'] == 'True (ground referenced to North)'):
        boat_data['TWD'] = fields['Wind Angle']
        boat_data['TWDS'] = fields['Wind Speed']
 
    if (boat_data): 
      if (boat_data['t'] > args.from_time):
        for k in all_data_keys:
          if k not in plot_data:
            plot_data[k] = {'x': [], 'y': []}
          if k in boat_data:
            plot_data[k]['x'].append(boat_data['t'])
            plot_data[k]['y'].append(boat_data[k])
      #print(json.dumps(boat_data))
   
    if (csv_writer):
      keys_values = []
      for k in csv_keys:
        if k in boat_data:
          keys_values.append(str(boat_data[k]))
        else:
          keys_values.append('')
      csv_writer.writerow(keys_values)
  except Exception as e:
      print("Error in line:", line, ":", e)
  line = sys.stdin.readline()

if (to_datetime == 0):
  to_datetime = boat_data['date'] + ' T' + boat_data['time']


from_datetime = datetime.datetime.strptime(from_datetime, "%Y.%m.%d T%H:%M:%S.%f")
to_datetime = datetime.datetime.strptime(to_datetime, "%Y.%m.%d T%H:%M:%S.%f")
total_seconds = (to_datetime-from_datetime).total_seconds()
total_time_tics = args.to_time - args.from_time

print("Show Plots:", from_datetime, to_datetime, total_seconds, total_time_tics)
dt = total_time_tics/total_seconds

ax1 = plt.subplot()
for k in data1_keys:
  print(k, "mean=", round(stat.mean(plot_data[k]['y']),2), "stdev=", str(round(stat.stdev(plot_data[k]['y']),2)))
  info = str(round(stat.mean(plot_data[k]['y']),2)) + " +/- " + str(round(stat.stdev(plot_data[k]['y']),2))
  ax1.plot(plot_data[k]['x'], plot_data[k]['y'], "-", label=k + " " + info)
  ax1.legend(loc='lower left')

plt.xlabel('Frames')
if (data2_keys):
    ax2 = ax1.twinx()
    for k in data2_keys:
      print(k, "mean=", round(stat.mean(plot_data[k]['y']),2), "stdev=", str(round(stat.stdev(plot_data[k]['y']),2)))
      info = str(round(stat.mean(plot_data[k]['y']),2)) + " +/- " + str(round(stat.stdev(plot_data[k]['y']),2))
      ax2.plot(plot_data[k]['x'], plot_data[k]['y'], '--', label=k + " " + info)
      ax2.legend(loc='lower right')

plt.ylabel('V(m/s)')


print(",".join(str(round(stat.mean(plot_data[k]['y']),2)) + "," + str(round(stat.stdev(plot_data[k]['y']),2)) for k in all_data_keys))


#Curve fitting
for k in fit_data_keys:
  print("Fit " + k)
  ax = (np.array(plot_data[k]['x'])- plot_data[k]['x'][0])/dt #Subtract the start so we have t=0
  ay = np.array(plot_data[k]['y'])
  ay1 = np.array(plot_data[k]['y'])
  optimizedParameters, pcov = opt.curve_fit(func, ax , ay)
  #     p0=[5.4, 2.7e-4, 9.37e-1, 1]) #, bounds=([0,0,0,0], [10, 20,100,20])) 
  print("Fuction: ", optimizedParameters)
  ay1 = func(ax, *optimizedParameters)
  ax3 = ax1.twiny()
  ax3.plot(ax, ay1 , label="fit - " + k + "eq: {}*exp(-{}*t)+{}".format(optimizedParameters[0], optimizedParameters[1], optimizedParameters[2]));
  plt.xlabel('Seconds')

boat_vel = 5.7
dt = 0.1
ax = np.empty(2200)
ay = np.empty(2200)
i=0
for time in np.arange(0,220,0.1):
    print(time, ",", boat_vel)
    ax[i] = time
    ay[i] = boat_vel
    i += 1
    drag = boat.get_total_drag(boat_vel)
    mass = 2* boat.Dc * boat.Rho
    added_mass = 0
    acceleration = (drag * 1000) / (mass+added_mass)
    boat_vel -= (acceleration * dt)

#ax3.plot(ax, ay , 'g',  label="model" )
#ax3.legend(loc='upper right')

if args.output_data:
  data_file.close()


plt.title(args.title + "\nData from " + str(from_datetime) + " to " + str(to_datetime) + " delta: " + str(total_seconds) + " sec")
#plt.savefig('image.png', dpi=500)
plt.show()
