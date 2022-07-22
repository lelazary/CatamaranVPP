#!/usr/bin/python3

import numpy as np
import math

from modules.Boat import Boat
from modules.Utils import * 


boat = Boat(
  Name = "Outremer51",
  Lwl = 15.41,   #Length of waterline in meters
  Bwl = 1.24,    #Beam of each hull
  Tc = 0.90,     #Maximum draft of haull body
  i = 5.44,      #half entry angle of the waterline at bow (neglecting the local rounded shape at the stem)
  Sw = 29.05,    #wetted area of one hull
  Dc = 6.8779,   #Displacment on one hull in m3
  Cp = 0.573,    #hull prismatic coefficient
  At = 0.00,     #immersed part of the transom area at zero speed
  LCB = 46.4,    # Longitudinal center of buoyancy, in % of Lwl from aft perpendicular
  S = 5.72,      #Space between hulls, transversal space between the 2 hulls axis
  S_aero = 11.35,  #Frontal Secion 
  Cx_aero = 0.40, #he Cx to take into account for the wind drag comptutation (to estimate due to the overall shape of the superstructure)
  Prop_eff = 0.54,  #to estimate with the one for the propeller and the one of the mechanical transmission. Example : 0.54 = 0.6 (propeller) x 0.9 (mechaniccal transmission)
  Cms = 0.70,     # midship section coefficient (section / Tc Bw)
  Cwp = 0.70      #waterplane area coefficient (Sf / Lw Bw)
)
 

#print(boat.get_total_drag(3.0))
Fn04 = 0.4*math.sqrt(9.81*boat.Lwl)*36/18.52
Fn03 = 0.3*math.sqrt(9.81*boat.Lwl)*36/18.52
Fn1 = 1*math.sqrt(9.81*boat.Lwl)*36/18.52
print(Fn03)
speed = Fn04
for i in range(0,15):
    Vb_ms = kts2ms(speed)
    print(speed, Vb_ms, boat.get_total_drag(Vb_ms))
    speed  += (Fn1-Fn04)/14
