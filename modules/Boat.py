#!/usr/bin/python

#Created by Lior Elazary
#liorelazary@gmail.com


import numpy as np
import math

class Boat(object):
    def __init__(self, Name, Lwl, Bwl, Tc, i, Sw, Dc, Cp,
                 At, LCB, S, S_aero, Cx_aero, Prop_eff, Cms, Cwp):

        self.Name = Name #Boat name for reports and comperison
        self.Lwl = Lwl   #Length of waterline in meters
        self.Bwl = Bwl    #Beam of each hull
        self.Tc = Tc     #Maximum draft of haull body
        self.i = i      #half entry angle of the waterline at bow (neglecting the local rounded shape at the stem)
        self.Sw = Sw    #wetted area of one hull
        self.Dc = Dc   #Displacment on one hull in m3
        self.Cp = Cp    #hull prismatic coefficient
        self.At = At     #immersed part of the transom area at zero speed
        self.LCB = LCB    # Longitudinal center of buoyancy, in % of Lwl from aft perpendicular
        self.S = S      #Space between hulls, transversal space between the 2 hulls axis
        self.S_aero = S_aero  #Frontal Secion 
        self.Cx_aero = Cx_aero #he Cx to take into account for the wind drag comptutation (to estimate due to the overall shape of the superstructure)
        self.Prop_eff = Prop_eff  #to estimate with the one for the propeller and the one of the mechanical transmission. Example : 0.54 = 0.6 (propeller) x 0.9 (mechaniccal transmission)
        self.Cms = Cms     # midship section coefficient (section / Tc Bw)
        self.Cwp = Cwp      #waterplane area coefficient (Sf / Lw Bw)

        self.Lw_Bw_ratio = self.Lwl/self.Bwl
        self.Lw_Dc = self.Lwl/math.pow(self.Dc, 1/3)
        self.Bw_Tc_ratio = self.Bwl/self.Tc
        self.Lw_Tc_ratio = self.Lwl/self.Tc
        self.S_Lw_ratio = self.S/self.Lwl
        self.Vb_ms = 0
        

        self.Rho = 1025 #Salt water density kg/m3
        self.v_1 = 0.000001  #Viscosity
        self.gravity = 9.81


    def set_vel(self, Vb_ms): #Need to set that for the velocity of the vessel
        self.Fn = Vb_ms/math.sqrt(self.gravity*self.Lwl) #Froude number
        self.Vb_ms = Vb_ms

    def get_total_drag(self, Vb_ms):

        self.set_vel(Vb_ms)

        self.Df = self.get_frictional_drag(Vb_ms)
        self.Dr = self.get_wave_drag(Vb_ms)
        self.Dr_aero = self.get_aero_drag(Vb_ms)
        self.Dtr = 0  #TODO: is the rear transoms drag
        
        self.D_total = self.Df + self.Dr + self.Dtr + self.Dr_aero  #Total Drag

        #print(self.Df, self.Dr, self.Dr_aero, self.D_total)
        
        return self.D_total

    def get_aero_drag(self, Vb_ms):

        R_air = 1.225  #Drag due to aerodynamical kg/m3
        Dr_aero =(0.5*R_air*math.pow(Vb_ms,2)*self.Cx_aero*self.S_aero)/1000

        return Dr_aero

    def get_frictional_drag(self, Vb_ms):
        #Df -- Hulls frictional Drag 

        Fr = 0.5 * self.Rho * math.pow(Vb_ms,2)
        #Compute the Reynolds number VL/v (v=viscosity)
        Re_hull = Vb_ms * (0.7 * self.Lwl)/self.v_1
        
        if Re_hull == 0:  #Handle the 0 case
            return 0
        Cf_hull = 0.075/math.pow((math.log10(Re_hull)-2),2)
        Df = Fr*Cf_hull*(2*self.Sw)/1000   #in kN
          
        return Df

    def get_wave_drag(self, Vb_ms):
    
        ##Compute Hull residuary drag Dr (Wave drag)
        #, where K is the interhulls factor due to waves interference :
        #Dr cata = (1+K) (Dr of the 2 hulls)
        
        Eq_1_a = 28.52563292378
        Eq_1_n = 3.0289127825466600
        
        Eq_2_a = -29.79735956867010
        Eq_2_b = 38.4738489436956
        Eq_2_c = -8.8720308419296900000
        
        Eq_3_a = 5.16744889765720000000
        Eq_3_b = 0.38475974097213300000
        
        #at 0.4 Froude number, the two waves at the bow and stern wave 
        #Cause the hull to drag sugnificatly (hull speed)
        Fn_hull_speed_ratio = 0.4    
        
        if (self.Fn < 0.4):
          Fn_1 = Eq_1_a*math.pow(self.Fn,Eq_1_n)
        elif(self.Fn >= 0.4 and self.Fn < 0.6):
          Fn_1 = Eq_2_a*math.pow(self.Fn,2)+Eq_2_b*self.Fn+Eq_2_c
        else:  #>0.6
          Fn_1 = Eq_3_a*self.Fn+Eq_3_b
        
        
        if self.Fn < 0.25:
          Dr_mg_ratio = self.get_Holtrop_Rw_A(Vb_ms)
        elif self.Fn>0.4: 
          Dr_mg_ratio = Fn_1
        else:
          holtrop_rw_a = self.get_Holtrop_Rw_A(Vb_ms)
          Dr_mg_ratio = (self.Fn-0.25)/0.15*Fn_1 +(0.4-self.Fn)/0.15*holtrop_rw_a
        
        #print(Dr_mg_ratio)
        
        #Interhull interaction
        Coeff_bis = max(0,0.45-0.05/0.2*(self.S_Lw_ratio-0.3))
        du_pic = min(1, 0.65+0.35*math.pow((10.5-self.Lw_Bw_ratio)/(10.5-7),2))
        du_corr = min(1, 0.3 + 0.6*(self.Lw_Dc-6.27)/(8.5-6.27))
        
        Coeff_ter_1 = 0.4*(du_corr*0.4)/(du_pic*0.6)*1.5
        Coeff_ter_2 = Coeff_ter_1 * 1.5
        
        #print("D:", Coeff_bis, du_pic, du_corr, Coeff_ter_1, Coeff_ter_2, self.S_Lw_ratio)
        #K_inter = math.max(
        d =  ( Coeff_ter_1+(0.6-Coeff_ter_1) 
              * math.exp(-0.5*(145.0-150.0*self.S_Lw_ratio)*math.pow((self.Fn-Coeff_bis),2))
             ) * (-2*self.S_Lw_ratio+1.6)*du_pic
        d1 =  ( Coeff_ter_2+(0.6-Coeff_ter_2) 
              * math.exp(-0.5*(145.0-150.0*self.S_Lw_ratio)*math.pow((self.Fn-Coeff_bis),2))
             ) * (-2*self.S_Lw_ratio+1.6)*du_pic
        
        if (self.Fn > Coeff_bis):
          K = d
        else:
          K = d1
 
        #print("K=", K, Dr_mg_ratio)
        Dr = (1+K)*(Dr_mg_ratio/100)*(2*self.Dc*self.Rho)*9.81/1000

        return Dr
    
    def get_Holtrop_Rw_A(self, Vb_ms):

        if Vb_ms == 0:
            return 0

        LwBw = self.Lw_Bw_ratio
        BwTc = self.Bw_Tc_ratio

        if (self.Lw_Bw_ratio < 12):
            lambda_1 = 1.446*self.Cp-0.3*self.Lw_Bw_ratio
        else:
            lambda_1 = 1.446*self.Cp-0.36

        if (self.Cp<0.8):
            c16 = 8.07981*self.Cp-13.8673*math.pow(self.Cp,2)+6.984388*math.pow(self.Cp,3)
        else:
            c16 = 1.73014-0.7067*self.Cp

        m1 = 0.0140407*(self.Lw_Tc_ratio)-1.75254/self.Lw_Dc-4.79323/self.Lw_Bw_ratio-c16

        LwDc_ratio = math.pow(self.Lw_Dc,3)

        if (LwDc_ratio < 512):
            c15 = -1.69385
        elif LwDc_ratio >1727:
            c15 = 0
        else:
            c15 = -1.69385+(self.Lw_Dc-8)/2.36

        d = -0.9

        c2 = 1.0
        Cm = self.Cms
        c5 = 1-0.8*self.At/(self.Bwl*self.Dc*self.Cms)

        c7 = 0.0991877700 #TODO, compute from eq

        c1 = 2223105*math.pow(c7,3.78613)*math.pow(1/self.Bw_Tc_ratio,1.07961)*math.pow(90-self.i,-1.37565)

        #print("Vb:", Vb_ms, "Fn:", self.Fn)
        #print(lambda_1, c16, m1, LwDc_ratio, c15, d, c2, Cm, c5, c7, c1)

        m4 = c15*0.4*math.exp(-0.034*math.pow(self.Fn,-3.29))

        #m1Fn_d = math.pow(self.Fn*m1, d)
        #m4cos = m4 * math.cos(lambda_1 * math.pow(self.Fn,-2))
        #math.exp(m1Fn_d + m4cos)a

        return c1*c2*c5*math.exp(m1*math.pow(self.Fn,d)+m4*math.cos(lambda_1*math.pow(self.Fn,-2)))*100

