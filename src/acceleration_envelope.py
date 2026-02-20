#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import sys

class AccelerationEnvelope:
    def __init__(self, session, mirror=True):
        velo = session.data['gps_velocity']
        alat = session.data['gps_lateral_acc']
        alon = session.data['gps_longitudinal_acc']
        if mirror:
            velo = np.append(velo, velo)
            alat = np.append(alat, -1.0*alat)
            alon = np.append(alon, alon)
        self.pts  = np.array([velo, alat, alon]).T
        self.hull = ConvexHull(self.pts)
        
        
    def get_accel(self, velo=None, alat=None, alon=None):
        """
        Find intersection of segment and triangle in 3D space following
        
        https://stackoverflow.com/questions/54143142/3d-intersection-between-segment-and-triangle

        Parameters
        ----------
        velo : TYPE, optional
            DESCRIPTION. Velocity
        alat : TYPE, optional
            DESCRIPTION. Lateral acceleration will search for longitudinal acc
        alon : TYPE, optional
            DESCRIPTION. Longitudinal acceleration will search for Lateral acc

        Returns
        -------
        ans : TYPE
            DESCRIPTION.

        """
        if velo is None:
            s1 = np.array([-100., alat, alon], dtype=np.float64)
            s2 = np.array([+100., alat, alon], dtype=np.float64)
            ix = 0         
        elif alat is None:
            s1 = np.array([velo, -100., alon], dtype=np.float64)
            s2 = np.array([velo, +100., alon], dtype=np.float64)
            ix = 1
        elif alon is None:
            s1 = np.array([velo, alat, -100.], dtype=np.float64)
            s2 = np.array([velo, alat, +100.], dtype=np.float64)
            ix = 2    
        
        ans_list = []
        for s in self.hull.simplices:
            v = self.pts[s, :]
            n = np.cross(v[1,:]-v[0,:], v[2,:]-v[0,:])
            n = n / np.linalg.norm(n, ord=2)
            d = -1.0 * np.dot(n, v[0,:])
            den = np.dot(n, (s2-s1))
            if( den == 0.0 ):
                continue
            t = -1.0 * (d + np.dot(n, s1)) / den
            p = s1 + t * (s2 - s1)
            
            n01 = np.cross(v[1,:] - v[0,:], n)
            d01 = np.dot(p - v[0,:], n01) / np.linalg.norm(n01, ord=2)
            if d01 > 0.0:
                continue
            
            n12 = np.cross(v[2,:] - v[1,:], n)
            d12 = np.dot(p - v[1,:], n12) / np.linalg.norm(n12, ord=2)
            if d12 > 0.0:
                continue
            
            n20 = np.cross(v[0,:] - v[2,:], n)
            d20 = np.dot(p - v[2,:], n20) / np.linalg.norm(n20, ord=2)
            if d20 > 0.0:
                continue
            
            ans_list.append(p[ix])

        if len(ans_list) == 0:
            return []
        return [min(ans_list), max(ans_list)]          
    
    
    def plot_gg(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')# Plot defining corner points
    
        ax.plot(self.pts.T[1], self.pts.T[2], self.pts.T[0], "ko")
    
#        for s in self.hull.simplices:
#            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
#            ax.plot(self.pts[s, 1], self.pts[s, 2], self.pts[s, 0], "r-")
        
        ax.set_xlabel('lateral_acc (m/s^2)')
        ax.set_ylabel('longitudinal_acc (m/s^2)')
        ax.set_zlabel('velocity (m/s)')
        plt.show()
        
    def plot_tg(self):
        fig = plt.figure()
        ax = fig.add_subplot()  # Plot defining corner points

        atot = np.sqrt(self.pts.T[1]**2 + self.pts.T[2]**2)
        ax.plot(self.pts.T[0], atot, "ko")
        ax.set_xlabel('Velocity (m/s)')
        ax.set_ylabel('Total Acceleration (m/s^2)')
        plt.grid()
        plt.show()

        