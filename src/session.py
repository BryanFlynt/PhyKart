

from mychron_csv import MychronCSV
import numpy as np
import pathlib
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

from transformer import Transformer

from scipy.spatial import ConvexHull

class Session:
    
    def __init__(self, 
                 mychron: MychronCSV):
        self.data        = {}    # Map of data fields
        self.transformer = None  # Coordinate Transformer
        self.lap_times   = []    # Lap Times (seconds)
               
        if "Beacon Markers" in mychron.header:
            marks = mychron.header["Beacon Markers"].split(',')
            marks.insert(0, '0.0')
            self.lap_times = np.diff( np.array(marks, dtype=np.float64) )        
        
        orig_key  = 'Time'
        if orig_key in mychron.data:
            new_key   = 'gps_time'
            func = lambda x: x
            self.data[new_key] = func( mychron.data[orig_key] )
            
        orig_key  = 'GPS Latitude'
        if orig_key in mychron.data:
            new_key   = 'gps_latitude'
            func = lambda x: x   # degrees
            self.data[new_key] = func( mychron.data[orig_key] ) 
            
        orig_key = 'GPS Longitude'
        if orig_key in mychron.data:
            new_key  = 'gps_longitude'
            func = lambda x: x   # degrees
            self.data[new_key] = func( mychron.data[orig_key] )  
            
        orig_key  = 'GPS Altitude'
        if orig_key in mychron.data:
            new_key   = 'gps_altitude'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'ft':
                func = lambda x: 0.3048 * x   # Feet -> meters
            elif orig_unit == 'm':
                func = lambda x: x
            self.data[new_key] = func( mychron.data[orig_key] ) 
            
        orig_key  = 'GPS Speed'
        if orig_key in mychron.data:
            new_key   = 'gps_velocity'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'mph':
                func = lambda x: 0.44704 * x      # mph -> meters/sec
            elif orig_unit == 'km/h':
                func = lambda x: x * 1000/3600.0  # km/hr -> meters/sec
            self.data[new_key] = func( mychron.data[orig_key] )
            
        orig_key = 'GPS LatAcc'
        if orig_key in mychron.data:
            new_key   = 'gps_lateral_acc'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'g':
                func = lambda x: 9.80665 * x  # g -> to meters/sec^2
            self.data[new_key] = func( mychron.data[orig_key] )
            
        orig_key = 'GPS LonAcc'    
        if orig_key in mychron.data:
            new_key  = 'gps_longitudinal_acc'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'g':
                func = lambda x: 9.80665 * x  # g -> to meters/sec^2
            self.data[new_key] = func( mychron.data[orig_key] ) 
               
        orig_key = 'GPS Slope'
        if orig_key in mychron.data:
            new_key  = 'gps_slope'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'deg':
                func = lambda x: x * np.pi / 180.0 # Deg -> Radians
            self.data[new_key] = func( mychron.data[orig_key] )             
            
        orig_key = 'RPM'
        if orig_key in mychron.data:
            new_key  = 'engine_rotation_rate'
            orig_unit =  mychron.units[orig_key]
            if orig_unit == 'rpm':
                func = lambda x: 2*np.pi/60.0 * x   # rpm -> radians/sec
            self.data[new_key] = func( mychron.data[orig_key] ) 
  
        # Build UTM & XYZ Coordinates
        if 'gps_latitude' in self.data:
            if 'gps_longitude' in self.data:
                if 'gps_altitude' in self.data:
                    self.transformer = Transformer(self.data['gps_latitude'],
                                                   self.data['gps_longitude'],
                                                   self.data['gps_altitude'])
                    xpj, ypj, zpj = self.transformer.gps_to_utm(self.data['gps_latitude'],
                                                                self.data['gps_longitude'],
                                                                self.data['gps_altitude'])
                    x, y, z = self.transformer.utm_to_xyz(xpj, ypj, zpj)
                    self.data['utm_x'] = xpj
                    self.data['utm_y'] = ypj
                    self.data['utm_z'] = zpj
                    self.data['x'] = x
                    self.data['y'] = y
                    self.data['z'] = z
                    
        # Build Total Acceleration
        if 'gps_lateral_acc' in self.data:
            if 'gps_longitudinal_acc' in self.data:
                alat = self.data['gps_lateral_acc']
                alon = self.data['gps_longitudinal_acc']
                self.data['total_acc'] = np.sqrt(alat**2 + alon**2)
                    
        # Build Distance
        if 'gps_latitude' in self.data:
            if 'gps_longitude' in self.data:
                ddist = self.transformer.get_path_distance(self.data['gps_latitude'], self.data['gps_longitude'])
                self.data['gps_distance'] = np.cumsum(ddist)
        #if 'gps_time' in self.data:
        #    if 'gps_velocity' in self.data:
        #        self.data['gps_speed_distance'] = self.transformer.get_speed_distance(self.data['gps_time'], self.data['gps_velocity'])
        
        # Add lap number
        if 'gps_time' in self.data:
            self.data['lap'] = np.zeros_like(self.data['gps_time'])
            stime = 0.0
            for i in range(0,len(self.lap_times)):
                etime = stime + self.lap_times[i]
                mask = (stime <= self.data['gps_time']) & (self.data['gps_time'] <= etime)
                self.data['lap'][mask] = i
                stime = etime
       
    def num_laps(self):
        return len(self.lap_times)
    
    
    def beacon_location(self, system='gps'):
        """
        Find the location that cooresponds to the lap time start/stop times.

        Parameters
        ----------
        system : TYPE, optional
            DESCRIPTION. Coordinate system to return results in.

        Returns
        -------
        None.

        """
        if len(self.lap_times) == 0:
            return None
            
        if system.lower()   == 'gps':
            cs_x = CubicSpline(self.data['gps_time'], self.data['gps_latitude'])
            cs_y = CubicSpline(self.data['gps_time'], self.data['gps_longitude'])
            cs_z = CubicSpline(self.data['gps_time'], self.data['gps_altitude'])
        elif system.lower() == 'utm':
            cs_x = CubicSpline(self.data['gps_time'], self.data['utm_x'])
            cs_y = CubicSpline(self.data['gps_time'], self.data['utm_y'])
            cs_z = CubicSpline(self.data['gps_time'], self.data['utm_z'])
        elif system.lower() == 'xyz':
            cs_x = CubicSpline(self.data['gps_time'], self.data['x'])
            cs_y = CubicSpline(self.data['gps_time'], self.data['y'])
            cs_z = CubicSpline(self.data['gps_time'], self.data['z'])
        else:
            print(f'ERROR: Coordinate system {system} not recognized')
        
        sum_times = np.cumsum(self.lap_times)
            
        # If 1 laptime then pick it
        # If 2 lap times pick the last
        # If multiple then drop 1st and last
        if len(self.lap_times) == 1:
            times = sum_times[0]
        elif len(self.lap_times) == 2:
            times = sum_times[-1]
        elif len(self.lap_times) > 2:
            times = sum_times[1:-1]
            
        x = []
        y = []
        z = []
        for t in times:
            x.append( cs_x(t) )
            y.append( cs_y(t) )
            z.append( cs_z(t) )
        return np.array([np.average(x), np.average(y), np.average(z)])        

    def get_intersections(self, xx, yy, system='gps', parameter='gps_time', tol=5.0):
        """
        Find the values of the parameter where the trajectory has it's closest
        approach to the provided point location
        
        Parameters
        ----------
        xx : Float
            Latitude, UTM_X, X Locations of intersection
        yy : Float
            Longitude, UTM_Y, Y Locations of intersection
        system : Float, optional
            The coodinate system of the provided point.
        parameter : Float, optional
            The parameter to approximate at intersection point.
        tol : Float, optional
            Closest approach must be within this tolerance to be counted
        
        Returns
        -------
        value of parameter at intersection point.
        
        """
        intersect_times = []
        if system.lower()   == 'gps':
            xp, yp = self.transform.gps_to_xyz(xx,yy)
        elif system.lower() == 'utm':
            xp, yp = self.transform.utm_to_xyz(xx,yy)
        elif system.lower() == 'xyz':
            xp, yp = xx, yy
        else:
            print(f'ERROR: Coordinate system {system} not recognized')
        
        # Calculate distance^2 to each point in time
        dist2 = np.zeros_like(self.data['gps_time'])
        for i, xy in enumerate(zip(self.data['x'],self.data['y'])):
            dist2[i] = (xy[0]-xp)**2 + (xy[1]-yp)**2
                
        # Build a list of entry/exit times from tol of point
        tol2  = tol**2
        stime = self.data['gps_time'][0]
        time_pairs = []
        for i in range(1,len(dist2)):
            if (dist2[i-1] > tol2) and (dist2[i] <= tol2):     # Now inside
                stime = self.data['gps_time'][i]
            elif (dist2[i-1] <= tol2) and (dist2[i] > tol2):   # Now outside
                time_pairs.append((stime, self.data['gps_time'][i]))
            pass
        if dist2[-1] <= tol2:
            time_pairs.append((stime, self.data['gps_time'][-1]))
        
        # Search for minimum within each interval
        ans = []
        cs  = CubicSpline(self.data['gps_time'], dist2)        
        for stime, etime in time_pairs:
            res = minimize_scalar(cs, bounds=(stime,etime))
            ans.append(res.x)
            
        return np.array(ans)
    
    
    def get_lap(self,num):
        ans = {}
        mask = self.data["lap"] == num
        for key, valarray in self.data.items():
            ans[key] = valarray[mask]
        return ans