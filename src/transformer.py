
import numpy as np
import pyproj


class Transformer:
    def __init__( self, 
                  lats: np.ndarray,
                  lons: np.ndarray,
                  alts: np.ndarray):        
        
        # Center of Region of Interest
        min_lat = np.min(lats)
        min_lon = np.min(lons) 
        min_alt = np.min(alts)
        max_lat = np.max(lats)
        max_lon = np.max(lons) 
        max_alt = np.max(alts) 
        med_lat = np.median(lats)
        med_lon = np.median(lons)     
        self.center_gps = (med_lat, med_lon, min_alt)
        
        # Transformer between coodinates
        aoi = pyproj.aoi.AreaOfInterest( west_lon_degree  = min_lon,
                                         south_lat_degree = min_lat,
                                         east_lon_degree  = max_lon,
                                         north_lat_degree = max_lat )
        utm_crs_list = pyproj.database.query_utm_crs_info( datum_name="WGS 84", area_of_interest=aoi)
        wgs84_crs = 'epsg:4326'
        utm_crs   = pyproj.CRS.from_epsg(utm_crs_list[0].code)  #f'epsg:{utm_crs_list[0].code}'
        
        self.gps_to_utm_transformer = pyproj.Transformer.from_crs(wgs84_crs, utm_crs, always_xy=False, only_best=True)
        
        # Center for UTM to XYZ
        xpj, ypj, zpj = self.gps_to_utm_transformer.transform(self.center_gps[0], self.center_gps[1], self.center_gps[2], direction='FORWARD')
        self.center_utm = np.asarray([xpj, ypj, zpj])
        
        # Ellipsoid Model for Distances
        self.geod = pyproj.Geod(ellps='WGS84')    
        
        
    def gps_to_utm( self, 
                    lats: np.ndarray,
                    lons: np.ndarray,
                    alts: np.ndarray):
        xpj, ypj, zpj = self.gps_to_utm_transformer.transform(lats, lons, alts, direction='FORWARD')
        return xpj, ypj, zpj
    
    def utm_to_xyz( self,
                    utm_x: np.ndarray,
                    utm_y: np.ndarray,
                    utm_z: np.ndarray ):
        x = utm_x - self.center_utm[0]
        y = utm_y - self.center_utm[1]
        z = utm_z - self.center_utm[2]
        return x, y, z
    
    def gps_to_xyz( self, 
                    lats: np.ndarray,
                    lons: np.ndarray,
                    alts: np.ndarray ):
        xpj, ypj, zpj = self.gps_to_utm(lats, lons, alts)
        x, y, z =  self.utm_to_xyz(xpj, ypj, zpj)
        return x, y, z
        
    def utm_to_gps( self, 
                    points: np.ndarray ):
        if points.ndim == 1:
            utm_points = points[np.newaxis,:]
        elif points.ndim == 2:
            utm_points = points  
        gps = self.gps_to_utm_transformer.transform(utm_points[:,0], utm_points[:,1], direction='INVERSE')
        return np.vstack( [gps[0], gps[1]] ).T
    
    def xyz_to_utm( self,
                    points: np.ndarray ):
        if points.ndim == 1:
            xyz_points = points[np.newaxis,:]
        elif points.ndim == 2:
            xyz_points = points
        utm_points = np.zeros( xyz_points.shape )
        for i in range(0, xyz_points.shape[0]):
            utm_points[i,:] = xyz_points[i,:] + self.center_utm
        return utm_points
        
    def xyz_to_gps( self,
                    points_xyz: np.ndarray ):
        utm = self.xyz_to_utm(points_xyz)
        gps = self.utm_to_gps(utm)
        return gps
    
    def gps_distance( self, 
                      from_lats: np.ndarray,
                      from_lons: np.ndarray,
                      to_lats:   np.ndarray,
                      to_lons:   np.ndarray ):
        az1, az2, dist_2d = self.geod.inv(from_lons, from_lats, to_lons, to_lats)
        return dist_2d
        
    def get_path_distance( self, 
                           lats: np.ndarray,
                           lons: np.ndarray ):
        dist_2d = self.gps_distance(lats[:-1], lons[:-1], lats[1:], lons[1:])
        return np.concatenate( ([0.0], dist_2d) )
    
    def get_speed_distance( self,
                            time:     np.ndarray,
                            velocity: np.ndarray ):
        dtime = np.diff(time)
        avelo = 0.5 * (velocity[1:] + velocity[:-1])
        ddist = dtime * avelo
        dist  = np.concatenate( ([0.0], np.cumsum(ddist)) )
        return dist
