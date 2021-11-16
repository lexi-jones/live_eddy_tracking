# transect_visualization.py
# Functions to visualize the eddy field and cruise transect. Use interactive visualizations with visualize_eddy_field_and_transect.ipynb Jupyter notebook. 

# Lexi Jones
# Date Created: 11/14/21
# Last Edited: 11/15/21

import math,datetime
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import shapely.geometry
import matplotlib.pyplot as plt

def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Returns the distance in kilometers between two coordinate points. Accepts negative (-180 to 180) or positive coordinate systems (0 to 360). 

    lat1,lon1: coordinates for point 1
    lat2,lon2: coordinates for point 2
    """

    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = math.radians(lat2-lat1),math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    
    a = math.sin(delta_lat/2) * math.sin(delta_lat/2) + math.cos(lat1_radians) * math.cos(lat2_radians) * ((math.sin(delta_lon/2))**2)
    dist = R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) # Distance in km
    
    return dist

def bearing_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Returns the bearing between two points.
    lat1,lon1: coordinates for point 1
    lat2,lon2: coordinates for point 2
    """
    
    delta_lon = math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    
    X = math.cos(lat2_radians) * math.sin(delta_lon)
    Y = math.cos(lat1_radians) * math.sin(lat2_radians) - (math.sin(lat1_radians) * math.cos(lat2_radians) * math.cos(delta_lon))

    return math.degrees(math.atan2(X,Y))

def find_lat_lon_end_point(lat1,lon1,dist,bearing):
    """
    Returns a latitude and longitude in degrees, the coordinate point some distance away from a starting point, travelling at a 'bearing' angle.

    lat1,lon1: starting coordinates
    dist: distance travelled in km
    bearing: angle clockwise from north (0 to 360) 
    """

    R = 6371
    bearing_radians = math.radians(bearing)
    lat1_radians,lon1_radians = math.radians(lat1),math.radians(lon1)
    
    lat2_radians = math.asin(math.sin(lat1_radians) * math.cos(dist/R) + math.cos(lat1_radians) * math.sin(dist/R) * math.cos(bearing_radians))
    lon2_radians = lon1_radians + math.atan2(math.sin(bearing_radians) * math.sin(dist/R) * math.cos(lat1_radians), math.cos(dist/R) - math.sin(lat1_radians)*math.sin(lat2_radians))
    
    return math.degrees(lat2_radians),math.degrees(lon2_radians)

def plot_eddy_field_transect(date,station_lats,station_lons,anti_eddy_data,cyc_eddy_data,cons,plms,LAVD_flag):
    """
    Plot SSH eddies, RCLVs, and expected cruise track. 

    date: Format 'yyyymmdd'
    station_lats,station_lons: array of the coordinates of the stations
    anti_eddy_data,cyc_eddy_data: SSH eddy data
    cons,plms: RCLV contours and local maxima from LAVD field
    LAVD_flag: 1 to plot LAVD in the background, 0 to not 

    """

    fig,ax = plt.subplots(1,1,figsize=(20,15))
    
    # PREP SSH DATA
    anti_day_inds = np.where(anti_eddy_data[:,0] == str(date))
    cyc_day_inds = np.where(cyc_eddy_data[:,0] == str(date))

    anti_eddy_bnds = np.squeeze(anti_eddy_data[anti_day_inds,9:])
    cyc_eddy_bnds = np.squeeze(cyc_eddy_data[cyc_day_inds,9:])
        
    anti_eddy_extremum_lon = np.squeeze(anti_eddy_data[anti_day_inds,5])
    anti_eddy_extremum_lat = np.squeeze(anti_eddy_data[anti_day_inds,6])
    cyc_eddy_extremum_lon = np.squeeze(cyc_eddy_data[cyc_day_inds,5])
    cyc_eddy_extremum_lat = np.squeeze(cyc_eddy_data[cyc_day_inds,6])

    #### PLOT LAVD IN BACKGROUND IF LAVD_FLAG == 1 ####
    date_dashed = '%s-%s-%s'%(date[0:4],date[4:6],date[6:8])
    lon_array = np.arange(-150+360,-117+360,0.03125)
    lat_array = np.arange(10,33,0.03125)

    if (LAVD_flag == 1):
        LAVD_dir = './LAVD/'
        LAVD_file_name = '%s_LAVD_8days_backward_runtime_20min_timestep_particle_start_lat_10_33_lon_-150_-117_spatial_step_0.03125_6hr_output_freq.npy'%(date_dashed)
        LAVD = np.load(LAVD_dir + LAVD_file_name)
        LAVD = np.ma.masked_where(np.isnan(LAVD),LAVD)
        LAVD_reshape = np.transpose(np.reshape(LAVD,(len(lon_array),len(lat_array))))
        
        ax.pcolormesh(lon_array,lat_array,LAVD_reshape,cmap='inferno',vmax=0.00001) 

    #### PLOT ANTICYCLONES ####
    i = 0 
    label_flag = 0
    for bnds in anti_eddy_bnds:
        x_values = [float(coord) for coord in bnds[0::2] if str(coord) != '']
        y_values = [float(coord) for coord in bnds[1::2] if str(coord) != '']
        
        # Crop x-axis
        if (np.max(x_values) < 210):
            i+=1 
            continue

        #Plot eddy bounds and extremum
        if LAVD_flag == 1:
            anti_color = 'w'
        else: 
            anti_color = '#991414'

        if label_flag == 0:
            ax.plot(x_values,y_values,c=anti_color,label='SSH Anticyclone')
            label_flag += 1
        else:
            ax.plot(x_values,y_values,c=anti_color)
        ax.scatter(float(anti_eddy_extremum_lon[i])+360,float(anti_eddy_extremum_lat[i]),color=anti_color)
            
        i+=1 

    #### PLOT CYCLONES ####
    j = 0 
    label_flag = 0
    for bnds in cyc_eddy_bnds:    
        x_values = [float(coord) for coord in bnds[0::2] if str(coord) != '']
        y_values = [float(coord) for coord in bnds[1::2] if str(coord) != '']
        
        # Crop x-axis
        if (np.max(x_values) < 210):
            j+=1
            continue
        
        #Plot eddy bounds and extremum
        if LAVD_flag == 1:
            cyc_color = 'w'
        else: 
            cyc_color = '#08088A'

        if label_flag == 0:    
            ax.plot(x_values,y_values,c=cyc_color,label='SSH Cyclone')
            label_flag += 1
        else:
            ax.plot(x_values,y_values,c=cyc_color)    
        ax.scatter(float(cyc_eddy_extremum_lon[j])+360,float(cyc_eddy_extremum_lat[j]),color=cyc_color)
            
        j+=1
                
    #### PLOT RCLVs and LAVD peaks ####
    label_flag = 0
    for con in cons:
        x_values = lon_array[[int(j) for j in con[:, 1]]]
        y_values = lat_array[[int(j) for j in con[:, 0]]]

        if (np.max(x_values) < 210):
            continue
        
        if (label_flag == 0):
            plt.plot(x_values,y_values,color='green',label='RCLV')
            label_flag += 1
        else:
            plt.plot(x_values,y_values,color='green')    

    for plm in plms:
        plm_lon = lon_array[plm[1]]
        plm_lat = lat_array[plm[0]]

        if (plm_lon < 210):
            continue    

        plt.scatter(plm_lon,plm_lat,color='green')
                        
    #### PLOT STATIONS ####
    if LAVD_flag == 1:
        station_color = 'w'
    else:
        station_color = 'k'

    if station_lons[0]<0:
        station_lons = [i+360 for i in station_lons]
    station_labels = ['San Diego','St. 1','St. 2','St. 3','St. 4','St. 5','St. 6']
    for k in np.arange(0,len(station_lons)):
        ax.scatter(station_lons[k],station_lats[k], c=station_color, marker=(5, 2), s=100)
        ax.text(station_lons[k]+0.75,station_lats[k]-0.2,station_labels[k],fontsize=18,color=station_color)        

        if k < len(station_lons)-1:
            ax.plot([station_lons[k],station_lons[k+1]],[station_lats[k],station_lats[k+1]],c=station_color)

    #### PLOT PARAMS ####
    ax.set_xlabel("Longitude",fontsize=24);
    ax.set_ylabel("Latitude",fontsize=24);
    ax.set_title("%s-%s-%s Eddy Field"%(date[4:6],date[6:8],date[0:4]), fontsize=26)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(fontsize=24,loc=4)

    return ax

def plot_eddy_field_transect_zoom(date,station_lats,station_lons,transect_lats,transect_lons,anti_eddy_data,cyc_eddy_data,cons,plms):
    if transect_lons[0]<0:
        transect_lons = [i+360 for i in transect_lons]

    fig,ax = plt.subplots(1,1,figsize=(7,7))
    
    # LOAD SSH DATA
    anti_day_inds = np.where(anti_eddy_data[:,0] == str(date))
    cyc_day_inds = np.where(cyc_eddy_data[:,0] == str(date))

    anti_eddy_bnds = np.squeeze(anti_eddy_data[anti_day_inds,9:])
    cyc_eddy_bnds = np.squeeze(cyc_eddy_data[cyc_day_inds,9:])
        
    anti_eddy_extremum_lon = np.squeeze(anti_eddy_data[anti_day_inds,5])
    anti_eddy_extremum_lat = np.squeeze(anti_eddy_data[anti_day_inds,6])
    cyc_eddy_extremum_lon = np.squeeze(cyc_eddy_data[cyc_day_inds,5])
    cyc_eddy_extremum_lat = np.squeeze(cyc_eddy_data[cyc_day_inds,6])

    #### PLOT ANTICYCLONES ####
    i = 0 
    for bnds in anti_eddy_bnds:
        x_values = [float(coord) for coord in bnds[0::2] if str(coord) != '']
        y_values = [float(coord) for coord in bnds[1::2] if str(coord) != '']

        min_lon,max_lon = np.min(x_values),np.max(x_values)
        min_lat,max_lat = np.min(y_values),np.max(y_values)
            
        #Crop the eddy field so we are only plotting with a 5 degree range from start & stop
        if (min_lon < transect_lons[1]-5) or (max_lon > transect_lons[0]+5) or (min_lat < transect_lats[1]-5) or (max_lat > transect_lats[0]+5):
            i+=1 
            continue
            
        ax.plot(x_values,y_values,c='#991414')
        ax.scatter(float(anti_eddy_extremum_lon[i])+360,float(anti_eddy_extremum_lat[i]),color='#991414')    
            
        i+=1 

    #### PLOT CYCLONES ####
    j = 0 
    for bnds in cyc_eddy_bnds:    
        x_values = [float(coord) for coord in bnds[0::2] if str(coord) != '']
        y_values = [float(coord) for coord in bnds[1::2] if str(coord) != '']
            
        min_lon,max_lon = np.min(x_values),np.max(x_values)
        min_lat,max_lat = np.min(y_values),np.max(y_values)

        if (min_lon < transect_lons[1]-3) or (max_lon > transect_lons[0]+3) or (min_lat < transect_lats[1]-3) or (max_lat > transect_lats[0]+3):
            j+=1
            continue

        ax.plot(x_values,y_values,c='#08088A')
        ax.scatter(float(cyc_eddy_extremum_lon[j])+360,float(cyc_eddy_extremum_lat[j]),color='#08088A')
            
        j+=1
        
    # LOAD RCLV DATA
    lon_array = np.arange(-150+360,-117+360,0.03125)
    lat_array = np.arange(10,33,0.03125)
        
    #### PLOT RCLVs and LAVD peaks ####
    for con in cons:
        x_values = lon_array[[int(j) for j in con[:, 1]]]
        y_values = lat_array[[int(j) for j in con[:, 0]]]
            
        min_lon,max_lon = np.min(x_values),np.max(x_values)
        min_lat,max_lat = np.min(y_values),np.max(y_values)

        if (min_lon < transect_lons[1]-3) or (max_lon > transect_lons[0]+3) or (min_lat < transect_lats[1]-3) or (max_lat > transect_lats[0]+3):
            continue

        ax.plot(x_values,y_values,color='green')
            
        for plm in plms:
            plm_lon = lon_array[plm[1]]
            plm_lat = lat_array[plm[0]]
            
        if (plm_lon < transect_lons[1]-3) or (plm_lon > transect_lons[0]+3) or (plm_lat < transect_lats[1]-3) or (plm_lat > transect_lats[0]+3):
            continue
            
        ax.scatter(plm_lon,plm_lat,marker=(5, 2),color='green')
            
            
    #### PLOT STATIONS ####
    station_labels = ['San Diego','St. 1','St. 2','St. 3','St. 4','St. 5','St. 6']
    if station_lons[0] < 360:
        station_lons = [i+360 for i in station_lons]
    for k in np.arange(0,len(station_lons)):
        station_lon = station_lons[k]
        station_lat = station_lats[k]
            
        if (station_lons[k] < transect_lons[1]-3) or (station_lons[k] > transect_lons[0]+3) or (station_lats[k] > transect_lats[0]+3) or (station_lats[k] < transect_lats[1]-3):
            continue
        ax.scatter(station_lon,station_lat,c='k', marker=(5, 2),s=100)        
        ax.text(station_lons[k]+0.5,station_lats[k]-0.1,station_labels[k],fontsize=18,color='k')

    #### PLOT TRANSECT ####
    ax.plot(transect_lons,transect_lats,c='k')
    ax.scatter(transect_lons,transect_lats,c='k')
    return ax

def get_locations_along_traj(date,start_hour,duration,station_lats,station_lons,station_arrivals,station_departures):
    """start date: 
       start hour: hour in the day (0 to 23)
       duration: number of hours
    """
    
    year,month,day = int(date[0:4]),int(date[4:6]),int(date[6:8])
    start_datetime = datetime.datetime(year,month,day,start_hour)
    
    ship_stop_starts = []
    for s in np.arange(0,len(station_arrivals)):
        ship_stop_starts.append(station_departures[s])
        ship_stop_starts.append(station_arrivals[s])
    
    transect_times = [start_datetime]
    for h in np.arange(0,duration):
        transect_times.append(start_datetime + datetime.timedelta(hours=int(h+1)))
    
    stop_datetime = transect_times[-1]
    
    print('Start time: %s'%(start_datetime))
    print('Stop time: %s\n'%(stop_datetime))
    
    transect_lats,transect_lons = [],[]
    for d in np.arange(0,len(ship_stop_starts)-1):
        if (start_datetime >= ship_stop_starts[d]) and (start_datetime < ship_stop_starts[d+1]):

            if (d%2 == 0): #in transit at start time
                
                #estimate lat/lon at this time
                past_station_lat,past_station_lon,past_station_time = station_lats[int(d/2)],station_lons[int(d/2)],ship_stop_starts[d]
                next_station_lat,next_station_lon,next_station_time = station_lats[int((d/2)+1)],station_lons[int((d/2)+1)],ship_stop_starts[d+1]
                                
                hrs_since_last_station = (start_datetime - past_station_time).total_seconds()/(60*60)
                hrs_between_stations = (next_station_time - past_station_time).total_seconds()/(60*60)
                
                start_lat = past_station_lat - ((past_station_lat - next_station_lat)/hrs_between_stations)*hrs_since_last_station
                start_lon = past_station_lon - ((past_station_lon - next_station_lon)/hrs_between_stations)*hrs_since_last_station
                
                transect_lats.append(start_lat)
                transect_lons.append(start_lon)
                
                print('Start during transit between station %s and %s at (%s,%s)'%(int(d/2),int((d/2)+1),round(start_lat,3),round(start_lon,3)))
                
                if stop_datetime < ship_stop_starts[d+1]: #still in transit 
                    hrs_since_last_station = (stop_datetime - past_station_time).total_seconds()/(60*60)
                    stop_lat = past_station_lat - ((past_station_lat - next_station_lat)/hrs_between_stations)*hrs_since_last_station
                    stop_lon = past_station_lon - ((past_station_lon - next_station_lon)/hrs_between_stations)*hrs_since_last_station
                
                    print('Stop during transit between station %s and %s at (%s,%s)'%(int(d/2),int((d/2)+1),round(stop_lat,3),round(stop_lon,3)))
            
                elif stop_datetime < ship_stop_starts[d+2]: #arrived at station
                    
                    print('Arrive at station %s (%s,%s) at %s'%(int((d/2)+1),next_station_lat,next_station_lon,ship_stop_starts[d+1]))
            
                    stop_lat = next_station_lat
                    stop_lon = next_station_lon
                    
                    print('Still at station %s at %s'%(int((d/2)+1),stop_datetime))
                    
                elif stop_datetime < ship_stop_starts[d+3]: #left new station
                    
                    print('Arrive at station %s (%s,%s) at %s'%(int((d/2)+1),next_station_lat,next_station_lon,ship_stop_starts[d+1]))
                    print('Leave station %s at %s'%(int((d/2)+1),ship_stop_starts[d+2]))

                    past_station_lat,past_station_lon,past_station_time = station_lats[int(d/2)+1],station_lons[int(d/2)+1],ship_stop_starts[d+2]
                    next_station_lat,next_station_lon,next_station_time = station_lats[int(d/2)+2],station_lons[int(d/2)+2],ship_stop_starts[d+3]
                
                    hrs_between_stations = (next_station_time - past_station_time).total_seconds()/(60*60)                                        
                    hrs_since_last_station = (stop_datetime - past_station_time).total_seconds()/(60*60)

                    stop_lat = past_station_lat - ((past_station_lat - next_station_lat)/hrs_between_stations)*hrs_since_last_station
                    stop_lon = past_station_lon - ((past_station_lon - next_station_lon)/hrs_between_stations)*hrs_since_last_station
                
                    print('Stop during transit between station %s and %s at (%s,%s)'%(int(d/2)+1,int((d/2)+2),round(stop_lat,3),round(stop_lon,3)))
                    
                else:
                    print('ERROR: Pick a shorter duration.')
                    
                transect_lats.append(stop_lat)
                transect_lons.append(stop_lon)
                                
            else: #on station at start time
                
                start_lat = station_lats[int((d+1)/2)]
                start_lon = station_lons[int((d+1)/2)]
                
                transect_lats.append(start_lat)
                transect_lons.append(start_lon)
                
                print('Start at station %s at (%s,%s)'%(int((d+1)/2),start_lat,start_lon))
                
                hours_left_on_station = (ship_stop_starts[d+1] - start_datetime).total_seconds()/(60*60)
                
                if hours_left_on_station >= duration:
                    print('Will not leave the station during this time')
                    
                    stop_lat,stop_lon = start_lat,start_lon

                    transect_lats.append(start_lat) #stop coords are same as start
                    transect_lons.append(start_lon)
                    
                elif hours_left_on_station < duration: 
                    print('Leave station %s at %s'%(int((d+1)/2),ship_stop_starts[d+1]))      
                    
                    remaining_time = duration - hours_left_on_station
                    
                    if stop_datetime < ship_stop_starts[d+2]: #still in transit 
                       #estimate lat/lon at this time

                        next_station_lat,next_station_lon,next_station_time = station_lats[int(((d+1)/2)+1)],station_lons[int(((d+1)/2)+1)],ship_stop_starts[d+2]
                        hrs_between_stations = (next_station_time - ship_stop_starts[d+1]).total_seconds()/(60*60)
                        
                        stop_lat = start_lat - ((start_lat - next_station_lat)/hrs_between_stations)*remaining_time
                        stop_lon = start_lon - ((start_lon - next_station_lon)/hrs_between_stations)*remaining_time

                        print('Stop during transit between station %s and %s at (%s,%s)'%(int((d+1)/2),int(((d+1)/2)+1),round(stop_lat,3),round(stop_lon,3)))

                    elif stop_datetime < ship_stop_starts[d+3]: #arrived at station

                        stop_lat = station_lats[int(((d+1)/2)+1)]
                        stop_lon = station_lons[int(((d+1)/2)+1)] 
                        
                        print('Arrive at station %s (%s,%s) at %s'%(int((d/2)+1),stop_lat,stop_lon,ship_stop_starts[d+1]))
                        print('Still at station %s at %s'%(int((d/2)+1),stop_datetime))   
                        
                    else:
                        print('ERROR: Pick a shorter duration.')
                        
                transect_lats.append(stop_lat)
                transect_lons.append(stop_lon)
    
    return transect_lats, transect_lons

def check_if_in_eddy(date,lon,lat,eddy_data):
    """
    Inputs
        date - Format YYYYMMDD
        lat - 10 to 70 (degrees north); latitude of the data point
        lon - input should be positive coordinates (110 to 250 degrees); longitude of the data point
        eddy_array - Name of the eddy dataset to look through (cyc or anti)
        
    Outputs
        0 if point is not in eddy, OR eddy ID if point is in an eddy
    """
    
    #Open .csv file with eddy data
    day_inds = np.where(eddy_data[:,0] == str(date)) #find all eddies from the requested date
    eddy_bnds = np.squeeze(eddy_data[day_inds,9:]) #get the eddy boundary coordinates & remove extraneous dimensions of the array
    eddy_ids = np.squeeze(eddy_data[day_inds,1]) #get the eddy ids
    
    #Iterate through the data for each eddy that is present on the requested date
    in_eddy = 0 #flag to determine if the point is in an eddy
    i = 0 #counter
    for eddy in eddy_bnds:
        bnds = [float(coord) for coord in eddy if str(coord) != ''] #Remove nan values from the eddy boundaries
        
        #Seperate the x & y values 
        x_values = bnds[0::2]
        y_values = bnds[1::2]
        
        #Check if the point is far away from the eddy we are checking right now, and if so, skip to the next eddy
        if lon > np.max(x_values) or lon < np.min(x_values) or lat > np.max(y_values) or lat < np.min(y_values):
            i += 1 #have to step up the counter here because the loop will skip the rest of the code here
            continue
        
        #If the point is close to the eddy, reformat the data so it is readable by the Polygon function
        poly_pts = [(x_values[pt],y_values[pt]) for pt in np.arange(0,len(x_values))]
        polygon = Polygon(poly_pts)
        data_pt = Point(lon,lat)
        
        if polygon.contains(data_pt):
            in_eddy = int(eddy_ids[i])
            break
                
        i += 1
    
    return in_eddy

def check_if_in_RCLV(date,lon,lat,cons):
    
    ## Open RCLV file
    lon_array = np.arange(-150+360,-117+360,0.03125)
    lat_array = np.arange(10,33,0.03125)

    #Iterate through the data for each eddy that is present on the requested date
    in_RCLV = 0 #flag to determine if the point is in an RCLV
    i = 0 #counter
    for con in cons:                  
        x_values = lon_array[[int(j) for j in con[:, 1]]]
        y_values = lat_array[[int(j) for j in con[:, 0]]]
        
        #Check if the point is far away from the eddy we are checking right now, and if so, skip to the next eddy
        if lon > np.max(x_values) or lon < np.min(x_values) or lat > np.max(y_values) or lat < np.min(y_values):
            i += 1 #have to step up the counter here because the loop will skip the rest of the code here
            continue
        
        #If the point is close to the eddy, reformat the data so it is readable by the Polygon function
        poly_pts = [(x_values[pt],y_values[pt]) for pt in np.arange(0,len(x_values))]
        polygon = Polygon(poly_pts)
        data_pt = Point(lon,lat)
        
        if polygon.contains(data_pt):
            in_RCLV = 1
            break
                
        i += 1
    
    return in_RCLV

def traj_eddy_intersection_pts(date,lon1,lat1,lon2,lat2,eddy_data):
    """
    Version 2 checks if the point in the middle of the line segment is in an eddy, rather than the border/edge points.
    
    Inputs
        date - Format YYYYMMDD
        lat - 10 to 70 (degrees north); latitude of the data point
        lon - input should be positive coordinates (110 to 250 degrees); longitude of the data point
        eddy_array - Name of the eddy dataset to look through (cyc or anti)
        
    Outputs
        0 if point is not in eddy, OR eddy ID if point is in an eddy
    """

    if lon1 < 0:
        lon1 = lon1 + 360
    if lon2 < 0:
        lon2 = lon2 + 360
    
    line = shapely.geometry.LineString([[lon1,lat1],[lon2,lat2]])
    
    #Open .csv file with eddy data
    day_inds = np.where(eddy_data[:,0] == str(date)) #find all eddies from the requested date
    eddy_bnds = np.squeeze(eddy_data[day_inds,9:]) #get the eddy boundary coordinates & remove extraneous dimensions of the array
    eddy_ids = np.squeeze(eddy_data[day_inds,1]) #get the eddy ids
    
    #Iterate through the data for each eddy that is present on the requested date
    eddy_id = []
    intersection_pts = []
    in_eddy = []
    i = 0 #counter
    for eddy in eddy_bnds:
        bnds = [float(coord) for coord in eddy if str(coord) != ''] #Remove nan values from the eddy boundaries
        
        #Seperate the x & y values 
        x_values = bnds[0::2]
        y_values = bnds[1::2]
        
        # Reformat the data so it is readable by the Polygon function
        poly_pts = [(x_values[pt],y_values[pt]) for pt in np.arange(0,len(x_values))]
        polygon = Polygon(poly_pts) 

        if line.intersects(polygon):
            try:
                pt_list = list(line.intersection(polygon).coords)

            except: #enter this loop if the transect intersects the same eddy multiple times (goes in, out, and back in again)
                pt_list = []    
                geoms = list(line.intersection(polygon).geoms)
                for g in np.arange(len(geoms)):
                    linestring = (list(line.intersection(polygon).geoms)[g])
                    pt_list = pt_list + list(linestring.coords)

            eddy_id.append(int(eddy_ids[i]))
            intersection_pts.append(pt_list)     
        
        i += 1
    
    return eddy_id,intersection_pts

def traj_RCLV_intersection_pts(date_dashed,lon1,lat1,lon2,lat2,cons):
    if lon1 < 0:
        lon1 = lon1 + 360
    if lon2 < 0:
        lon2 = lon2 + 360
    
    line = shapely.geometry.LineString([[lon1,lat1],[lon2,lat2]])
    
    ## Open RCLV file
    lon_array = np.arange(-150+360,-117+360,0.03125)
    lat_array = np.arange(10,33,0.03125)

    
    #Iterate through the data for each eddy that is present on the requested date
    intersection_pts,in_eddy = [],[]
    i = 0 #counter
    for con in cons:                  
        x_values = lon_array[[round(j) for j in con[:, 1]]]
        y_values = lat_array[[round(j) for j in con[:, 0]]]
        
        poly_pts = [(x_values[pt],y_values[pt]) for pt in np.arange(0,len(x_values))]
        polygon = Polygon(poly_pts) 
        
        # Reformat the data so it is readable by the Polygon function
        poly_pts = [(x_values[pt],y_values[pt]) for pt in np.arange(0,len(x_values))]
        polygon = Polygon(poly_pts) 

        if line.intersects(polygon):
            try:
                pt_list = list(line.intersection(polygon).coords)

            except: #enter this loop if the transect intersects the same eddy multiple times (goes in, out, and back in again)
                pt_list = []    
                geoms = list(line.intersection(polygon).geoms)
                for g in np.arange(len(geoms)):
                    linestring = (list(line.intersection(polygon).geoms)[g])
                    pt_list = pt_list + list(linestring.coords)
            intersection_pts.append(pt_list)     
        i += 1
    
    return intersection_pts

def transect_eddy_action(date,start_lat,start_lon,stop_lat,stop_lon,anti_eddy_data,cyc_eddy_data,cons):
    
    date_dashed = '%s-%s-%s'%(date[0:4],date[4:6],date[6:8])
    anti_IDs,cyc_IDs = [],[]        
    
    anti_start = check_if_in_eddy(date,start_lon,start_lat,anti_eddy_data)
    cyc_start = check_if_in_eddy(date,start_lon,start_lat,cyc_eddy_data)
    RCLV_start = check_if_in_RCLV(date_dashed,start_lon,start_lat,cons)
        
    anti_stop = check_if_in_eddy(date,stop_lon,stop_lat,anti_eddy_data)
    cyc_stop = check_if_in_eddy(date,stop_lon,stop_lat,cyc_eddy_data)
    RCLV_stop = check_if_in_RCLV(date_dashed,stop_lon,stop_lat,cons)
        
    if (anti_start == 0) and (cyc_start == 0) and (RCLV_start == 0): #starts in background waters
        start_flag = 0
        start_str = 'in background waters'
    elif ((anti_start != 0) or (cyc_start != 0)) and (RCLV_start == 0): #starts in SSH eddy
        start_flag = 1
        start_str = 'in SSH eddy'
    elif (anti_start == 0) and (cyc_start == 0) and (RCLV_start == 1): #starts in RCLV
        start_flag = 2
        start_str = 'in RCLV'
    else:
        start_flag = 3
        start_str = 'in overlapping SSH eddy and RCLV'
            
    if (anti_stop == 0) and (cyc_stop == 0) and (RCLV_stop == 0): #starts in background waters
        stop_flag = 0
        stop_str = 'in background waters'
    elif ((anti_stop != 0) or (cyc_stop != 0)) and (RCLV_stop == 0): #starts in SSH eddy
        stop_flag = 1
        stop_str = 'in SSH eddy'
    elif (anti_stop == 0) and (cyc_stop == 0) and (RCLV_stop == 1): #starts in RCLV
        stop_flag = 2
        stop_str = 'in RCLV'
    else:
        stop_flag = 3
        stop_str = 'in overlapping SSH eddy and RCLV'

    intersected_anti_eddy_ID,intersected_anti_eddy_pts = traj_eddy_intersection_pts(date,start_lon,start_lat,stop_lon,stop_lat,anti_eddy_data)
    intersected_cyc_eddy_ID,intersected_cyc_eddy_pts = traj_eddy_intersection_pts(date,start_lon,start_lat,stop_lon,stop_lat,cyc_eddy_data)
    intersected_RCLV_pts = traj_RCLV_intersection_pts(date_dashed,start_lon,start_lat,stop_lon,stop_lat,cons)
        
    eddy_flag = 0
    intersected_lats,intersected_lons = [],[]
        
    if intersected_anti_eddy_ID: # interacted with an anticyclone
        for i in np.arange(0,len(intersected_anti_eddy_ID)):
            anti_IDs.append(intersected_anti_eddy_ID[i])
            for j in np.arange(0,len(intersected_anti_eddy_pts[i])):
                intersected_lons.append(intersected_anti_eddy_pts[i][j][0])
                intersected_lats.append(intersected_anti_eddy_pts[i][j][1])
            eddy_flag += 1
            
    if intersected_cyc_eddy_ID: # interacted with a cyclone
        for i in np.arange(0,len(intersected_cyc_eddy_ID)):
            cyc_IDs.append(intersected_cyc_eddy_ID[i])
            for j in np.arange(0,len(intersected_cyc_eddy_pts[i])):
                intersected_lons.append(intersected_cyc_eddy_pts[i][j][0])
                intersected_lats.append(intersected_cyc_eddy_pts[i][j][1])
            eddy_flag += 1
                
    if intersected_RCLV_pts: # interacted with an RCLV
        for i in np.arange(0,len(intersected_RCLV_pts)):
            for j in np.arange(0,len(intersected_RCLV_pts[i])):
                intersected_lons.append(intersected_RCLV_pts[i][j][0])
                intersected_lats.append(intersected_RCLV_pts[i][j][1])
            eddy_flag += 1
                
    if (eddy_flag == 0): # did not pass through an eddy
        dist = distance_from_lat_lon(start_lat,start_lon,stop_lat,stop_lon)
        print('%s km in background waters from (%s,%s) to (%s,%s)'%(round(dist,2),round(start_lat,3),round(start_lon-360,3),round(stop_lat,3),round(stop_lon-360,3)))
            
    else:
        # Sort intersections by latitude (trajectory is only headed south) so that the events are in the correct order through time
        intersected_lons = [x for _,x in sorted(zip(intersected_lats,intersected_lons))][::-1]
        intersected_lats = sorted(intersected_lats)[::-1]

        if intersected_lats[0] != start_lat: # need to add starting point
            intersected_lons = [start_lon] + intersected_lons
            intersected_lats = [start_lat] + intersected_lats
        elif intersected_lats[1] == start_lat: # start double counted b/c started in overlapping SSH eddy & RCLV
            intersected_lons = intersected_lons[1:]
            intersected_lats = intersected_lats[1:]

        if intersected_lats[-1] != stop_lat: # need to add end point
            intersected_lons = intersected_lons + [stop_lon]
            intersected_lats = intersected_lats + [stop_lat]
        elif intersected_lats[-2] == stop_lat: # stop double counted b/c stopped in overlapping SSH eddy & RCLV
            intersected_lons = intersected_lons[:-1]
            intersected_lats = intersected_lats[:-1]                

        in_eddy = []
        for p in np.arange(0,len(intersected_lons)-1):
            midpt_lon = intersected_lons[p] - (intersected_lons[p] - intersected_lons[p+1])/2
            midpt_lat = intersected_lats[p] - (intersected_lats[p] - intersected_lats[p+1])/2            

            anti_check = check_if_in_eddy(date,midpt_lon,midpt_lat,anti_eddy_data)
            cyc_check = check_if_in_eddy(date,midpt_lon,midpt_lat,cyc_eddy_data)
            RCLV_check = check_if_in_RCLV(date_dashed,midpt_lon,midpt_lat,cons)

            if (anti_check == 0) and (cyc_check == 0) and (RCLV_check == 0): #starts in background waters
                in_eddy.append(0)
            elif ((anti_check != 0) or (cyc_check != 0)) and (RCLV_check == 0): #starts in SSH eddy
                in_eddy.append(1)
            elif (anti_check == 0) and (cyc_check == 0) and (RCLV_check == 1): #starts in RCLV
                in_eddy.append(2)
            else:
                in_eddy.append(3)

        for k in np.arange(0,len(in_eddy)): 

            lat1,lon1 = intersected_lats[k],intersected_lons[k]
            lat2,lon2 = intersected_lats[k+1],intersected_lons[k+1]
            dist = distance_from_lat_lon(lat1,lon1,lat2,lon2)

            if (in_eddy[k] == 0):
                condition_str = 'in background waters'
            elif (in_eddy[k] == 1):
                condition_str = 'in SSH eddy'
            elif (in_eddy[k] == 2):
                condition_str = 'in RCLV'
            elif (in_eddy[k] == 3):
                condition_str = 'in overlapping SSH eddy and RCLV'

            if k == 0:
                print('%s km %s from (%s,%s) to (%s,%s)'%(round(dist,2),condition_str,round(lat1,3),round(lon1-360,3),round(lat2,3),round(lon2-360,3)))
            else:
                print('%s km %s from (%s,%s) to (%s,%s)'%(round(dist,2),condition_str,round(lat1,3),round(lon1-360,3),round(lat2,3),round(lon2-360,3)))

