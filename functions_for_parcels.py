# Lexi Jones
# Last edited: 10/19/21

from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, DiffusionUniformKh,AdvectionRK4_3D,ErrorCode,Field
import numpy as np
import math, os, time
import xarray as xr
import itertools
from glob import glob
from datetime import timedelta

############# FUNCTIONS ##############

def CalcVort(particle, fieldset, time):  
    """
    Calculates vorticity & velocity at the OceanParcels particle location.
    
    Utilizes auxillary coords surrounding the particle (*):
    
        |------x[1]------|
        |                |
        |                |
        |                |
        x[2]     *      x[0]
        |                |
        |                |
        |------x[3]------|
    
    """ 
    
    ## Calculate vorticity
    rho = (0.5)*0.020830 #appx num degrees per cell

    # Find auxillary coords where we will calculate velocities surrounding the particle
    lat1_aux = particle.lat + rho
    lat3_aux = particle.lat - rho
    lon0_aux = particle.lon + rho
    lon2_aux = particle.lon - rho

    # Calculate velocities at auxillary coords
    # U & V velocities are converted to degree/sec under the hood, so we have to convert them back to m/s here
    u1_aux = fieldset.U[time,particle.depth,lat1_aux,particle.lon]*1852*60*math.cos(particle.lat*(math.pi/180))
    u3_aux = fieldset.U[time,particle.depth,lat3_aux,particle.lon]*1852*60*math.cos(particle.lat*(math.pi/180))
    v0_aux = fieldset.V[time,particle.depth,particle.lat,lon0_aux]*1852*60
    v2_aux = fieldset.V[time,particle.depth,particle.lat,lon2_aux]*1852*60

    # Calculate distance in meters between the auxillary coords
    earth_radius = 6.3781*(10**6) #meters
    deg_to_radians = math.pi/180
    lat_dist = earth_radius*(2*rho)*deg_to_radians
    lon_dist = lat_dist*math.cos(particle.lat*deg_to_radians) #depends on the latitude of the particle

    #calculate dv/dx and du/dy
    grad_xv = ((v0_aux-v2_aux)/lon_dist)
    grad_yu = ((u1_aux-u3_aux)/lat_dist)
    
    ## Add the custom data to the output  
    particle.vort = grad_xv - grad_yu
    
    #U & V velocities are converted to degree/sec under the hood, so we have to convert them back to m/s here
    particle.u = fieldset.U[time,particle.depth,particle.lat,particle.lon]*1852*60*math.cos(particle.lat*(math.pi/180))
    particle.v = fieldset.V[time,particle.depth,particle.lat,particle.lon]*1852*60

def DeleteParticle(particle, fieldset, time):
    """
    Stop a particle from running through the simulation. Can be used if a particle runs on land or out of bounds of the grid.
    """
    particle.delete()

def particle_grid2d(fieldset,custom_particle,lat_params,lon_params,time):
    """
    Create lats and lons arrays to create a 2D grid of particles to feed into ParticleSet.from_list() method.
    
    Input
        custom_particle: The class for custom particle outputs 
        _params: list of the form [_start,_stop,_step]        
    Ouput
        pset_grid: Particle set to run through OceanParcels
    """
    lon_array = np.arange(lon_params[0],lon_params[1],lon_params[2])
    lat_array = np.arange(lat_params[0],lat_params[1],lat_params[2])
    
    if (not len(lon_array)):
        print('Warning: Lon_array is empty.')
    if (not len(lat_array)):
        print('Warning: Lat_array is empty.')
    
    #Format desired coords for reading into `from_list` function
    lats = list(itertools.chain.from_iterable([list(lat_array)]*len(lon_array)))
    lons = list(itertools.chain.from_iterable([[l]*len(lat_array) for l in lon_array]))
    
    print('Number of particles: %s'%(len(lats)))
    
    pset_grid = ParticleSet.from_list(fieldset = fieldset, pclass = custom_particle, lon = lons, lat = lats, depth = [0]*len(lons), time = time)
    
    return pset_grid,len(lats)

def simulate_particles2d(pset,output_file_path,runtime,runtime_unit,timestep_mins,output_hrs):
    """
    Simulate particles through the trajectory field.
    
    Input
        runtime_unit: 'hours' or 'days'
        runtime: amount of time (with the given units) to simulate particle trajectories
        timestep_mins: Time between when particle location is calculated
        output_hrs: Time between when particle location is outputed
        num_particles: number of particles that will be simulated
        
    Output
        output_file will write the data if requested, otherwise the pset will now contain particles in their new locations
    """

    ## Create output file & remove contents of the existing file (if it already exists) to ensure the code writes the new data to the file
    if os.path.exists(output_file_path):
        f = open(output_file_path, 'w')
        f.close()      
    output_file = pset.ParticleFile(name=output_file_path,outputdt=timedelta(hours=output_hrs))

    ## Assign the custom kernel behavior to the particle set  
    custom_kernel = pset.Kernel(CalcVort)    
        
    ## Run the parcels simulation
    if (runtime_unit == 'hours'):
        pset.execute(AdvectionRK4 + custom_kernel,
                    runtime=timedelta(hours=runtime),
                    dt=timedelta(minutes=timestep_mins),
                    output_file=output_file,
                    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        
    elif (runtime_unit == 'days'):
        pset.execute(AdvectionRK4 + custom_kernel,
                    runtime=timedelta(days=runtime),
                    dt=timedelta(minutes=timestep_mins),
                    output_file=output_file,
                    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}) 
    else:
        print('Warning: Output should be hours or days.')
        
    #Write the trajectory information to the output file
    output_file.close() #deletes the Parcels output folder with npy files

def simulate_particles2d_backwards(pset,output_file_path,runtime,runtime_unit,timestep_mins,output_hrs):
    """
    Simulate particles through the trajectory field backwards in time.

    Input
	runtime_unit: 'hours' or 'days'
        runtime: amount of time (with the given units) to simulate particle trajectories
        timestep_mins: Time between when particle location is calculated
        output_hrs: Time between when particle location is outputed
        num_particles: number of particles that will be simulated

    Output
	output_file will write the data if requested, otherwise the pset will now contain particles in their new locations
    """

    ## Create output file & remove contents of the existing file (if it already exists) to ensure the code writes the new data to the file
    if os.path.exists(output_file_path):
        f = open(output_file_path, 'w')
        f.close()
    output_file = pset.ParticleFile(name=output_file_path,outputdt=timedelta(hours=output_hrs))

    ## Assign the custom kernel behavior to the particle set
    custom_kernel = pset.Kernel(CalcVort)

    ## Run the parcels simulation
    if (runtime_unit == 'hours'):
        pset.execute(AdvectionRK4 + custom_kernel,
                    runtime=timedelta(hours=runtime),
                    dt=-timedelta(minutes=timestep_mins),
                    output_file=output_file,
                    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

    elif (runtime_unit == 'days'):
        pset.execute(AdvectionRK4 + custom_kernel,
                    runtime=timedelta(days=runtime),
                    dt=-timedelta(minutes=timestep_mins),
                    output_file=output_file,
                    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    else:
        print('Warning: Output should be hours or days.')

    #Write the trajectory information to the output file
    output_file.close() #deletes the Parcels output folder with npy files


