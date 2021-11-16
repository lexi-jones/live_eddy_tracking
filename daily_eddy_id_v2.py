# Download CMEMS data of a requested date directly to local computer. 
# Lexi Jones 
# Date created: 10/12/21
# Last edited: 11/14/21

################ Import Required Packages #####################

import numpy as np
import math,os,time,sys,subprocess
from datetime import timedelta,datetime
import xarray as xr
from glob import glob
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, DiffusionUniformKh,ErrorCode,Field
from functions_for_parcels import CalcVort, DeleteParticle,particle_grid2d,simulate_particles2d,simulate_particles2d_backwards
from floater import rclv
from skimage.feature import peak_local_max

# User inputs the date of the data they want to download / run particle sim from
date_input = str(input("Enter the date in form yyyy-mm-dd: "))
date_input_no_dash = date_input.replace('-','')
start_year = int(date_input_no_dash[0:4])
start_month = int(date_input_no_dash[4:6])
start_day = int(date_input_no_dash[6:8])

##################### Download CMEMS data #####################


def download_CMEMS_data():

    # Get the current date & time
    now = datetime.now()
    download_datetime = now.strftime("%Y-%m-%d_%H:%M")

    # First check if output dir exists & create one if it does not
    CMEMS_data_dir = os.getcwd() + '/CMEMS_data/'
    CMEMS_data_exists = os.path.exists(CMEMS_data_dir)
    if not CMEMS_data_exists:
        os.makedirs(CMEMS_data_dir)

    # NOTE: If the download command is not working, try the following steps:
    #       1) Go to dataset URL: https://resources.marine.copernicus.eu/product-detail/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/INFORMATION
    #       2) Click ‘Data access’ on right (& log in if prompted)
    #       3) Click ‘Download options’
    #       4) Select ‘View script’ dropdown
    #       5) Check that the command is aligned with the ‘template command’ (product names sometimes change)

    CMEMS_user = str(input("CMEMS Username: "))
    CMEMS_pass = str(input("CMEMS Password: "))

    subprocess.run(["python","-m","motuclient","--motu","https://nrt.cmems-du.eu/motu-web/Motu","--service-id","SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046-TDS", \
    "--product-id","dataset-duacs-nrt-global-merged-allsat-phy-l4","--longitude-min","-165","--longitude-max","-115","--latitude-min","-10","--latitude-max","35", \
    "--date-min",date_input+" 00:00:00","--date-max",date_input+" 00:00:00","--out-dir","CMEMS_data","--out-name", \
    "CMEMS_%s_lat_-10_35_lon_-165_-115_downloaded_%s.nc"%(date_input,download_datetime),"--variable","sla","--variable","ugos","--variable","vgos",\
    "--user",str(CMEMS_user),"--pwd",str(CMEMS_pass)])


#################### Run Lagrangian Simulation ####################


def run_parcels(lat_start,lat_stop,lon_start,lon_stop,spatial_step,time_step,output_freq,runtime,runtime_unit,backwards):

    # Run an OceanParcels simulation with functions from 'functions_for_parcels'. Code comes from run_parcels_CMEMS_v3.py on Engaging.
    # Format: python run_parcels_CMEMS_v3.py lat_start lat_stop lon_start lon_stop spatial_step date_input time_step output_freq runtime runtime_unit backwards
    #       lat_start,lat_stop,lon_start,lon_stop: latitude and longitude bounds to intialize particles
    #       spatial_step: lat/lon spacing between initial particles
    #       date_input: start date with format YYYYMMDD
    #       time_step: (minutes) the amount of time that passes in the fieldset before calculating new particle positions & attributes
    #       output_freq: (hours) frequency to output information about the particle positions & attributes
    #       runtime: how long to run the simulation
    #       runtime_unit: 'days' or 'hours', corresponds with runtime
    #       backwards: 'y' if the particles should be run backwards in time, 'n' if the particles should be run forwards in time

    class SpinnyParticle(JITParticle):
        u = Variable('u',dtype=np.float64)
        v = Variable('v',dtype=np.float64)
        vort = Variable('vort',dtype=np.float64)

    # Prep the list of file names by getting the most recently downloaded CMEMS files per day. File names should be in the
    # format from the download_CMEMS_data() function.

    input_dir = os.getcwd() + '/CMEMS_data/'
    all_input_files = sorted(glob(input_dir+'*.nc'))
    latest_downloads = []
    date = ''
    for file in all_input_files[::-1]:
        if (file.split('/'))[-1][:16] != date:
            latest_downloads.append(file)
            date = (file.split('/'))[-1][:16]
    parcels_input_files = latest_downloads[::-1]

    start_date = datetime(start_year,start_month,start_day)

    # Create fieldset 
    filenames = {'U': parcels_input_files,'V': parcels_input_files}
    variables = {'U': 'ugos','V': 'vgos'}
    dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
                  'V': {'lon':'longitude','lat':'latitude','time':'time'}}
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
    print('Fieldset created.')

    pset_dynamic,num_particles = particle_grid2d(fieldset,SpinnyParticle,[lat_start,lat_stop,spatial_step],[lon_start,lon_stop,spatial_step],start_date)
    print('Particles initialized.')

    parcels_trajs_dir = os.getcwd() + '/parcels_trajs/'
    parcels_trajs_exists = os.path.exists(parcels_trajs_dir)
    if not parcels_trajs_exists:
        os.makedirs(parcels_trajs_dir)

    # Execute particle simulation
    if backwards == 'y':
        print('Particles are being advected backward in time...')
        output_file_path = '%s%s_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)        
        simulate_particles2d_backwards(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq)
    else:
        print('Particles are being advected forward in time...')
        output_file_path = '%s%s_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
        simulate_particles2d(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq)

    print('Output file: %s'%(output_file_path))


########################### Calculate LAVD ##########################


def calc_LAVD(lat_start,lat_stop,lon_start,lon_stop,spatial_step,time_step,output_freq,runtime,runtime_unit,backwards):

    # Open trajectory file
    parcels_trajs_dir = os.getcwd() + '/parcels_trajs/'
    if backwards == 'y':
        parcels_trajs_path = '%s%s_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        parcels_trajs_path = '%s%s_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    
    # Read in data with xarray
    traj_ds = xr.open_dataset(parcels_trajs_path)
    vort_premask = traj_ds.variables["vort"]
    vort = np.array(vort_premask.where(vort_premask != 0)) #filters out land values
    vort_avg_t = np.nanmean(vort,axis=0)[1:] #Find average vorticity over the entire spatial domain at each time step
    LAVD = np.trapz(np.absolute(vort[:,1:] - vort_avg_t), dx=output_freq*60*60, axis=1)/(runtime*24*60*60-output_freq*60*60) #trapz does the integration; convert data to be in seconds units

    LAVD_dir = os.getcwd() + '/LAVD/'
    LAVD_path_exists = os.path.exists(LAVD_dir)
    if not LAVD_path_exists:
        os.makedirs(LAVD_dir)

    if backwards == 'y':
        LAVD_file_path = '%s%s_LAVD_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:    
        LAVD_file_path = '%s%s_LAVD_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    np.save(LAVD_file_path,LAVD)


def calc_LAVD_spatial_subset(lon_stop_ind,lat_stop_ind,lat_start,lat_stop,lon_start,lon_stop,spatial_step,time_step,output_freq,runtime,runtime_unit,backwards):

    # Open trajectory file
    parcels_trajs_dir = os.getcwd() + '/parcels_trajs/'
    if backwards == 'y':
        parcels_trajs_path = '%s%s_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        parcels_trajs_path = '%s%s_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)

    # Read in data with xarray
    traj_ds = xr.open_dataset(output_file_path)
    vort_premask = traj_ds.variables["vort"]
    vort = np.array(vort_premask.where(vort_premask != 0)) #filters out land values
    vort_reshape = vort.reshape((len(np.arange(lon_start,lon_stop,spatial_step)),len(np.arange(lat_start,lat_stop,spatial_step)), vort.shape[1]))
    vort_transpose = np.empty((vort_reshape.shape[2],vort_reshape.shape[1],vort_reshape.shape[0]))

    for i in np.arange(0,vort_reshape.shape[2]):
        vort_transpose[i,:,:] = np.transpose(vort_reshape[:,:,i])

    vort_subset = np.array(vort_transpose[:,:lat_stop_ind,:lon_stop_ind])
    vort_flat = np.empty((vort_subset.shape[0],vort_subset.shape[1]*vort_subset.shape[2]))

    for i in np.arange(0,vort_subset.shape[0]):
        vort_flat[i] = vort_subset[i,:,:].flatten()

    vort_flipped = np.swapaxes(vort_flat,0,1)
    vort_avg_t = np.nanmean(np.array(vort_flipped),axis=0)[1:]
    LAVD = np.trapz(np.absolute(vort_flipped[:,1:] - vort_avg_t), dx=output_freq*60*60, axis=1)/(runtime*24*60*60-output_freq*60*60)


    ######## NOTE: Need to adjust title to reflect the subsetted data
    
    LAVD_dir = os.getcwd() + '/LAVD/'
    LAVD_path_exists = os.path.exists(LAVD_dir)
    if not LAVD_path_exists:
        os.makedirs(LAVD_dir)

    if backwards == 'y':
        LAVD_file_path = '%s%s_LAVD_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        LAVD_file_path = '%s%s_LAVD_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    np.save(LAVD_file_path,LAVD)

    ########

def calc_LAVD_temporal_subset(new_runtime,lat_start,lat_stop,lon_start,lon_stop,spatial_step,time_step,output_freq,runtime,runtime_unit,backwards):
    
    #  Open trajectory file
    parcels_trajs_dir = os.getcwd() + '/parcels_trajs/'
    if backwards == 'y':
        parcels_trajs_path = '%s%s_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        parcels_trajs_path = '%s%s_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc' \
            %(parcels_trajs_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)

    # Read in data with xarray
    traj_ds = xr.open_dataset(output_file_path)
    vort_premask = traj_ds.variables["vort"]
    vort = np.array(vort_premask.where(vort_premask != 0)) #filters out land values

    vort_avg_t = np.nanmean(np.array(vort),axis=0)[1:16] #Find average vorticity over the entire spatial domain at each time step
    LAVD = np.trapz(np.absolute(vort[:,1:new_runtime+1] - vort_avg_t), dx=output_freq*60*60, axis=1)/(new_runtime*24*60*60-output_freq*60*60)


    ######## NOTE: Need to adjust title to reflect the subsetted data

    LAVD_dir = os.getcwd() + '/LAVD/'
    LAVD_path_exists = os.path.exists(LAVD_dir)
    if not LAVD_path_exists:
        os.makedirs(LAVD_dir)

    if backwards == 'y':
        LAVD_file_path = '%s%s_LAVD_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        LAVD_file_path = '%s%s_LAVD_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    np.save(LAVD_file_path,LAVD)

    ########


############################ Find RCLVs #############################


def find_RCLVs(min_dist,min_area,def_tol,target_cd,lat_start,lat_stop,lon_start,lon_stop,spatial_step,time_step,output_freq,runtime,runtime_unit,backwards):
    
    ### Load in LAVD file and reshape it to run through floater
    LAVD_dir = os.getcwd() + '/LAVD/'
    if backwards == 'y':
        LAVD_file_path = '%s%s_LAVD_%s%s_backward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    else:
        LAVD_file_path = '%s%s_LAVD_%s%s_forward_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(LAVD_dir,date_input,runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)

    LAVD = np.load(LAVD_file_path)
    LAVD = np.ma.masked_where(np.isnan(LAVD),LAVD)
    LAVD_reshape = np.transpose(np.reshape(LAVD,(len(np.arange(lon_start,lon_stop,spatial_step)),len(np.arange(lat_start,lat_stop,spatial_step)))))

    # Run floater package
    print("Finding local maxima in the LAVD field...")
    plm = peak_local_max(LAVD_reshape, min_distance=min_dist) #number of pixels between local maxima
    
    print("Finding RCLVs...")
    cons,areas,cds = [],[],[]
    for ji in plm:
        con, area, cd = rclv.convex_contour_around_maximum(LAVD_reshape, ji, convex_def=target_cd, convex_def_tol=def_tol)
        try:
            if area > min_area:
                cons.append(con)
                areas.append(area)
                cds.append(cd)
        except:
            pass
    cons = np.array(cons,dtype=object)
    areas = np.array(areas,dtype=object)
    cds = np.array(cds,dtype=object)

    # Save data
    RCLV_dir = os.getcwd() + '/RCLVs/'
    RCLV_path_exists = os.path.exists(RCLV_dir)
    if not RCLV_path_exists:
        os.makedirs(RCLV_dir)

    if backwards == 'y':
        direction = "backward"
    else: 
        direction = "forward"

    con_file_path = '%s%s_CONS_min_dist_%s_min_area_%s_def_tol_%s_targetCD_%s_%s%s_%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(RCLV_dir,date_input,min_dist,min_area,def_tol,target_cd,runtime,runtime_unit,direction,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    area_file_path = '%s%s_AREAS_min_dist_%s_min_area_%s_def_tol_%s_targetCD_%s_%s%s_%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(RCLV_dir,date_input,min_dist,min_area,def_tol,target_cd,runtime,runtime_unit,direction,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    cd_file_path = '%s%s_CDS_min_dist_%s_min_area_%s_def_tol_%s_targetCD_%s_%s%s_%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(RCLV_dir,date_input,min_dist,min_area,def_tol,target_cd,runtime,runtime_unit,direction,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)
    plm_file_path = '%s%s_PLM_min_dist_%s_min_area_%s_def_tol_%s_targetCD_%s_%s%s_%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.npy' \
            %(RCLV_dir,date_input,min_dist,min_area,def_tol,target_cd,runtime,runtime_unit,direction,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)

    np.save(con_file_path,cons)
    np.save(area_file_path,areas)
    np.save(cd_file_path,cds)
    np.save(plm_file_path,plm)

############################## EXECUTE ##############################

download_CMEMS_data()
run_parcels(10,33,-150,-117,0.03125,20,6,8,'days','y')
calc_LAVD(10,33,-150,-117,0.03125,20,6,8,'days','y')
find_RCLVs(40,24,0.005,0.01,10,33,-150,-117,0.03125,20,6,8,'days','y')
