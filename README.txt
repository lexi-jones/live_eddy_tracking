# Instructions to setup and run 'live_eddy_tracking' software (for Mac/Linux)

# Lexi Jones
# Date Created: 10/13/21
# Last Updated: 11/16/21

###### SET UP CONDA ENV ######
1) If you do not have Anaconda installed, do that from https://www.anaconda.com/products/individual

2) Create a new conda environment with command: 
	conda create --name <ENV NAME>

	To activate this environment:
		conda activate <ENV NAME>
	To deactivate an active environment:
		conda deactivate


###### INSTALL PYTHON PACKAGES ######
# By installing many required packages at once, it helps to avoid package conflicts.
# All of these installations should happen inside the activated conda environment. 

1) Activate new conda environment (see above)

2) Package installations:
	conda install python numpy scipy netcdf4 xarray git motuclient 

	NOTE: If the above command does not work, try running the following one at a time
		conda install python numpy scipy netcdf4 xarray git
		conda install -c conda-forge motuclient
	
3) Install packages for parcels (see https://oceanparcels.org/#installing for more info):
	conda install -c conda-forge parcels jupyter cartopy ffmpeg shapely

4) Download parcels package using git:
	git clone https://github.com/ocean-transport/floater.git

5) Navigate into floater directory and run the command
	python setup.py install

6) Test code by running: 
	python daily_eddy_id_v2.py
	
	NOTE: You will need to download at least 8 days of satellite data before running the full script. Comment out all functions at the
	bottom of the script other than 'download_CMEMS_data()'.
	

###### MATLAB / OceanEddies SETUP #####

1) Install the following add-on requirements in MATLAB:
	- Mapping Toolbox
	- Image Processing Toolbox
	- Statistics and Machine Learning Toolbox

	NOTE: On Mac, you need to agree to the Xcode Terms

2) Install OceanEddies software from Github:
	git clone https://github.com/ifrenger/OceanEddies.git

3) Run the following command inside the '/OceanEddies/mha/' directory (for more details, see https://github.com/jfaghm/OceanEddies):
	python setup.py build_ext -b mht

4) When running for the first time, you must run the 'cell_area.py' script after at least one CMEMS file
is downloaded in the 'CMEMS_data' directory.

5) Alias the MATLAB application so that it can be run from the Terminal with the command
	alias matlab='/Applications/MATLAB_R2021b.app/bin/matlab -nodesktop -nosplash $*'
   where R2019b is replaced with whatever MATLAB version you have on your computer.
	
6) Run MATLAB script from command line (or open MATLAB & run it from the interface)
	matlab -r "G4_CMEMS_SSH_eddies; exit;"

	NOTE: You will either need to download 30 days of satellite data OR change that parameter that is hardcoded in the MATLAB script.
