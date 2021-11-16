% Lexi Jones
% Created: 10/12/21
% Last Edited: 10/19/21

% Run OceanEddies with data in seperate netCDF files for each time step.

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% If you get an error message that says a function in not found in 
% the current folder or on the MATLAB path:
%   Go to the 'Current Folder' menu, right click the /OceanEddies/ 
%   directory, select 'Add to Path' > 'Selected Folders and Subfolders'
addpath(genpath(pwd))

%% Settings (paths, switches, definition of subregions...)

% Switches
if_detect = 1;   % = 1 to identify features (closed contours of quantity)
if_track = 1;    % = 1 to track features over time
if_reformat = 1; % = 1 reformat data to ("chelton-like") matrix-style 
if_plot = 0;     % = 1 to create quick plot of results

% Path of satellite data
pathd = [extractBefore(pwd,'OceanEddies'),'CMEMS_data/'];

% Output path 
formatOut = 'yyyymmdd'; 
today = datestr(date,formatOut);
pathout = ['./OceanEddies/OceanEddiesOutput/',today,'/'];

% File containing information on the area of grid boxes; needs to be
% adjusted to your input files (e.g., calculated based on lon/lat data)!
areafile = [pathd,'CMEMS_area_file_lat_-10_35_lon_-165_-115.csv'];

temfrec = 'days';

% Time stamps of input files (must have 8 digits, e.g., YYYYMMDD)
% Use last 30 days of data
datevect = flip(datestr((datetime(date) - days(0:30)),formatOut));

% Counter for input files
timevec = [1:size(datevect,1)]';

% Sample filename
file_pattern = fullfile([pathd,'CMEMS_',datevect(timevec(end),1:4),'-',datevect(timevec(end),5:6),'-',datevect(timevec(end),7:8),'_lat_-10_35_lon_-165_-115_downloaded_*.nc']);
actual_names = dir(file_pattern);

% This loop will end up choosing the latest downloaded file for the given
% date
for k = 1:length(actual_names)
    filen = actual_names(k).name;
end

% Read in longitude/latitude/data
lon = ncread([pathd,filen],'longitude');
lat = ncread([pathd,filen],'latitude');

regstr = {'samplesubregion','all'};
r = 1; % r = 1 exemplary subselects a region; r = 2: all -here global- domain
      
lonsel = [min(lon),min(lon); max(lon),max(lon)];
latsel = [10,min(lat); max(lat),max(lat)];

% Input for detection routine scan_single.m (see more detail in function)
scan_type='v2'; % Detection method 
cycs = {'anticyc','cyclonic'}; % Local minima (cyclones) or maxima (anticyclones)

% Minimum size of detected features
% I think this is number of pixels?
minsz = 15; % Mostly set > 1, e.g,. = 4
minage = 3; % minimum lifespan of eddies to plot/save
   
% Input for tracking routine tolerance_track_lnn.m (see more detail in function)
time_frequency = 1; % Number of days between time steps
tolno = 1; %  Number of TIMESTEPS (here 1 time step is = time_frequency) an eddy can "disappear"

rotsig = {1,-1}; % =1 anticyclones, =-1 cyclones

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
%% Load data
%
% Get lon/lat corners of subregion 
xmin= lonsel(1,r);
xmax= lonsel(2,r);
ymin= latsel(1,r);
ymax= latsel(2,r);

% Find indices of (sub)region
idxx = find(lon>=xmin & lon<=xmax);
idxy = find(lat>=ymin & lat<=ymax);
lona = lon(idxx);
lata = lat(idxy);

% Area matrix (needs to be adjusted to your grid/domain!)
load(areafile)
arn = CMEMS_area_file_lat__10_35_lon__165__115(idxy,idxx)'; 
%
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
%% Identify features

% ----------------------------------------------------------------------
if if_detect ==1
    display('Switch on to detect eddies.')

    % Create subdirectories for output data
    if ~exist(pathout, 'dir')
     mkdir(pathout)
    end
    eddies_save_path = [pathout,'Data/Detection/',regstr{r},'/'];
    system(['rm -rf ',eddies_save_path]);
    if ~exist(eddies_save_path, 'dir')
        mkdir(eddies_save_path);
    end

    clear eddies
    tic % Measures computation time
    for t=1:length(timevec) % loop over each time step

        % Date suffix
        datestrv=datevect(int32(timevec(t)),:);
        display(datestrv(1,:));

        % Filename
        file_pattern = fullfile([pathd,'CMEMS_',datestrv(1:4),'-',datestrv(5:6),'-',datestrv(7:8),'_lat_-10_35_lon_-165_-115_downloaded_*.nc']);
        actual_names = dir(file_pattern);

        % This loop will end up choosing the latest downloaded file for the given
        % date
        for k = 1:length(actual_names)
            filen = actual_names(k).name;
        end
        
        % Longitude should take values -180...180
        lonn = wrapTo180(lona);
        i = find(min(lonn)==lonn);
        lonn = double(circshift(lonn,-(i-1)));
        latn = double(lata);

        % Read in (subregion of) sea level anomaly
        sla = ncread([pathd,filen],'sla',[idxx(1),idxy(1),1],[length(idxx),length(idxy),Inf]);
        % sla data in [cm]
        m2cm = 100; % convert from m to cm
        field = m2cm*double(circshift(sla,[-(i-1) 0])); 

        for rot=1:2 % cyclones and anticyclones

            % cyclones or anticyclones
            cyc = cycs{rot};

            % Feature identification
            eddies = scan_single(field',latn,lonn,timevec(t),cyc,...
                scan_type,arn','minimumArea',minsz,...
                'thresholdStep',.01,'isPadding',0,'sshUnits','centimeters');
            %##### Required arguments
            %1. `field` - input variable, 2D map, with nans for land, size should be `[length(lat) length(lon)]`
            %2. `latn` - 1D array of the latitudes of field grid
            %3. `lonn` - 1D array of the longitudes of field grid
            %4. 'timevec' - used to create output file name
            %4. `cyc` - `'anticyc'` or `'cyclonic'` (eddy polarity, will be attributed based
            %            on local minimum/or maximum of field
            %5. `scan_type` -  `'v1'`, `'v2'`, `'hybrid'`
            %  - here v2: Will find outermost closed contour enclosing a single extremum (labelled "bottom-up)
            %6. `areamap` - A 2D array that refer to the area of each pixel in SSH data (should have same size as ssh), or 1D array that refer to area of each pixel for a specific lat in a regular grid (pixel have same area for the same latitude)

            %##### Optional parameters (only applicable for v2 eddyscan):
            %1. `'minimumArea'` - minimum number of pixels for an eddy, used for validating eddies, default value is `9`
            %2. `'thresholdStep'` - the minimum step for thresholding, the unit is SSH's unit, default value is `0.05`
            %3. `'isPadding'` - whether or not to pad SSH data, should be true when scanning SSH data with the longitudes expanding the whole world dmap. Set to false if only partial SSH data is used. Default value is true
            %4. `'sshUnits'` -  The units the field data; bottom_up_single is built to work natively on centimeter SSH data.  Valid parameters are `'meters'` and `'centimeters'`. If the parameter passed in is `'meters'`, the SSH data will be multiplied by 100. No changes will be made if the parameter passed in is `'centimeters'`.  The default value of 'sshUnits' is centimeters.

            save([eddies_save_path,cyc,'_',num2str(round(timevec(t)),'%08.0f')],'-v7.3','eddies');

            clear sla eddies
        end

    end
    u=toc; display(['Time elapsed in min: ',num2str(u/60)])

end % switch for detection

clear rot t m2cm u i ans cyc s field latn lonn
%
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
%% Track features

% Select time period features shall be tracked over
t1=timevec(1);
t2=timevec(end);

% Directory names for output
eddies_save_path = [pathout,'Data/Detection/',regstr{r},'/'];
tracks_save_path = [pathout,'Data/Tracking/',regstr{r},'/'];
viewer_data_save_path = [pathout,'Viewer/',regstr{r},'/'];
modified_eddies_path = [viewer_data_save_path, 'modified_eddies/'];
modified_tracks_path = [viewer_data_save_path, 'modified_tracks/'];

% ----------------------------------------------------------------------
if if_track == 1
    display('Switch on to track eddies.')

    % Clean up/create subdirectories
    system(['rm -rf ',tracks_save_path]);
    system(['rm -rf ',viewer_data_save_path]);
    if ~exist(tracks_save_path, 'dir')
        mkdir(tracks_save_path);
    end
    if ~exist(viewer_data_save_path, 'dir')
        mkdir(viewer_data_save_path);
    end
    if exist(modified_eddies_path, 'dir')
          system(['rm -rf ',modified_eddies_path]);
    end
    mkdir(modified_eddies_path);
    if exist(modified_tracks_path, 'dir')
        system(['rm -rf ',modified_tracks_path]);
    end
    mkdir(modified_tracks_path);

    clear eddies_track
    tic
    for rot=1:2

        % Feature tracking
        [tmp tmprev tmpdrop tmpdropv]= tolerance_track_lnn(eddies_save_path,cycs{rot},time_frequency,tolno,minsz);

        eval(['!rm -f ',modified_eddies_path,'*',cycs{rot},'*'])
        if ~exist(modified_eddies_path, 'dir')
            mkdir(modified_eddies_path);
        end
        eval(['!rm -rf ',modified_tracks_path,'*',cycs{rot},'*'])
        if ~exist(modified_tracks_path, 'dir')
            mkdir(modified_tracks_path);
        end

        disp('Modifying eddy data to include fake eddies (for the sake of the tracks) - takes time');
        process_eddies_and_tracks_tolerance(cycs{rot}, timevec(t1:t2), eddies_save_path, tmp, modified_eddies_path, modified_tracks_path);
        vars = load([modified_tracks_path, cycs{rot},'_tracks_processed.mat']);
        s = fieldnames(vars);
        tmp = vars.(s{1});

        eddies_track{rot,:} = tmp;
    end
    u=toc; display(['Time elapsed in min: ',num2str(u/60)])
    save([pathout,'Data/eddy-detection_eddies_tracked_',num2str(t1),'-',num2str(t2),'_timenostr_',regstr{r}],'-v7.3','eddies_track','tmprev','tmpdrop','tmpdropv');

else % tracksw=0
    load([pathout,'Data/eddy-detection_eddies_tracked_',num2str(t1),'-',num2str(t2),'_timenostr_',regstr{r}],'eddies_track');
end

clear u ans datestrv rot tmp tmprev tmpdrop tmpdropv s vars
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
%% Reformat "chelton like), i.e. eddy positions + properties in one matrix

eddy_dir = modified_eddies_path;

if if_reformat == 1
    
    reformatted_data = reformat_track_data_to_chelton(eddy_dir, timevec(t1:t2), eddies_track{1,:}, eddies_track{2,:});
    save([pathout,'Data/eddy-detection_eddies_tracked_reform_',num2str(t1),'-',num2str(t2),'timenostr_',regstr{r}],'-v7.3','reformatted_data');

else % tracksw=0
    load([pathout,'Data/eddy-detection_eddies_tracked_reform_',num2str(t1),'-',num2str(t2),'timenostr_',regstr{r}],'reformatted_data');
end

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
%% FIGURE: Plot with tracks, only eddies of certain lifespans
% Background shows input field for a given time step, outer boundaries of 
% tracked features in the same time step are shown in red (positive anomalies, i.e., anticyclones) 
% and blue (negative anomalies, i.e. cyclones); dots mark the centers (extrema) 
% of tracked features, lines the pathways of features taken until present time step. 
%

colorstr = {[139/255 0 0],[0 0 139/255]}; % colors for anticyclones (red)/cyclones (blue)

pathplot = [pathout,'Plots/',regstr{r}];
if ~exist([pathplot], 'dir')
 mkdir([pathplot])
end

figure(1)
clf;hold on
set(gcf,'Position',[45   224   609   481],'Color','w')

boundsw = 1; % Plot eddy boundaries
ddt = 1; 
dates = timevec(t1:ddt:t2); % What time steps to plot

% Load tracked eddies
load([pathout,'Data/eddy-detection_eddies_tracked_reform_',num2str(t1),'-',num2str(t2),'timenostr_',regstr{r}],'reformatted_data');
for t=1:length(dates) % Loop over time steps
    clf

    % Date string
    datestrv=datevect(int32(timevec(t)),:);

    % File name
    file_pattern = fullfile([pathd,'CMEMS_',datestrv(1:4),'-',datestrv(5:6),'-',datestrv(7:8),'_lat_-10_35_lon_-165_-115_downloaded_*.nc']);
    actual_names = dir(file_pattern);
    for k = 1:length(actual_names)
        filen = actual_names(k).name;
    end
    
    tf = knnsearch(timevec,t);

    % Read in field (here sla)
    var = ncread([pathd,filen],'sla',[idxx(1),idxy(1),1],[length(idxx),length(idxy),Inf]);

    % Background: sla in colors
    contourf(lona,lata,var','levelstep',.01,'linestyle','none');hold on
    title(['Sea level height anomaly [m] on ',datevect(t,:)],'fontsize',18)
    caxis([-0.3 0.3])
    colormap((colormap_redblue(24)));colorbar('linewidth',1,'fontsize',18);

    % Add approx land contour
    topoprox=var;
    topoprox(~isnan(topoprox)) = 1;
    topoprox(isnan(topoprox)) = 0;
    contour(wrapTo180(lona),lata,topoprox',[0:.1:.1],'k-','linewidth',2)

    for rot=1:2;

        % determine indices of eddies of certain rotation
        [trl didxun]=hist(reformatted_data.id(...
            reformatted_data.cyc==rotsig{rot}),...
            [1:max(reformatted_data.id)]);
        % trl: number of occurrence of certain eddy id, i.e. lifespan of eddy
        % didxun: id of this eddy

        % Indices with eddies of minimum lifespan
        idx=find(trl*time_frequency > minage);

        % Load detected eddy data file for the respective step (contains
        % eddy properties)
        filen = ls([eddy_dir,cycs{rot},'*',num2str(round(timevec(t)),'%08.0f'),'*']);
        filen = strtrim(filen);
        load(filen,'eddies');
        eddiesvs=eddies;

        % Loop over tracked eddies 
        for e = 1:length(idx)

            % Find time steps the eddy existed
            datv = reformatted_data.track_day(reformatted_data.id==didxun(idx(e)));

            % Check if eddy existed in present time step, if not, skip
            tidx=find(datv==t,1);
            if isempty(tidx);continue;end

            % lon/lat vector of eddy
            lone = reformatted_data.x(reformatted_data.id==didxun(idx(e)));
            late = reformatted_data.y(reformatted_data.id==didxun(idx(e)));

            % Get eddy id (search only in present time step)
            ide=reformatted_data.ide(...
                reformatted_data.track_day==t&reformatted_data.id==didxun(idx(e)));

            if boundsw==1
                
                % prepare boundary
                dotidx=eddiesvs(ide).Stats.PixelIdxList;
                amask=zeros(size(var'));
                amask(dotidx)=1;
                apoly=mask2poly(amask);
                [xx yy]=find(amask==1);
                bdry = bwtraceboundary(amask,[xx(1),yy(1)],'N');
                xx=bdry(:,2);
                yy=bdry(:,1);
                
                % Plot boundary
                h6=plot3((lona(xx)),lata(yy),ones(1,length(bdry(:,1)))*10^20,'w-','linewidth',1,'Color','k');
            
            end

            % Plot track
            h=plot(wrapTo180(lone),late,'k-','Linewidth',3,'color','w');
            h=plot(wrapTo180(lone),late,'k-','Linewidth',2,'color',colorstr{rot});

            % Plot present location of eddy center
            h2=plot3((lone(tidx)),late(tidx),10^20,'wo','Markersize',5,'Markerface',colorstr{rot},'Linewidth',1);
 
            % Adjust plot aesthetics
            set(gca,'fontsize',18,'linewidth',1);box on
            xlim([min(wrapTo180(lona)) max(wrapTo180(lona))])
            ylim([min(lata) max(lata)])
            c=colorbar;
            set(c,'fontsize',18,'linewidth',2)
            
        end % Eddy indices

    end

    % Save figure
    set(gcf,'visible','on')
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'InvertHardcopy','off');
    print(gcf,'-dpng',[pathplot,'/eddy-detection_with_tracks_',num2str(round(timevec(t)),'%08.0f'),'_',regstr{r}]);

end


%% Convert output to readable format

% CSV File (one for cycs & one for antis): Date | Eddy ID | Age | Area | Ls | Center Lon | Center Lat | Amplitude | Translation Speed | Boundary Coords
% Age: Number of days eddy has been in dataset
% Lifetime: Max age in the dataset (always the same for the same ID)
% Area: square kilometers
% Ls: radius (km)
% Center Lat/Lon: location of eddy extremum


% Shift longitude
lat = lata;
lon = wrapTo180(lona);
i = find(min(lon)==lon);
lon = circshift(lon,-(i-1));

clear i 

%Output save path
lj_format_save_path = [pathout,'Data/LJ_format/'];
if ~exist(lj_format_save_path, 'dir')
    mkdir(lj_format_save_path);
end

% Load tracked eddies
viewer_data_save_path = [pathout,'Viewer/',regstr{r},'/'];
modified_eddies_path = [viewer_data_save_path, 'modified_eddies/'];
eddy_dir = modified_eddies_path;
load([pathout,'/Data/eddy-detection_eddies_tracked_reform_',num2str(timevec(1)),'-',num2str(timevec(end)),'timenostr_samplesubregion.mat'],'reformatted_data');

% Dates
output_inds = (timevec(1):timevec(end));
datestrv = datevect(int32(timevec(output_inds)),:);

% Load eddy ids
trl = {};
for rot=1:2 
    didxun = (1:max(reformatted_data.id));
    d = diff(didxun)/2;
    edges = [didxun(1)-d(1), didxun(1:end-1)+d, didxun(end)+d(end)];
    edges(2:end) = edges(2:end)+eps(edges(2:end));
    [counts, edges] = histcounts(reformatted_data.id(reformatted_data.cyc==rotsig{rot}),edges);
    trl{rot}=counts;
end

%Loop through cyclonic and anticyclonic eddies
disp('Reformatting eddies to save in CSV...')
for rot=1:2     
    idx=find(trl{rot}*time_frequency > minage); % Indices with eddies of minimum lifespan

    % Loop over tracked eddies in parallel
    parfor e = 1:length(idx)

        if isfile([lj_format_save_path,cycs{rot},'_',num2str(idx(e),'%010.0f')])
            continue
        end
        % Date | Eddy ID | Age | Area | Ls | Center Lon | Center Lat | Amplitude | Translation Speed | Boundary Coords
        file1_mat = [{'Date','Eddy ID','Age (days)','Area (km^2)','Ls (km)','Center Lon','Center Lat','Amplitude (cm)','Translation Speed (m/s)','Boundary Coords'}];

        % Find time steps the eddy existed
        datv = reformatted_data.track_day(reformatted_data.id==didxun(idx(e)));

        for t = datv(1):datv(end)
            % Load detected eddy data file for the respective step (contains eddy properties)
            filen = ls([eddy_dir,cycs{rot},'*',num2str(round(timevec(t)),'%08.0f'),'*']);
            filen = strtrim(filen);
            eddiesload = load(filen,'eddies');
            eddiesvs = eddiesload.eddies;

            % Get eddy age at current time step
            age = t + 1 - datv(1);

            % Get eddy attributes
            id = didxun(idx(e));
            area = reformatted_data.area(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e))); 
            Ls = reformatted_data.Ls(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e))); 
            center_lon = reformatted_data.x(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e))); 
            center_lat = reformatted_data.y(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e))); 
            amplitude = reformatted_data.amp(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e)));
            translation_speed =reformatted_data.u(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e)));

            % Get boundary coordinates
            ide = reformatted_data.ide(reformatted_data.track_day==t & reformatted_data.id==didxun(idx(e))); % Eddy ID only in present time
            dotidx = eddiesvs(ide).Stats.PixelIdxList;
            amask = zeros(size(arn')); %CHECK IF THIS SHOULD BE TRANSPOSED
            amask(dotidx) = 1;
            [xx, yy] = find(amask==1);
            bdry = bwtraceboundary(amask,[xx(1),yy(1)],'N');
            xx = bdry(:,2);
            yy = bdry(:,1);

            % Converting negative longitudes to positive
            bdry_coords = [];
            temp_lon = lon(xx);
            temp_lat = lat(yy);
            for c = 1:length(xx)
                x = temp_lon(c);
                if x < 0
                    x = x + 360; 
                end
                y = temp_lat(c);
                bdry_coords = [bdry_coords,x,y];
            end      

            % boundary coords are turned into a cell here to keep the rows
            % the same length in matlab, but they're written out at the end 
            % of the row when saved as a CSV
            % Date | Eddy ID | Age | Area | Ls | Center Lon | Center Lat | Amplitude | Translation Speed | Boundary Coords
            file1_mat = [file1_mat; str2double(datestrv(t,:)), id, age, area, Ls, center_lon, center_lat, amplitude, translation_speed, {bdry_coords}]; 
        end

        parsave([lj_format_save_path,cycs{rot},'_',num2str(idx(e),'%010.0f')],file1_mat);
    end
end

clear i temp_array freq rot unique_ids

%% Convert eddy matrix files to cyc and anticyc CSVs

lj_format_save_path = [pathout,'Data/LJ_format/'];

disp('Saving data & finishing up...')
for rot=1:2     
    eddy_files = dir(fullfile(lj_format_save_path,[cycs{rot},'*.mat'])); %gets all wav files in struct
    file2_mat = [{'Date','Eddy ID','Age (days)','Area (km^2)','Ls (km)','Center Lon','Center Lat','Amplitude (cm)','Translation Speed (m/s)','Boundary Coords'}];

    parfor k = 1:length(eddy_files)    
        single_eddy_data=load([lj_format_save_path,eddy_files(k).name]).x;
        file2_mat = [file2_mat;single_eddy_data(2:end,1:end)];
    end
    
    writecell(file2_mat,strcat(pathout,'/Data/',cycs{rot},'_eddy_data_minage',int2str(minage),'_',string(datevect(1,:)),'_to_',string(datevect(end,:)),'.csv'));
end
%%

clear e rot t c xx yy ide tidx datv idx trl didxun var lone late dotidx amask
clear apoly bdry bdry_coords id area amplitude translation_speed age output_inds
clear temp_lon temp_lat 

%%

function parsave(fname, x)
  save(fname, 'x')
end