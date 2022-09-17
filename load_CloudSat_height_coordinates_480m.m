% Script that loads CloudSat 2B-GEOPROF-LIDAR, 2B-FLXHR-LIDAR, and ECMWF-AUX data
% and calculates zonal-mean pendad-mean values. 
% The script assumes that the user has saved CloudSat data
% products 2B-FLXHR-LIDAR, 2B-GEOPROF-LIDAR, and ECMWF-AUX in separate
% directories, and for each data product, all data files are saved in the same directory. 
% The user will need to change the directory sturcture in the script below to match their own machine.
% The script then proceeds by loading data files from one pentad at a time,
% averaging data accordingly, and then saving the values in a grid. Data files are only
% loaded if data are available for all three products (2B-FLXHR-LIDAR, 
% 2B-GEOPROF-LIDAR, and ECMWF-AUX).
% Note: this script needs to be run before "LW_CRH_decomposition.m" and
% "compute_LW_CRH_climatology_and_SAM_regression.m"
clearvars; clc;
addpath ./functions/

% parameters for latitude-time binning
lat_bins=(-82.5):2.5:2.5; lat_mdpts=(lat_bins(1:(end-1))+lat_bins(2:end))./2;
time_bins=datenum(2007,1,1):5:datenum(2010,12,31); time_mdpts=(time_bins(1:(end-1))+time_bins(2:end))./2;
z_bins=0:(0.24.*2):(18+(.24.*2).*2); z_mdpts=(z_bins(1:(end-1))+z_bins(2:end))./2;
lat_min=-83; lat_max=3;

% constants
g=9.81; % graviatational constant (m/s^2)
cp=1004; % specific heat capacity of dry air at constant pressure (J/K kg)

% initialize arrays
CloudFraction_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_hr_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_hr_clr_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_ACRE_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_up_clr_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_down_clr_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
Temperature_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
Specific_humidity_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
U_velocity_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
V_velocity_gridded=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
num_obs=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
H_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
M_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
L_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
LM_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
LH_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
MH_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
LMH_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
clear_frac=NaN*ones(numel(lat_mdpts),numel(time_mdpts));

for t=1:numel(time_mdpts)
    % identify files with data in current time range
    current_time_range=(time_bins(t)-1):(time_bins(t)+4);
    
    fnames_CloudSat={}; % file names for 2B-GEOPROF-LIDAR
    fnames_ECMWF={}; % files names for ECMWF-AUX
    fnames_FLXHR={}; % file names for 2B-FLXHR-LIDAR
    for n=1:numel(current_time_range)
        current_year=year(current_time_range(n));
        current_day=current_time_range(n)-datenum(current_year,1,1)+1;
        
        % convert to string
        current_year=num2str(current_year);
        if current_day<10
            current_day=['00' num2str(current_day)];
        elseif current_day>=10 & current_day<100
            current_day=['0' num2str(current_day)];
        else
            current_day=num2str(current_day);
        end
        
        % check for files
        fnamesa=dir(['/data/cawall/CloudSat/2B_GEOPROF_LIDAR/' current_year current_day '*.hdf']); % the user will need to change this path
        for j=1:numel(fnamesa)
            % check for corresponding ECMWF-AUX data
            current_name=fnamesa(j).name;
            granule_id=current_name(1:19);           
            fname_ECMWF=dir(['/data/cawall/CloudSat/ECMWF_AUX/' granule_id '*.hdf']); % the user will need to change this path
            
            % check for corresponding 2B-FLXHR-LIDAR data
            fname_FLXHR=dir(['/data/cawall/CloudSat/2B_FLXHR_LIDAR/' granule_id '*.hdf']); % the user will need to change this path
            
            if ~isempty(fname_ECMWF) & ~isempty(fname_FLXHR)
                fnames_CloudSat{end+1}=fnamesa(j).name;
                fnames_ECMWF{end+1}=fname_ECMWF.name;
                fnames_FLXHR{end+1}=fname_FLXHR.name;
            end
            clearvars current_name granule_number fname_ECMWF fname_FLXHR          
        end
        clearvars fnamesa j current_year current_day                
    end
    clearvars n
    
    %%%%%% load data %%%%%%%
    nfiles=numel(fnames_CloudSat);
    if nfiles==0 % no data for current pentad
        CloudFraction_gridded(:,:,t)=NaN;
        LW_hr_gridded(:,:,t)=NaN;
        LW_hr_clr_gridded(:,:,t)=NaN;
        LW_ACRE_gridded(:,:,t)=NaN;        
        LW_up_clr_gridded(:,:,t)=NaN;
        LW_down_clr_gridded(:,:,t)=NaN;
        Temperature_gridded(:,:,t)=NaN;
        Specific_humidity_gridded(:,:,t)=NaN;
        U_velocity_gridded(:,:,t)=NaN;
        V_velocity_gridded(:,:,t)=NaN;
        num_obs(:,t)=0;
        H_frac(:,t)=NaN;
        M_frac(:,t)=NaN;
        L_frac(:,t)=NaN;
        LM_frac(:,t)=NaN;
        LH_frac(:,t)=NaN;
        MH_frac(:,t)=NaN;
        LMH_frac(:,t)=NaN;
        clear_frac(:,t)=NaN;
    else % data available
        
        for n=1:nfiles
            % load 2B-GEOPROF-LIDAR
            path=['/data/cawall/CloudSat/2B_GEOPROF_LIDAR/' fnames_CloudSat{n}]; % the user will need to change this path            
            
            % load time array and select data from current pendad
            TAI_start = hdfread(path, '/2B-GEOPROF-LIDAR/Geolocation Fields/TAI_start');
            TAI_start = double(TAI_start{1})./86400 + datenum(1993,1,1);
            Profile_time = hdfread(path, '/2B-GEOPROF-LIDAR/Geolocation Fields/Profile_time');
            timea = double(Profile_time{1})./86400 + TAI_start;
            time_ind=find(timea>=time_bins(t) & timea<time_bins(t+1));
            clearvars TAI_start Profile_time timea
            
            if ~isempty(time_ind)
                
                % load GEOPROF-LIDAR
                path=['/data/cawall/CloudSat/2B_GEOPROF_LIDAR/' fnames_CloudSat{n}]; % the user will need to change this path
                Latitude = hdfread(path, '/2B-GEOPROF-LIDAR/Geolocation Fields/Latitude');
                Latitude = Latitude{1};
                Latitude = Latitude(time_ind);
                Height = double(hdfread(path, '/2B-GEOPROF-LIDAR/Geolocation Fields/Height'))./1000;
                Height = Height(time_ind,:);
                CloudFraction = double(hdfread(path, '/2B-GEOPROF-LIDAR/Data Fields/CloudFraction'));
                CloudFraction = CloudFraction(time_ind,:);
                CloudLayers = hdfread(path, '/2B-GEOPROF-LIDAR/Data Fields/CloudLayers', 'Fields', 'CloudLayers');
                CloudLayers=CloudLayers{1};
                CloudLayers=CloudLayers(time_ind);
                LayerTop = single(hdfread(path, '/2B-GEOPROF-LIDAR/Data Fields/LayerTop'));
                LayerTop = LayerTop(time_ind,:);
                LayerBase = single(hdfread(path, '/2B-GEOPROF-LIDAR/Data Fields/LayerBase'));
                LayerBase = LayerBase(time_ind,:);
                
                % load FLXHR-LIDAR
                path=['/data/cawall/CloudSat/2B_FLXHR_LIDAR/' fnames_FLXHR{n}]; % the user will need to change this path
                FD = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FD');
                LW_down = double(squeeze(FD(2,time_ind,:)))./10;
                FU = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FU');
                LW_up = double(squeeze(FU(2,time_ind,:)))./10;
                FD_NC = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FD_NC');
                LW_down_clr = double(squeeze(FD_NC(2,time_ind,:)))./10;
                FU_NC = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FU_NC');
                LW_up_clr = double(squeeze(FU_NC(2,time_ind,:)))./10;
                QR = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/QR');
                LW_hr = double(squeeze(QR(2,time_ind,:)))./100;
                clearvars FD FU FD_NC FU_NC QR                                                
                
                % fill in missing/bad data with NaN
                CloudFraction(CloudFraction<0)=NaN; LayerTop(LayerTop<0)=NaN; LayerBase(LayerBase<0)=NaN;     
                LW_down(LW_down==-999)=NaN; LW_up(LW_up==-999)=NaN; LW_down_clr(LW_down_clr==-999)=NaN; 
                LW_up_clr(LW_up_clr==-999)=NaN; LW_hr(LW_hr==-999./100)=NaN;
                
                % load ECMWF-AUX data
                path=['/data/cawall/CloudSat/ECMWF_AUX/' fnames_ECMWF{n}]; % the user will need to change this path
                Temperature = double(hdfread(path, '/ECMWF-AUX/Data Fields/Temperature'));
                Temperature = Temperature(time_ind,:);
                Specific_humidity = double(hdfread(path, '/ECMWF-AUX/Data Fields/Specific_humidity'));
                Specific_humidity = Specific_humidity(time_ind,:);
                Specific_humidity(Specific_humidity<0)=NaN;
                U_velocity = double(hdfread(path, '/ECMWF-AUX/Data Fields/U_velocity'));
                U_velocity = U_velocity(time_ind,:);
                V_velocity = double(hdfread(path, '/ECMWF-AUX/Data Fields/V_velocity'));
                V_velocity = V_velocity(time_ind,:);
                Pressure = double(hdfread(path, '/ECMWF-AUX/Data Fields/Pressure'));    
                Pressure = Pressure(time_ind,:);
                Pressure(Pressure<0)=NaN;
                Pressure_half_int=NaN*ones(size(Pressure,1),size(Pressure,2)+1);
                Pressure_half_int(:,2:(end-1))=(Pressure(:,1:(end-1))+Pressure(:,2:end))./2;
                Pressure_half_int(:,1)=Pressure(:,1)-abs(Pressure_half_int(:,2)-Pressure(:,1));
                Pressure_half_int(:,end)=Pressure(:,end)+abs(Pressure(:,end)-Pressure_half_int(:,end-1));
                dp=abs(Pressure_half_int(:,1:(end-1))-Pressure_half_int(:,2:end));
                
                % calculate clear-sky LW heating rate
                LW_hr_clr=86400.*(g./cp).*(LW_down_clr(:,1:(end-1)) - LW_up_clr(:,1:(end-1)) - LW_down_clr(:,2:end) + LW_up_clr(:,2:end))./dp;
                clearvars Pressure_half_int path time_ind LW_up LW_down
                
                % average clear-sky upward and downward LW fluxes to
                % midpoints between their native levels for consistency
                % with vertical resolution of other variables
                LW_up_clr=(LW_up_clr(:,1:(end-1))+LW_up_clr(:,2:end))./2;
                LW_down_clr=(LW_down_clr(:,1:(end-1))+LW_down_clr(:,2:end))./2;
                
                % select data from SH
                lat_ind=find(Latitude>=lat_min & Latitude<=lat_max);
                Height=Height(lat_ind,:);
                Latitude=Latitude(lat_ind);
                CloudFraction=CloudFraction(lat_ind,:);                  
                CloudLayers=CloudLayers(lat_ind);
                LayerTop=LayerTop(lat_ind,:);
                LayerBase=LayerBase(lat_ind,:);
                LW_hr=LW_hr(lat_ind,:);
                LW_hr_clr=LW_hr_clr(lat_ind,:);                
                LW_up_clr=LW_up_clr(lat_ind,:);
                LW_down_clr=LW_down_clr(lat_ind,:);
                Temperature=Temperature(lat_ind,:);
                Specific_humidity=Specific_humidity(lat_ind,:);
                U_velocity=U_velocity(lat_ind,:);
                V_velocity=V_velocity(lat_ind,:);
                dp=dp(lat_ind,:);
                clearvars lat_ind                             
                
                % L/M/H cloud categorization
                cloud_labels = cloud_categorization(CloudLayers,LayerTop,LayerBase);
                
                % save data for later use
                if ~exist('Latitude_save')
                    nelements=5e4.*nfiles; % random large number used initially to create array
                    Latitude_save=NaN*ones(1,nelements);
                    Height_save=NaN*ones(nelements,size(Height,2));
                    CloudFraction_save=NaN*ones(nelements,size(Height,2));
                    cloud_labels_save=NaN*ones(nelements,size(cloud_labels,2));
                    dp_save=NaN*ones(nelements,size(Height,2));
                    LW_hr_save=NaN*ones(nelements,size(Height,2));
                    LW_hr_clr_save=NaN*ones(nelements,size(Height,2));                   
                    LW_up_clr_save=NaN*ones(nelements,size(Height,2));
                    LW_down_clr_save=NaN*ones(nelements,size(Height,2));
                    Temperature_save=NaN*ones(nelements,size(Height,2));
                    Specific_humidity_save=NaN*ones(nelements,size(Height,2));
                    U_velocity_save=NaN*ones(nelements,size(Height,2));
                    V_velocity_save=NaN*ones(nelements,size(Height,2));
                    
                    current_index_start=1;
                    current_index_end=numel(Latitude);
                    Latitude_save(1,current_index_start:current_index_end)=Latitude;
                    Height_save(current_index_start:current_index_end,:)=Height;
                    CloudFraction_save(current_index_start:current_index_end,:)=CloudFraction;
                    cloud_labels_save(current_index_start:current_index_end,:)=cloud_labels;
                    dp_save(current_index_start:current_index_end,:)=dp;
                    LW_hr_save(current_index_start:current_index_end,:)=LW_hr;
                    LW_hr_clr_save(current_index_start:current_index_end,:)=LW_hr_clr;                   
                    LW_up_clr_save(current_index_start:current_index_end,:)=LW_up_clr;
                    LW_down_clr_save(current_index_start:current_index_end,:)=LW_down_clr;
                    Temperature_save(current_index_start:current_index_end,:)=Temperature;
                    Specific_humidity_save(current_index_start:current_index_end,:)=Specific_humidity;
                    U_velocity_save(current_index_start:current_index_end,:)=U_velocity;
                    V_velocity_save(current_index_start:current_index_end,:)=V_velocity;
                else
                    current_index_start=current_index_end+1;
                    current_index_end=current_index_start+numel(Latitude)-1;
                    Latitude_save(1,current_index_start:current_index_end)=Latitude;
                    Height_save(current_index_start:current_index_end,:)=Height;
                    CloudFraction_save(current_index_start:current_index_end,:)=CloudFraction;
                    cloud_labels_save(current_index_start:current_index_end,:)=cloud_labels;
                    dp_save(current_index_start:current_index_end,:)=dp;
                    LW_hr_save(current_index_start:current_index_end,:)=LW_hr;
                    LW_hr_clr_save(current_index_start:current_index_end,:)=LW_hr_clr;                    
                    LW_up_clr_save(current_index_start:current_index_end,:)=LW_up_clr;
                    LW_down_clr_save(current_index_start:current_index_end,:)=LW_down_clr;
                    Temperature_save(current_index_start:current_index_end,:)=Temperature;
                    Specific_humidity_save(current_index_start:current_index_end,:)=Specific_humidity;
                    U_velocity_save(current_index_start:current_index_end,:)=U_velocity;
                    V_velocity_save(current_index_start:current_index_end,:)=V_velocity;
                end
                                            
            end
            
            if t==1
                n
            end
        end
              
        if ~exist('current_index_end') % no available data for this pentad
            CloudFraction_gridded(:,:,t)=NaN;
            LW_hr_gridded(:,:,t)=NaN;
            LW_hr_clr_gridded(:,:,t)=NaN;
            LW_ACRE_gridded(:,:,t)=NaN;            
            LW_up_clr_gridded(:,:,t)=NaN;
            LW_down_clr_gridded(:,:,t)=NaN;
            Temperature_gridded(:,:,t)=NaN;
            Specific_humidity_gridded(:,:,t)=NaN;
            U_velocity_gridded(:,:,t)=NaN;
            V_velocity_gridded(:,:,t)=NaN;
            num_obs(:,t)=0;
            H_frac(:,t)=NaN;
            M_frac(:,t)=NaN;
            L_frac(:,t)=NaN;
            LM_frac(:,t)=NaN;
            LH_frac(:,t)=NaN;
            MH_frac(:,t)=NaN;
            LMH_frac(:,t)=NaN;
            clear_frac(:,t)=NaN;          
        else
        
            Latitude_save=Latitude_save(1:current_index_end);
            Height_save=Height_save(1:current_index_end,:);
            CloudFraction_save=CloudFraction_save(1:current_index_end,:);
            cloud_labels_save=cloud_labels_save(1:current_index_end,:);
            dp_save=dp_save(1:current_index_end,:);
            LW_hr_save=LW_hr_save(1:current_index_end,:);
            LW_hr_clr_save=LW_hr_clr_save(1:current_index_end,:);           
            LW_up_clr_save=LW_up_clr_save(1:current_index_end,:);
            LW_down_clr_save=LW_down_clr_save(1:current_index_end,:);
            Temperature_save=Temperature_save(1:current_index_end,:);
            Specific_humidity_save=Specific_humidity_save(1:current_index_end,:);
            U_velocity_save=U_velocity_save(1:current_index_end,:);
            V_velocity_save=V_velocity_save(1:current_index_end,:);
            clearvars CloudFraction Height Latitude current_index_start current_index_end nelements n nfiles ...
                fnames_CloudSat cloud_labels dp LW_hr LW_hr_clr Temperature Specific_humidity LW_up_clr LW_down_clr ...
                U_velocity V_velocity
            
            % bin by latitude and height
            Latitude_save=transpose(Latitude_save);
            Latitude_save=repmat(Latitude_save,1,size(Height_save,2));                       

            for l=1:numel(lat_mdpts)
                % calculate number of profiles
                lata=Latitude_save(:,1);
                ind=find(lata>=lat_bins(l) & lata<lat_bins(l+1));
                num_obs(l,t)=numel(ind);                              
                
                % calculate frequency of occurrence of L/M/H categories
                cloud_labelsa=cloud_labels_save(ind,:);
                H_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==1 & cloud_labelsa(:,2)==0 & cloud_labelsa(:,3)==0))./size(cloud_labelsa,1);
                M_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==0 & cloud_labelsa(:,2)==1 & cloud_labelsa(:,3)==0))./size(cloud_labelsa,1);
                L_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==0 & cloud_labelsa(:,2)==0 & cloud_labelsa(:,3)==1))./size(cloud_labelsa,1);
                LM_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==0 & cloud_labelsa(:,2)==1 & cloud_labelsa(:,3)==1))./size(cloud_labelsa,1);
                MH_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==1 & cloud_labelsa(:,2)==1 & cloud_labelsa(:,3)==0))./size(cloud_labelsa,1);
                LH_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==1 & cloud_labelsa(:,2)==0 & cloud_labelsa(:,3)==1))./size(cloud_labelsa,1);
                LMH_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==1 & cloud_labelsa(:,2)==1 & cloud_labelsa(:,3)==1))./size(cloud_labelsa,1);
                clear_frac(l,t)=numel(cloud_labelsa(cloud_labelsa(:,1)==0 & cloud_labelsa(:,2)==0 & cloud_labelsa(:,3)==0))./size(cloud_labelsa,1);       
                clearvars cloud_labelsa lata ind
                
                % calculate mean cloud fraction and heating rates
                for k=1:numel(z_mdpts)
                    ind=find(Latitude_save>=lat_bins(l) & Latitude_save<lat_bins(l+1) & ...
                        Height_save>=z_bins(k) & Height_save<z_bins(k+1));
                    
                    % don't use pressure weighting for cloud fraction since it is defined by volume
                    CloudFraction_gridded(l,k,t)=nanmean(CloudFraction_save(ind));                    
                    
                    % use pressure weighting for heating rates and temperature
                    LW_hr_gridded(l,k,t)=nansum(LW_hr_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    LW_hr_clr_gridded(l,k,t)=nansum(LW_hr_clr_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    LW_ACRE_gridded(l,k,t)=nansum((LW_hr_save(ind) - LW_hr_clr_save(ind)).*dp_save(ind))./nansum(dp_save(ind));
                    LW_up_clr_gridded(l,k,t)=nansum(LW_up_clr_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    LW_down_clr_gridded(l,k,t)=nansum(LW_down_clr_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    Temperature_gridded(l,k,t)=nansum(Temperature_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    Specific_humidity_gridded(l,k,t)=nansum(Specific_humidity_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    U_velocity_gridded(l,k,t)=nansum(U_velocity_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                    V_velocity_gridded(l,k,t)=nansum(V_velocity_save(ind).*dp_save(ind))./nansum(dp_save(ind));
                end
                                                                                                         
            end
            clearvars l k CloudFraction_save Height_save Latitude_save cloud_labels_save LW_hr_save LW_hr_clr_save ...
                dp_save Temperature_save Specific_humidity_save U_velocity_save V_velocity_save
        end
    end
                
    t
end
clearvars t

%%%%%% save data %%%%%%
time=time_mdpts;
lat=lat_mdpts;
z=z_mdpts;
time_edges=time_bins;
lat_edges=lat_bins;
z_edges=z_bins;
cloud_fraction=CloudFraction_gridded;
LW_hr=LW_hr_gridded;
LW_hr_clr=LW_hr_clr_gridded;
LW_ACRE=LW_ACRE_gridded;
LW_up_clr=LW_up_clr_gridded;
LW_down_clr=LW_down_clr_gridded;
T=Temperature_gridded;
Specific_humidity=Specific_humidity_gridded;
u=U_velocity_gridded;
v=U_velocity_gridded;
clearvars time_mdpts lat_mdpts z_mdpts time_bins lat_bins z_bins CloudFraction_gridded Temperature_gridded ...
    LW_hr_gridded LW_hr_clr_gridded LW_ACRE_gridded LW_up_clr_gridded LW_down_clr_gridded Specific_humidity_gridded ...
    U_velocity_gridded U_velocity_gridded

save('./data/CloudSat_pentad_zonal_mean_height_coords_480m.mat','time','time_edges','lat','lat_edges','z','z_edges','num_obs','cloud_fraction',...
    'T','Specific_humidity','u','v','LW_hr','LW_hr_clr','LW_ACRE','H_frac','M_frac','L_frac','LM_frac','LH_frac','MH_frac','LMH_frac','clear_frac',...
    'LW_up_clr','LW_down_clr')


%% load SAM index
load('./data/CloudSat_pentad_zonal_mean_height_coords_480m.mat','time','time_edges')

% load SAM
path='./data/ERA5_AM_19790101_20181231.nc';
time_SAM=ncread(path,'time') + datenum(1900,1,1);
SAM=double(ncread(path,'SAM'));
clearvars path

% compute pendad-mean values
SAM_p=NaN*ones(size(time));

for t=1:numel(time)
    time_ind=find(time_SAM>=time_edges(t) & time_SAM<time_edges(t+1));
    SAM_p(t)=mean(SAM(time_ind));
end
SAM=SAM_p;
clearvars t time_ind SAM_p time_SAM

save('./data/CloudSat_pentad_zonal_mean_height_coords_480m.mat','SAM','-append')






