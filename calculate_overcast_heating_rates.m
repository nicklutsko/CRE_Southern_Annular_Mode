% This script calculates the atmospheric LW cloud radiative heating
% rate averaged over each cloud regime (L, M, H, LM, etc.).
% The script assumes that the user has saved CloudSat data
% products 2B-FLXHR-LIDAR, 2B-GEOPROF-LIDAR, and ECMWF-AUX in separate
% directories, and for each data product, all data files are saved in the same directory. 
% The user will need to change the directory sturcture in the script below to match their own machine.
% Note: this script needs to be run before "LW_CRH_decomposition.m"
clearvars; clc;
addpath ./functions/

% parameters for latitude-time binning
lat_bins=(-82.5):2.5:2.5; lat_mdpts=(lat_bins(1:(end-1))+lat_bins(2:end))./2;
time_bins=datenum(2007,1,1):5:datenum(2010,12,31); time_mdpts=(time_bins(1:(end-1))+time_bins(2:end))./2;
z_bins=0:0.48:(18+.48.*2); z_mdpts=(z_bins(1:(end-1))+z_bins(2:end))./2;
lat_min=-90; lat_max=2.5;

% constants
g=9.81; % graviatational constant (m/s^2)
cp=1004; % specific heat capacity of dry air at constant pressure (J/K kg)

% initialize arrays
LW_CRH_L=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_M=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_H=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_LM=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_LH=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_MH=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
LW_CRH_LMH=NaN*ones(numel(lat_mdpts),numel(z_mdpts),numel(time_mdpts));
num_obs=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_L=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_M=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_H=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_LM=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_LH=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_MH=NaN*ones(numel(lat_mdpts),numel(time_mdpts));
num_obs_LMH=NaN*ones(numel(lat_mdpts),numel(time_mdpts));

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
        LW_CRH_L(:,:,t)=NaN;
        LW_CRH_M(:,:,t)=NaN;
        LW_CRH_H(:,:,t)=NaN;
        LW_CRH_LM(:,:,t)=NaN;
        LW_CRH_LH(:,:,t)=NaN;
        LW_CRH_MH(:,:,t)=NaN;
        LW_CRH_LMH(:,:,t)=NaN;
        num_obs(:,t)=0;
        num_obs_L(:,t)=0;
        num_obs_M(:,t)=0;
        num_obs_H(:,t)=0;
        num_obs_LM(:,t)=0;
        num_obs_LH(:,t)=0;
        num_obs_MH(:,t)=0;
        num_obs_LMH(:,t)=0;          
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
                FD_NC = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FD_NC');
                LW_down_clr = double(squeeze(FD_NC(2,time_ind,:)))./10;
                FU_NC = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/FU_NC');
                LW_up_clr = double(squeeze(FU_NC(2,time_ind,:)))./10;
                QR = hdfread(path, '/2B-FLXHR-LIDAR/Data Fields/QR');
                LW_hr = double(squeeze(QR(2,time_ind,:)))./100;
                clearvars FD_NC FU_NC QR                                                
                
                % fill in missing/bad data with NaN
                CloudFraction(CloudFraction<0)=NaN; LayerTop(LayerTop<0)=NaN; LayerBase(LayerBase<0)=NaN; LW_down_clr(LW_down_clr==-999)=NaN;     
                LW_up_clr(LW_up_clr==-999)=NaN; LW_hr(LW_hr==-999./100)=NaN; Height(Height==-9999./1000)=NaN;
                
                % load ECMWF-AUX data
                path=['/data/cawall/CloudSat/ECMWF_AUX/' fnames_ECMWF{n}];  % the user will need to change this path              
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
                clearvars Pressure_half_int path time_ind LW_up_clr LW_down_clr Pressure

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
                    
                    current_index_start=1;
                    current_index_end=numel(Latitude);
                    Latitude_save(1,current_index_start:current_index_end)=Latitude;
                    Height_save(current_index_start:current_index_end,:)=Height;  
                    CloudFraction_save(current_index_start:current_index_end,:)=CloudFraction;
                    cloud_labels_save(current_index_start:current_index_end,:)=cloud_labels;
                    dp_save(current_index_start:current_index_end,:)=dp;
                    LW_hr_save(current_index_start:current_index_end,:)=LW_hr;
                    LW_hr_clr_save(current_index_start:current_index_end,:)=LW_hr_clr;                                       
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
                end
                                            
            end
            
            if t==1
                n
            end
        end
              
        if ~exist('current_index_end') % no available data for this pentad
            LW_CRH_L(:,:,t)=NaN;
            LW_CRH_M(:,:,t)=NaN;
            LW_CRH_H(:,:,t)=NaN;
            LW_CRH_LM(:,:,t)=NaN;
            LW_CRH_LH(:,:,t)=NaN;
            LW_CRH_MH(:,:,t)=NaN;
            LW_CRH_LMH(:,:,t)=NaN;
            num_obs(:,t)=0;
            num_obs_L(:,t)=0;
            num_obs_M(:,t)=0;
            num_obs_H(:,t)=0;
            num_obs_LM(:,t)=0;
            num_obs_LH(:,t)=0;
            num_obs_MH(:,t)=0;
            num_obs_LMH(:,t)=0;        
        else
        
            Latitude_save=Latitude_save(1:current_index_end);
            Height_save=Height_save(1:current_index_end,:);  
            CloudFraction_save=CloudFraction_save(1:current_index_end,:);  
            cloud_labels_save=cloud_labels_save(1:current_index_end,:);
            dp_save=dp_save(1:current_index_end,:);
            LW_hr_save=LW_hr_save(1:current_index_end,:);
            LW_hr_clr_save=LW_hr_clr_save(1:current_index_end,:);                       
            clearvars Height CloudFraction Latitude current_index_start current_index_end nelements n nfiles ...
                fnames_CloudSat cloud_labels dp LW_hr LW_hr_clr   
            
            % calculate flag for available data
            available_data_flag=NaN*ones(size(Latitude_save));
            for k1=1:numel(Latitude_save)
                if all(isnan(Height_save(k1,:)))
                    available_data_flag(k1)=0;
                else
                    available_data_flag(k1)=1;
                end
            end
            clearvars k1
            
            % bin by latitude and height
            for l=1:numel(lat_mdpts)
                % calculate number of profiles
                ind=find(Latitude_save>=lat_bins(l) & Latitude_save<lat_bins(l+1) & available_data_flag==1);
                cloud_labelsa=cloud_labels_save(ind,:);
                Heighta=Height_save(ind,:);
                CloudFractiona=CloudFraction_save(ind,:);
                LW_hra=LW_hr_save(ind,:);
                LW_hr_clra=LW_hr_clr_save(ind,:);
                dpa=dp_save(ind,:);
                num_obs(l,t)=numel(ind); 
                clearvars ind
                
                [LW_CRH_La,LW_CRH_Ma,LW_CRH_Ha,LW_CRH_LMa,LW_CRH_LHa,LW_CRH_MHa,LW_CRH_LMHa,num_obs_La,...
                    num_obs_Ma,num_obs_Ha,num_obs_LMa,num_obs_LHa,num_obs_MHa,num_obs_LMHa] = calculate_overcast_CRH_pentad(cloud_labelsa,...
                    Heighta,LW_hra,LW_hr_clra,dpa,z_bins,z_mdpts);
                LW_CRH_L(l,:,t)=LW_CRH_La;
                LW_CRH_M(l,:,t)=LW_CRH_Ma;
                LW_CRH_H(l,:,t)=LW_CRH_Ha;
                LW_CRH_LM(l,:,t)=LW_CRH_LMa;
                LW_CRH_LH(l,:,t)=LW_CRH_LHa;
                LW_CRH_MH(l,:,t)=LW_CRH_MHa;
                LW_CRH_LMH(l,:,t)=LW_CRH_LMHa;
                num_obs_L(l,t)=num_obs_La;
                num_obs_M(l,t)=num_obs_Ma;
                num_obs_H(l,t)=num_obs_Ha;
                num_obs_LM(l,t)=num_obs_LMa;
                num_obs_LH(l,t)=num_obs_LHa;
                num_obs_MH(l,t)=num_obs_MHa;
                num_obs_LMH(l,t)=num_obs_LMHa;
                clearvars LW_CRH_La LW_CRH_Ma LW_CRH_Ha LW_CRH_LMa LW_CRH_LHa LW_CRH_MHa LW_CRH_LMHa num_obs_La num_obs_Ma ...
                    num_obs_Ha num_obs_LMa num_obs_LHa num_obs_MHa num_obs_LMHa
                
            end                  
            clearvars l k Height_save Latitude_save cloud_labels_save LW_hr_save LW_hr_clr_save ...
                dp_save
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
clearvars time_mdpts lat_mdpts z_mdpts time_bins lat_bins z_bins

save('./data/CloudSat_overcast_CRH_480m.mat','time','time_edges','lat','lat_edges','z','z_edges','num_obs',...
    'LW_CRH_L','LW_CRH_M','LW_CRH_H','LW_CRH_LM','LW_CRH_LH','LW_CRH_MH','LW_CRH_LMH',...
    'num_obs_L','num_obs_M','num_obs_H','num_obs_LM','num_obs_LH','num_obs_MH','num_obs_LMH','-v7.3')


                
                
                
                
                
                

