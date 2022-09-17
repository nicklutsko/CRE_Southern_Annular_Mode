% Compute climatology of atmospheric LW cloud radaitive heating (R) and
% cloud fraction, and compute regressions of R and cloud fraction against
% the SAM index. Note: "load_CloudSat_heigt_coordinates_480m.m" needs to be run
% before this script.
clearvars; clc;
addpath ./functions/

load('./data/CloudSat_pentad_zonal_mean_height_coords_480m.mat','lat','lat_edges','LW_ACRE','cloud_fraction','z','z_edges','SAM','num_obs','time')

% remove times with missing data
time_ind=[];
for n=1:numel(time)
    if any(num_obs(:,n)~=0)
        time_ind(end+1)=n;
    end
end
clearvars n

LW_ACRE=LW_ACRE(:,:,time_ind);
cloud_fraction=cloud_fraction(:,:,time_ind);
SAM=SAM(time_ind);
num_obs=num_obs(:,time_ind);
time=time(time_ind);
clearvars time_ind

% compute climatology
LW_ACRE_mean=NaN.*ones(numel(lat),numel(z));
cloud_fraction_mean=NaN.*ones(numel(lat),numel(z));

for n=1:numel(lat)
    for m=1:numel(z)
        LW_ACREa=squeeze(LW_ACRE(n,m,:));
        cloud_fractiona=squeeze(cloud_fraction(n,m,:));
        Wa=squeeze(num_obs(n,:));
        Wa=Wa./max(Wa);
        
        LW_ACRE_mean(n,m)=nansum(LW_ACREa(:).*Wa(:))./nansum(Wa(:));
        cloud_fraction_mean(n,m)=nansum(cloud_fractiona(:).*Wa(:))./nansum(Wa(:));
        clearvars LW_ACREa cloud_fractiona Wa
    end
end
clearvars n m

% compute regression of cloud properties on SAM index
SAM=(SAM-mean(SAM))./std(SAM); % standardize SAM index

% remove means
LW_ACRE_save=NaN*ones(size(LW_ACRE));
cloud_fraction_save=NaN*ones(size(cloud_fraction));
for l=1:numel(time)
    LW_ACRE_save(:,:,l)=LW_ACRE(:,:,l)-nanmean(LW_ACRE,3);
    cloud_fraction_save(:,:,l)=cloud_fraction(:,:,l)-nanmean(cloud_fraction,3);
end
LW_ACRE=LW_ACRE_save;
cloud_fraction=cloud_fraction_save;
clearvars LW_ACRE_save cloud_fraction_save l

% compute regressions
dLW_ACRE_dSAM=NaN*ones(numel(lat),numel(z)); % regression slope central estimate
dLW_ACRE_dSAM_CI=NaN*ones(numel(lat),numel(z)); % regression slope central estimate
dcloud_fraction_dSAM=NaN*ones(numel(lat),numel(z)); % regression slope 95% CI
dcloud_fraction_dSAM_CI=NaN*ones(numel(lat),numel(z)); % regression slope 95% CI
for n=1:numel(lat)
    for m=1:numel(z)
        LW_ACREa=squeeze(LW_ACRE(n,m,:));
        cloud_fractiona=squeeze(cloud_fraction(n,m,:));
        num_obsa=squeeze(num_obs(n,:));
        
        ind=find(num_obsa>0);
        LW_ACREa=LW_ACREa(ind);
        cloud_fractiona=cloud_fractiona(ind);
        SAMa=transpose(SAM(ind));
        
        % LW CRH
        mdl=fitlm(SAMa,LW_ACREa);
        dLW_ACRE_dSAM(n,m)=mdl.Coefficients.Estimate(2);
        dLW_ACRE_dSAM_CI(n,m) = compute_regression_slope_95_CI(LW_ACREa,SAMa,mdl.Coefficients.SE(2));
        clearvars mdl
        
        % cloud fraction
        mdl=fitlm(SAMa,cloud_fractiona);
        dcloud_fraction_dSAM(n,m)=mdl.Coefficients.Estimate(2);
        dcloud_fraction_dSAM_CI(n,m) = compute_regression_slope_95_CI(cloud_fractiona,SAMa,mdl.Coefficients.SE(2));
        clearvars mdl
        
        clearvars LW_ACREa cloud_fractiona SAMa num_obsa ind
    end
    n
end
clearvars n m

