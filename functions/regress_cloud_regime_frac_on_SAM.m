% function to compute regression of cloud-category fractional occurrence
% on SAM index. These values are used to decompose the overall
% cloud radaitive heating anomaly associated with SAM anomaly into
% contributions from different cloud regimes.
%
% Inputs:
% (1) lat, time: vetors of latitude and time
% (2) cloud_frac: fractional occurrence of cloud regime (units: %,
% dimensions: lat, time)
% (3) SAM: standardized SAM index (units: sigma SAM index, dimensions: time)
% (4) num_obs: number of observations (dimensions: lat, time)
%
% Outputs:
% (1) dcloud_frac_dSAM: regression coefficient for cloud_frac regressed on
% SAM (single linear regression model, coefficients represent total derivatives) 
% (units: % sigma_SAM-1; dimensions: lat)
% (2) dcloud_frac_dSAM_CI: 95% confidence interval for regression
% coefficients (units: % sigma_SAM-1; dimensions: lat)

function [dcloud_frac_dSAM, dcloud_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,cloud_frac,SAM,num_obs)

% remove times with missing data
time_ind=[];
for n=1:numel(time)
    if any(num_obs(:,n)~=0)
        time_ind(end+1)=n;
    end
end
clearvars n

cloud_frac=cloud_frac(:,time_ind);
num_obs=num_obs(:,time_ind);
SAM=SAM(time_ind);
time=time(time_ind);
clearvars time_ind

% standardize SAM index
SAM=(SAM-mean(SAM))./std(SAM);

for l=1:numel(time)
    cloud_frac_save(:,l)=cloud_frac(:,l)-nanmean(cloud_frac,2);
end
cloud_frac=cloud_frac_save;
clearvars cloud_frac_save l

% regress cloud_frac on SAM index
dcloud_frac_dSAM=NaN*ones(size(lat));
dcloud_frac_dSAM_CI=NaN*ones(size(lat));

for n=1:numel(lat)
    cloud_fraca=cloud_frac(n,:);
    num_obsa=num_obs(n,:);
    
    ind=find(num_obsa>0);
    cloud_fraca=transpose(cloud_fraca(ind));
    SAMa=transpose(SAM(ind));
    clearvars ind num_obsa
    
    mdl=fitlm(SAMa,cloud_fraca);
    dcloud_frac_dSAM(n)=mdl.Coefficients.Estimate(2);
    dcloud_frac_dSAM_CI(n)=compute_regression_slope_95_CI(cloud_fraca,SAMa,mdl.Coefficients.SE(2));
    clearvars mdl cloud_fraca SAMa
end
clearvars n








