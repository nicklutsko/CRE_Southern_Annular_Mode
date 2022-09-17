% compute 95% confidence intervals for regime-average cloud fraction and LW
% CRH using bootstrapping.
% 
% Inputs:
% (1) LW_CRH_regime: LW CRH for cloud regime (units: K d-1; dimension: lat, z,
% time)
% (2) f_regime: fractional occurrence of cloud regime (units: 0-1;
% dimensions: lat,time)
% (3) num_obs_regime: number of obervations of regime (dimensions: lat,
% time)
% (4) lat, z, time: vectors of latitude, height, and time
%
% Outputs:
% (1) LW_CRH_mean_CI: 95% confidence interval regime-avearge LW CRH (i.e. confidence interval is LW_CRH_mean +/- LW_CRH_mean_CI)(units: K d-1)
% (2) f_mean_CI: 95% confidence interval for mean fractional occurrence of cloud regime (units: 0-1)

function [LW_CRH_mean_CI,f_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_regime,f_regime,num_obs_regime,lat,z,time)

% select random samples with replacement
n_bootstraps=1000; % number of bootstrap samples
LW_CRH_bootstraps=NaN*ones(numel(lat),numel(z),n_bootstraps);
f_bootstraps=NaN*ones(numel(lat),n_bootstraps);

for n=1:n_bootstraps
    ind=randi(numel(time),1,numel(time)); % random sample of time points chosen with replacement
    [LW_CRH_bootstraps(:,:,n),f_bootstraps(:,n)]=regime_average_CRH_cloud_frac(LW_CRH_regime(:,:,ind),f_regime(:,ind),num_obs_regime(:,ind),lat,z);
    clearvars ind
    
    if mod(n,10)==0
        n
    end
end
clearvars n

% compute confidence intervals
LW_CRH_mean_CI=NaN*ones(numel(lat),numel(z));
f_mean_CI=NaN*ones(size(lat));

for n=1:numel(lat)
    fa=squeeze(f_bootstraps(n,:));
    f_mean_CI(n)=(quantile(fa,0.975)-quantile(fa,0.025))./2;
    
    for m=1:numel(z)
        LW_CRHa=squeeze(LW_CRH_bootstraps(n,m,:));
        LW_CRH_mean_CI(n,m)=(quantile(LW_CRHa,0.975)-quantile(LW_CRHa,0.025))./2;
    end
end
clearvars fa LW_CRHa n
        
        
        



    
    

