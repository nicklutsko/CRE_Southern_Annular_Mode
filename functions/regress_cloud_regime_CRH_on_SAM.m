% function to compute regression of cloud-category LW CRH profile
% on SAM index. These values are used to decompose the overall
% cloud radaitive heating anomaly associated with SAM anomaly into
% contributions from different cloud regimes.
%
% Inputs:
% (1) lat, z, time: vetors of latitude, height, and time
% (2) LW_CRH_regime: regime-average LW CRH (units: K d-1; dimension: lat, z
% time)
% (3) SAM: standardized SAM index (units: sigma SAM index, dimensions: time)
% (4) num_obs: number of observations (dimensions: lat, time)
%
% Outputs:
% (1) dLW_CRH_regime_dSAM: regression coefficient for LW_CRH_regime regressed on
% SAM (single linear regression model, coefficients represent total derivatives) 
% (units: K d-1 sigma_SAM-1; dimensions: lat, z)
% (2) dLW_CRH_regime_dSAM_CI: 95% confidence interval for regression
% coefficients (units: K d-1 sigma_SAM-1; dimensions: lat, z)

function [dLW_CRH_regime_dSAM, dLW_CRH_regime_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_regime,SAM,num_obs)

% remove times with missing data
time_ind=[];
for n=1:numel(time)
    if any(num_obs(:,n)~=0)
        time_ind(end+1)=n;
    end
end
clearvars n

LW_CRH_regime=LW_CRH_regime(:,:,time_ind);
num_obs=num_obs(:,time_ind);
SAM=SAM(time_ind);
time=time(time_ind);
clearvars time_ind

% standardize SAM index
SAM=(SAM-mean(SAM))./std(SAM);

for l=1:numel(time)
    LW_CRH_regime_save(:,:,l)=LW_CRH_regime(:,:,l)-nanmean(LW_CRH_regime,3);
end
LW_CRH_regime=LW_CRH_regime_save;
clearvars LW_CRH_regime_save l

% regress LW_CRH_regime on SAM index
dLW_CRH_regime_dSAM=NaN*ones(numel(lat),numel(z));
dLW_CRH_regime_dSAM_CI=NaN*ones(numel(lat),numel(z));

for n=1:numel(lat)
    for l=1:numel(z)
        LW_CRH_regimea=transpose(squeeze(LW_CRH_regime(n,l,:)));
        num_obsa=num_obs(n,:);

        ind=find(num_obsa>0 & ~isnan(LW_CRH_regimea));
        LW_CRH_regimea=LW_CRH_regimea(ind);
        SAMa=SAM(ind);
        LW_CRH_regimea=reshape(LW_CRH_regimea,numel(LW_CRH_regimea),1);
        SAMa=reshape(SAMa,numel(SAMa),1);
        clearvars ind num_obsa
        
        if isempty(LW_CRH_regimea) % no available data for this cloud regime
            dLW_CRH_regime_dSAM(n,l)=NaN;
            dLW_CRH_retime_dSAM_CI(n,l)=NaN;
        else
            mdl=fitlm(SAMa,LW_CRH_regimea);
            dLW_CRH_regime_dSAM(n,l)=mdl.Coefficients.Estimate(2);
            dLW_CRH_regime_dSAM_CI(n,l)=compute_regression_slope_95_CI(LW_CRH_regimea,SAMa,mdl.Coefficients.SE(2));
            clearvars mdl LW_CRH_regimea SAMa
        end
    
    end
end
clearvars n




















