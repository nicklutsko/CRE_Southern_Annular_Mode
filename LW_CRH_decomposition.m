% Decompose cloud radiative heating anomalies associated with SAM variability
% into contributions from different cloud regimes (L, M, H, LM, etc.)
% Note: "calculate_overcast_heating_rates.m" and
% "load_CloudSat_height_coordinates_480m.mat" need to be run before this
% script.
clearvars; clc;
addpath ./functions/

%%%%%%%%%%%%%%% calculte CRH decomposition %%%%%%%%%%%%%%%%%
% load data
load('./data/CloudSat_pentad_zonal_mean_height_coords_480m.mat','lat','lat_edges','LW_ACRE','z','z_edges','SAM',...
    'num_obs','time','L_frac','M_frac','H_frac','LM_frac','LH_frac','MH_frac','LMH_frac','clear_frac')
d=load('./data/CloudSat_overcast_CRH_480m.mat','lat','z','time');
if any(d.lat(:)~=lat(:)) | any(d.z(:)~=z(:)) | any(d.time(:)~=time(:))
    'error: dimensions do not match'
else
    load('./data/CloudSat_overcast_CRH_480m.mat','LW_CRH_L','LW_CRH_LM','LW_CRH_LMH','LW_CRH_LH','LW_CRH_M','LW_CRH_H','LW_CRH_MH',...
        'num_obs_L','num_obs_LM','num_obs_LMH','num_obs_LH','num_obs_M','num_obs_H','num_obs_MH')
    num_obs_clear=num_obs - (num_obs_L + num_obs_LM + num_obs_LMH + num_obs_MH + num_obs_LH + num_obs_M + num_obs_H);
    LW_CRH_clear=zeros(size(LW_CRH_L));
end

% compute regime averages
[LW_CRH_L_mean,L_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_L,L_frac,num_obs_L,lat,z);
[LW_CRH_LM_mean,LM_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_LM,LM_frac,num_obs_LM,lat,z);
[LW_CRH_LMH_mean,LMH_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_LMH,LMH_frac,num_obs_LMH,lat,z);
[LW_CRH_LH_mean,LH_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_LH,LH_frac,num_obs_LH,lat,z);
[LW_CRH_MH_mean,MH_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_MH,MH_frac,num_obs_MH,lat,z);
[LW_CRH_M_mean,M_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_M,M_frac,num_obs_M,lat,z);
[LW_CRH_H_mean,H_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_H,H_frac,num_obs_H,lat,z);
[LW_CRH_clear_mean,clear_frac_mean]=regime_average_CRH_cloud_frac(LW_CRH_clear,clear_frac,num_obs_clear,lat,z);
'regime averages done'

% compute confidence intervals for regime averages
[LW_CRH_L_mean_CI,L_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_L,L_frac,num_obs_L,lat,z,time);
[LW_CRH_LM_mean_CI,LM_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_LM,LM_frac,num_obs_LM,lat,z,time);
[LW_CRH_LMH_mean_CI,LMH_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_LMH,LMH_frac,num_obs_LMH,lat,z,time);
[LW_CRH_LH_mean_CI,LH_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_LH,LH_frac,num_obs_LH,lat,z,time);
[LW_CRH_MH_mean_CI,MH_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_MH,MH_frac,num_obs_MH,lat,z,time);
[LW_CRH_M_mean_CI,M_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_M,M_frac,num_obs_M,lat,z,time);
[LW_CRH_H_mean_CI,H_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_H,H_frac,num_obs_H,lat,z,time);
[LW_CRH_clear_mean_CI,clear_frac_mean_CI]=regime_average_CRH_cloud_frac_CI(LW_CRH_clear,clear_frac,num_obs_clear,lat,z,time);
'regime average CIs done'

% regress cloud-regime fraction on SAM index
[dL_frac_dSAM, dL_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,L_frac,SAM,num_obs);
[dM_frac_dSAM, dM_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,M_frac,SAM,num_obs);
[dH_frac_dSAM, dH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,H_frac,SAM,num_obs);
[dLM_frac_dSAM, dLM_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LM_frac,SAM,num_obs);
[dLH_frac_dSAM, dLH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LH_frac,SAM,num_obs);
[dMH_frac_dSAM, dMH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,MH_frac,SAM,num_obs);
[dLMH_frac_dSAM, dLMH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LMH_frac,SAM,num_obs);
[dclear_frac_dSAM, dclear_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,clear_frac,SAM,num_obs);
'cloud fraction regression done'

% regress LW CRH profile on SAM index
[dLW_CRH_L_dSAM, dLW_CRH_L_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_L,SAM,num_obs);
[dLW_CRH_LM_dSAM, dLW_CRH_LM_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_LM,SAM,num_obs);
[dLW_CRH_LMH_dSAM, dLW_CRH_LMH_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_LMH,SAM,num_obs);
[dLW_CRH_LH_dSAM, dLW_CRH_LH_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_LH,SAM,num_obs);
[dLW_CRH_MH_dSAM, dLW_CRH_MH_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_MH,SAM,num_obs);
[dLW_CRH_M_dSAM, dLW_CRH_M_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_M,SAM,num_obs);
[dLW_CRH_H_dSAM, dLW_CRH_H_dSAM_CI] = regress_cloud_regime_CRH_on_SAM(lat,z,time,LW_CRH_H,SAM,num_obs);
dLW_CRH_clear_dSAM=zeros(size(dLW_CRH_L_dSAM)); dLW_CRH_clear_dSAM_CI=zeros(size(dLW_CRH_L_dSAM_CI));
'CRH regression done'

% organize data into cells
regime_label={'L','M','H','LM','LH','MH','LMH','clear'};
dLW_CRH_regime_dSAM={dLW_CRH_L_dSAM,dLW_CRH_M_dSAM,dLW_CRH_H_dSAM,dLW_CRH_LM_dSAM,dLW_CRH_LH_dSAM,dLW_CRH_MH_dSAM,dLW_CRH_LMH_dSAM,dLW_CRH_clear_dSAM};
dLW_CRH_regime_dSAM_CI={dLW_CRH_L_dSAM_CI,dLW_CRH_M_dSAM_CI,dLW_CRH_H_dSAM_CI,dLW_CRH_LM_dSAM_CI,dLW_CRH_LH_dSAM_CI,...
    dLW_CRH_MH_dSAM_CI,dLW_CRH_LMH_dSAM_CI,dLW_CRH_clear_dSAM_CI};
df_regime_dSAM={dL_frac_dSAM,dM_frac_dSAM,dH_frac_dSAM,dLM_frac_dSAM,dLH_frac_dSAM,dMH_frac_dSAM,dLMH_frac_dSAM,dclear_frac_dSAM};
df_regime_dSAM_CI={dL_frac_dSAM_CI,dM_frac_dSAM_CI,dH_frac_dSAM_CI,dLM_frac_dSAM_CI,dLH_frac_dSAM_CI,dMH_frac_dSAM_CI,dLMH_frac_dSAM_CI,dclear_frac_dSAM_CI};
LW_CRH_regime_mean={LW_CRH_L_mean,LW_CRH_M_mean,LW_CRH_H_mean,LW_CRH_LM_mean,LW_CRH_LH_mean,LW_CRH_MH_mean,LW_CRH_LMH_mean,LW_CRH_clear_mean};
LW_CRH_regime_mean_CI={LW_CRH_L_mean_CI,LW_CRH_M_mean_CI,LW_CRH_H_mean_CI,LW_CRH_LM_mean_CI,LW_CRH_LH_mean_CI,LW_CRH_MH_mean_CI,LW_CRH_LMH_mean_CI,LW_CRH_clear_mean_CI};
f_regime_mean={L_frac_mean,M_frac_mean,H_frac_mean,LM_frac_mean,LH_frac_mean,MH_frac_mean,LMH_frac_mean,clear_frac_mean};
f_regime_mean_CI={L_frac_mean_CI,M_frac_mean_CI,H_frac_mean_CI,LM_frac_mean_CI,LH_frac_mean_CI,MH_frac_mean_CI,LMH_frac_mean_CI,clear_frac_mean_CI};

%%%%% compute contribution of cloud regime to dR/ds %%%%%%
% This example shows how to compute the contribution of L clouds to dR/ds,
% where R is atmospheric LW cloud radaitive heating and s is the SAM index.
% You can change 'L' below to 'M', 'H', etc. to compute the contribution of
% another cloud regime to dR/ds.
regime_ind=find(strcmp(regime_label,regime_to_plot{case_num},'L')==1);

df_dSAM=df_regime_dSAM{regime_ind};
df_dSAM_CI=df_regime_dSAM_CI{regime_ind};
LW_CRH_mean=LW_CRH_regime_mean{regime_ind};
LW_CRH_mean_CI=LW_CRH_regime_mean_CI{regime_ind};
f_mean=f_regime_mean{regime_ind};
f_mean_CI=f_regime_mean_CI{regime_ind};
dLW_CRH_dSAM=dLW_CRH_regime_dSAM{regime_ind};
dLW_CRH_dSAM_CI=dLW_CRH_regime_dSAM_CI{regime_ind};

dCRH_dSAM_tot=NaN*ones(numel(lat),numel(z)); % contribution of cloud regime to dR/ds
dCRH_dSAM_tot_CI=NaN*ones(numel(lat),numel(z)); % 95% confidence interval for contribution of cloud regime to dR/ds

for n=1:numel(lat)
    for m=1:numel(z)
        dCRH_dSAM_tot(n,m) = df_dSAM(n).*LW_CRH_mean(n,m) + f_mean(n).*dLW_CRH_dSAM(n,m);

        dCRH_dSAM_tot_CI(n,m) = sqrt((df_dSAM_CI(n).^2).*LW_CRH_mean(n,m).^2 + ...
            (LW_CRH_mean_CI(n,m).^2).*df_dSAM(n).*2 + ...
            (f_mean_CI(n).^2).*dLW_CRH_dSAM(n,m).^2 + ...
            (dLW_CRH_dSAM_CI(n,m).^2).*f_mean(n).^2);
    end
end
clearvars n m




