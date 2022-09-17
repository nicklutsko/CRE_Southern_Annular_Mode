% calculate overcast heating rates for cloud regimes L, M, H, etc.
% assumes that Height array is spaced by 240 m in the vertical
% Inputs:
% (1) cloud_labels: in second dimension, 1 in first entry mean H, 1 in
% second entry mean M, 1 in third entry means L (dimensions: footprint #, 3)
% (2) Height: array of heights (dimensions: footprint #, height; units: m)
% (3) LW_hr, LW_hr_clr: LW heating rates (dimensions: footprint #, height;
% units: Kd-1)
% 
% Outputs
% (1) LW_CRH_L, LW_CRH_M, LW_CRH_H, LW_CRH_LM, LW_CRH_LH, LW_CRH_MH,
% LW_CRH_LMH: overcast heating rates for L, M, H, etc. cloud regimes
% (units: Kd-1; NaN if cloud regime is not observed)
% (2) num_obs_L, num_obs_M, num_obs_H, etc.: number of observations for
% each cloud category

function [LW_CRH_L,LW_CRH_M,LW_CRH_H,LW_CRH_LM,LW_CRH_LH,LW_CRH_MH,LW_CRH_LMH,num_obs_L,...
    num_obs_M,num_obs_H,num_obs_LM,num_obs_LH,num_obs_MH,num_obs_LMH] = calculate_overcast_CRH_pentad(cloud_labels,Height,LW_hr,LW_hr_clr,...
    dp,z_bins,z_mdpts)

% change arrays so that the final index is the lowest level above 0 km elevation
Height_save=NaN*ones(size(Height,1),100);
LW_hr_save=NaN*ones(size(Height,1),100);
LW_hr_clr_save=NaN*ones(size(Height,1),100);
dp_save=NaN*ones(size(Height,1),100);

for n=1:size(Height,1)
    Heighta=Height(n,:);
    z0_ind=find(Heighta==min(Heighta(Heighta>0)));
    Height_save(n,:)=Height(n,(z0_ind-100+1):z0_ind);
    LW_hr_save(n,:)=LW_hr(n,(z0_ind-100+1):z0_ind);
    LW_hr_clr_save(n,:)=LW_hr_clr(n,(z0_ind-100+1):z0_ind);
    dp_save(n,:)=dp(n,(z0_ind-100+1):z0_ind);
end
Height=Height_save; LW_hr=LW_hr_save; LW_hr_clr=LW_hr_clr_save; dp=dp_save;
clearvars Heighta n z0_ind Height_save LW_hr_save LW_hr_clr_save dp_save

% indices for vertical averaging - assumes same vertical spacing for all
% profiles
for n=1:numel(z_mdpts)
    Heighta=Height(1,:);
    z_ind{n}=find(Heighta>=z_bins(n) & Heighta<z_bins(n+1));
end
clearvars Heighta

% indices for cloud categories
H_ind=find(cloud_labels(:,1)==1 & cloud_labels(:,2)==0 & cloud_labels(:,3)==0);
M_ind=find(cloud_labels(:,1)==0 & cloud_labels(:,2)==1 & cloud_labels(:,3)==0);
L_ind=find(cloud_labels(:,1)==0 & cloud_labels(:,2)==0 & cloud_labels(:,3)==1);
LM_ind=find(cloud_labels(:,1)==0 & cloud_labels(:,2)==1 & cloud_labels(:,3)==1);
LH_ind=find(cloud_labels(:,1)==1 & cloud_labels(:,2)==0 & cloud_labels(:,3)==1);
MH_ind=find(cloud_labels(:,1)==1 & cloud_labels(:,2)==1 & cloud_labels(:,3)==0);
LMH_ind=find(cloud_labels(:,1)==1 & cloud_labels(:,2)==1 & cloud_labels(:,3)==1);

%%%%% overcast heating rates %%%%%
[LW_CRH_H, num_obs_H] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,H_ind,z_ind,z_mdpts);
[LW_CRH_M, num_obs_M] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,M_ind,z_ind,z_mdpts);
[LW_CRH_L, num_obs_L] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,L_ind,z_ind,z_mdpts);
[LW_CRH_LM, num_obs_LM] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,LM_ind,z_ind,z_mdpts);
[LW_CRH_LH, num_obs_LH] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,LH_ind,z_ind,z_mdpts);
[LW_CRH_MH, num_obs_MH] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,MH_ind,z_ind,z_mdpts);
[LW_CRH_LMH, num_obs_LMH] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,LMH_ind,z_ind,z_mdpts);

