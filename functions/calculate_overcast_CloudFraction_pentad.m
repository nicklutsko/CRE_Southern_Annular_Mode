% calculate cloud-fraction vertical profile for cloud regimes L, M, H, etc.
% assumes that Height array is spaced by 240 m in the vertical
% Inputs:
% (1) cloud_labels: in second dimension, 1 in first entry mean H, 1 in
% second entry mean M, 1 in third entry means L (dimensions: footprint #, 3)
% (2) Height: array of heights (dimensions: footprint #, height; units: m)
% (3) CloudFraction: cloud-fraction verticla profile (dimensions: footprint #, height;
% units: %)
% 
% Outputs
% (1) CloudFraction_L, CloudFraction_M, CloudFraction_H, etc.: 
% cloud-fraction profiles for L, M, H, etc. cloud regimes
% (units: %; NaN if cloud regime is not observed)
% (2) num_obs_L, num_obs_M, num_obs_H, etc.: number of observations for
% each cloud category

function [CloudFraction_L,CloudFraction_M,CloudFraction_H,CloudFraction_LM,CloudFraction_LH,CloudFraction_MH,CloudFraction_LMH,num_obs_L,...
    num_obs_M,num_obs_H,num_obs_LM,num_obs_LH,num_obs_MH,num_obs_LMH] = calculate_overcast_CloudFraction_pentad(cloud_labels,Height,CloudFraction,...
    z_bins,z_mdpts)

% change arrays so that the final index is the lowest level above 0 km elevation
Height_save=NaN*ones(size(Height,1),100);
CloudFraction_save=NaN*ones(size(Height,1),100);

for n=1:size(Height,1)
    Heighta=Height(n,:);
    z0_ind=find(Heighta==min(Heighta(Heighta>0)));
    Height_save(n,:)=Height(n,(z0_ind-100+1):z0_ind);
    CloudFraction_save(n,:)=CloudFraction(n,(z0_ind-100+1):z0_ind);
end
Height=Height_save; CloudFraction=CloudFraction_save;
clearvars Heighta n z0_ind Height_save CloudFraction_save

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

%%%%% mean cloud fraction profiles %%%%%
[CloudFraction_H, num_obs_H] = calculate_non_weighted_mean(CloudFraction,H_ind,z_ind,z_mdpts);
[CloudFraction_M, num_obs_M] = calculate_non_weighted_mean(CloudFraction,M_ind,z_ind,z_mdpts);
[CloudFraction_L, num_obs_L] = calculate_non_weighted_mean(CloudFraction,L_ind,z_ind,z_mdpts);
[CloudFraction_LM, num_obs_LM] = calculate_non_weighted_mean(CloudFraction,LM_ind,z_ind,z_mdpts);
[CloudFraction_LH, num_obs_LH] = calculate_non_weighted_mean(CloudFraction,LH_ind,z_ind,z_mdpts);
[CloudFraction_MH, num_obs_MH] = calculate_non_weighted_mean(CloudFraction,MH_ind,z_ind,z_mdpts);
[CloudFraction_LMH, num_obs_LMH] = calculate_non_weighted_mean(CloudFraction,LMH_ind,z_ind,z_mdpts);
















