% function to calculate non-weighted mean cloud fraction
% Inputs:
% (1) CloudFraction : CloudFraction vertical profile (dimensions: footprint #, height;
% units: %)
% (2) lat_ind: index indntifying footprints with a particular cloud
% category L, M, H, etc. (dimensions: footprint #)
% (3) z_ind: index for height averaging
% (4) z_mdpts: midpoints for height averaging
% 
% Outputs:
% (1) CloudFraction_mean: avearge LW cloud radiative heating for footprints from
% index "lat_ind" and vertical intervals corresponding to indices "z_ind"
% (units: %; dimensions: z_mdpts)
% (2) num_obs: number of observations used in average

function [CloudFraction_mean, num_obs] = calculate_non_weighted_mean(CloudFraction,lat_ind,z_ind,z_mdpts)

for n=1:numel(z_mdpts)
    if isempty(lat_ind)
        CloudFraction_mean(n)=NaN;
    else
        CloudFractiona=CloudFraction(lat_ind,z_ind{n});
        
        ind=find(~isnan(CloudFractiona));
        if isempty(ind)
            CloudFraction_mean(n)=NaN;
        else
            CloudFraction_mean(n)=nansum(CloudFractiona(ind));  
        end
        clearvars CloudFractiona ind
    end
end
clearvars n

if isempty(lat_ind)
    num_obs=0;
else
    num_obs=numel(lat_ind);
end
