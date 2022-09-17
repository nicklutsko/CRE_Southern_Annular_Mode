% compute mean LW cloud radaitive heating weighted by number of
% observations
%
% Inputs:
% (1) lat, z - vectors of latitude and height
% (2) LW_CRH - overcast LW cloud radiative heating (dimensions: lat, z,
% time; units: K d-1)
% (3) num_obs - number of observations (dimensions: lat, time)
%
% Outputs:
% (1) LW_CRH_mean - mean LW cloud radative heating weighted by num_obs
% (dimensions: lat, z; units: K d-1)

function [LW_CRH_mean] = compute_number_weighted_mean(lat,z,LW_CRH,num_obs)

for n=1:numel(lat)
    for m=1:numel(z)
        LW_CRHa=squeeze(LW_CRH(n,m,:));
        num_obsa=squeeze(num_obs(n,:));
        
        ind=find(num_obsa(:)>0 & ~isnan(LW_CRHa(:)));
        if isempty(ind)
            LW_CRH_mean(n,m)=NaN;
        else
            LW_CRHa=LW_CRHa(ind);
            num_obsa=num_obsa(ind);           
            num_obsa=num_obsa./max(num_obsa);           
            LW_CRH_mean(n,m)=nansum(LW_CRHa(:).*num_obsa(:))./nansum(num_obsa(:));
        end
        clearvars ind LW_CRHa num_obsa
    end
end
clearvars n m
            