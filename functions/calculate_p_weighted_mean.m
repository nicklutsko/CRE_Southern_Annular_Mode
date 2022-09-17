% function to calculate pressure-weighted mean heating rate
% Inputs:
% (1) LW_hr, LW_hr_clr : LW heating rates (dimensions: footprint #, height;
% units: Kd-1)
% (2) dp: pressure difference for layer (dimensions: footprint #, height;
% units: Pa)
% (3) lat_ind: index indntifying footprints with a particular cloud
% category L, M, H, etc. (dimensions: footprint #)
% (4) z_ind: index for height averaging
% (5) z_mdpts: midpoints for height averaging
% 
% Outputs:
% (1) LW_CRH_mean: avearge LW cloud radiative heating for footprints from
% index "lat_ind" and vertical intervals corresponding to indices "z_ind"
% (units: Kd-1; dimensions: z_mdpts)

function [LW_CRH_mean, num_obs] = calculate_p_weighted_mean(LW_hr,LW_hr_clr,dp,lat_ind,z_ind,z_mdpts)

for n=1:numel(z_mdpts)
    if isempty(lat_ind)
        LW_CRH_mean(n)=NaN;
    else
        LW_CRHa=LW_hr(lat_ind,z_ind{n}) - LW_hr_clr(lat_ind,z_ind{n});
        dpa=dp(lat_ind,z_ind{n});
        
        ind=find(~isnan(dpa) & ~isnan(LW_CRHa));
        if isempty(ind)
            LW_CRH_mean(n)=NaN;
        else
            LW_CRH_mean(n)=nansum(LW_CRHa(ind).*dpa(ind))./nansum(dpa(ind));  
        end
        clearvars dpa LW_CRHa ind
    end
end
clearvars n

if isempty(lat_ind)
    num_obs=0;
else
    num_obs=numel(lat_ind);
end




