% function to compute LW CRH climatology and cloud fraction climatology for
% different cloud regimes (L, LM, LMH, etc.)
%
% Inputs:
% (1) LW_CRH_regime: LW CRH for cloud regime (units: K d-1; dimension: lat, z,
% time)
% (2) f_regime: fractional occurrence of cloud regime (units: 0-1;
% dimensions: lat,time)
% (3) num_obs_regime: number of obervations of regime (dimensions: lat,
% time)
% (4) lat, z: vectors of latitude and height
%
% Outputs:
% (1) LW_CRH_mean: regime-avearge LW CRH (units: K d-1)
% (2) f_mean: mean fractional occurrence of cloud regime (units: 0-1)

function [LW_CRH_mean,f_mean]=regime_average_CRH_cloud_frac(LW_CRH_regime,f_regime,num_obs_regime,lat,z)

% compute LW CRH climatology
LW_CRH_mean=NaN.*ones(numel(lat),numel(z));

for n=1:numel(lat)
    for m=1:numel(z)
        LW_CRHa=squeeze(LW_CRH_regime(n,m,:));
        Wa=squeeze(num_obs_regime(n,:));
        Wa=Wa./max(Wa);
        
        LW_CRH_mean(n,m)=nansum(LW_CRHa(:).*Wa(:))./nansum(Wa(:));
        clearvars LW_CRHa Wa
    end
end
clearvars n m

% compute cloud fraction climatology
f_mean=NaN.*ones(size(lat));

for n=1:numel(lat)
    fa=f_regime(n,:);
    Wa=squeeze(num_obs_regime(n,:));
    Wa=Wa./max(Wa);
    
    f_mean(n)=nansum(fa(:).*Wa(:))./nansum(Wa(:));
    clearvars fa Wa
end
clearvars n









