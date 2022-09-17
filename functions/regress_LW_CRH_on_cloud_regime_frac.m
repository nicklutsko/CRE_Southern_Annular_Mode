% function to compute regression of LW cloud radiative heating on cloud-category
% frequency of occurrence. These values are used to decompose the overall
% cloud radaitive heating anomaly associated with SAM anomaly into
% contributions from different cloud regimes
%
% Inputs:
% (1) lat, time, z: vetors of latitude, time, and height
% (2) L_frac, M_frac, etc.: fractional occurrence of cloud regimes (units: %,
% dimensions: lat, time)
% (3) LW_CRH: LW cloud radiative heating (units: K d-1, dimensions: lat,
% z, time)
% (4) num_obs: number of observations (dimensions: lat, time)
%
% Outputs:
% (1) dCRH_dL, dCRH_dM, etc.: regression coefficients for LW_CRH regressed
% on L_frac, M_frac, etc. with multiple linear regression model (i.e. represents partial derivatives)
% (units: K d-1 %-1; dimensions: lat, z)
% (2) dCRH_dL_CI, dCRH_dM_CI, etc.: 95% confidence intervals for regression
% coefficients (units: K d-1 %-1; dimensions: lat,z)

function [dCRH_dL,dCRH_dM,dCRH_dH,dCRH_dLM,dCRH_dLH,dCRH_dMH,dCRH_dLMH,...
    dCRH_dL_CI,dCRH_dM_CI,dCRH_dH_CI,dCRH_dLM_CI,dCRH_dLH_CI,dCRH_dMH_CI,dCRH_dLMH_CI] = ...
    regress_LW_CRH_on_cloud_regime_frac(lat,time,z,L_frac,M_frac,H_frac,LM_frac,...
    LH_frac,MH_frac,LMH_frac,LW_CRH,num_obs)

% remove times with missing data
time_ind=[];
for n=1:numel(time)
    if any(num_obs(:,n)~=0)
        time_ind(end+1)=n;
    end
end
clearvars n

LW_CRH=LW_CRH(:,:,time_ind);
L_frac=L_frac(:,time_ind);
M_frac=M_frac(:,time_ind);
H_frac=H_frac(:,time_ind);
LM_frac=LM_frac(:,time_ind);
LH_frac=LH_frac(:,time_ind);
MH_frac=MH_frac(:,time_ind);
LMH_frac=LMH_frac(:,time_ind);
num_obs=num_obs(:,time_ind);
time=time(time_ind);
clearvars time_ind

% remove mean
LW_CRH_save=NaN*ones(size(LW_CRH));
L_frac_save=NaN*ones(size(L_frac));
M_frac_save=NaN*ones(size(M_frac));
H_frac_save=NaN*ones(size(H_frac));
LM_frac_save=NaN*ones(size(LM_frac));
LH_frac_save=NaN*ones(size(LH_frac));
MH_frac_save=NaN*ones(size(MH_frac));
LMH_frac_save=NaN*ones(size(L_frac));

for l=1:numel(time)
    LW_CRH_save(:,:,l)=LW_CRH(:,:,l)-nanmean(LW_CRH,3);
    L_frac_save(:,l)=L_frac(:,l)-nanmean(L_frac,2);
    M_frac_save(:,l)=M_frac(:,l)-nanmean(M_frac,2);
    H_frac_save(:,l)=H_frac(:,l)-nanmean(H_frac,2);
    LM_frac_save(:,l)=LM_frac(:,l)-nanmean(LM_frac,2);
    LH_frac_save(:,l)=LH_frac(:,l)-nanmean(LH_frac,2);
    MH_frac_save(:,l)=MH_frac(:,l)-nanmean(MH_frac,2);
    LMH_frac_save(:,l)=LMH_frac(:,l)-nanmean(LMH_frac,2);
end
LW_CRH=LW_CRH_save;
L_frac=L_frac_save;
M_frac=M_frac_save;
H_frac=H_frac_save;
LM_frac=LM_frac_save;
LH_frac=LH_frac_save;
MH_frac=MH_frac_save;
LMH_frac=LMH_frac_save;
clearvars LW_CRH_save l L_frac_save M_frac_save H_frac_save LM_frac_save LH_frac_save MH_frac_save LMH_frac_save


% regress LW_CRH on L_frac, M_frac, etc. (coefficient is proportional to mean
% overcast cloud radiative heating rate for cloud regime)   
dCRH_dL=NaN*ones(numel(lat),numel(z));
dCRH_dM=NaN*ones(numel(lat),numel(z));
dCRH_dH=NaN*ones(numel(lat),numel(z));
dCRH_dLM=NaN*ones(numel(lat),numel(z));
dCRH_dLH=NaN*ones(numel(lat),numel(z));
dCRH_dMH=NaN*ones(numel(lat),numel(z));
dCRH_dLMH=NaN*ones(numel(lat),numel(z));

dCRH_dL_CI=NaN*ones(numel(lat),numel(z));
dCRH_dM_CI=NaN*ones(numel(lat),numel(z));
dCRH_dH_CI=NaN*ones(numel(lat),numel(z));
dCRH_dLM_CI=NaN*ones(numel(lat),numel(z));
dCRH_dLH_CI=NaN*ones(numel(lat),numel(z));
dCRH_dMH_CI=NaN*ones(numel(lat),numel(z));
dCRH_dLMH_CI=NaN*ones(numel(lat),numel(z));

for n=1:numel(lat)
    for m=1:numel(z)
        LW_CRHa=squeeze(LW_CRH(n,m,:));
        L_fraca=squeeze(L_frac(n,:));
        M_fraca=squeeze(M_frac(n,:));
        H_fraca=squeeze(H_frac(n,:));
        LM_fraca=squeeze(LM_frac(n,:));
        LH_fraca=squeeze(LH_frac(n,:));
        MH_fraca=squeeze(MH_frac(n,:));
        LMH_fraca=squeeze(LMH_frac(n,:));
        num_obsa=squeeze(num_obs(n,:));
        
        ind=find(num_obsa>0);
        LW_CRHa=LW_CRHa(ind); 
        L_fraca=transpose(L_fraca(ind));
        M_fraca=transpose(M_fraca(ind));
        H_fraca=transpose(H_fraca(ind));
        LM_fraca=transpose(LM_fraca(ind));
        LH_fraca=transpose(LH_fraca(ind));
        MH_fraca=transpose(MH_fraca(ind));
        LMH_fraca=transpose(LMH_fraca(ind));
        clearvars ind num_obsa
        
        X=NaN*ones(numel(time),7);
        X(:,1)=L_fraca;
        X(:,2)=M_fraca;
        X(:,3)=H_fraca;
        X(:,4)=LM_fraca;
        X(:,5)=LH_fraca;
        X(:,6)=MH_fraca;
        X(:,7)=LMH_fraca;
        
        mdl=fitlm(X,LW_CRHa);
        dCRH_dL(n,m)=mdl.Coefficients.Estimate(2);
        dCRH_dL_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,L_fraca,mdl.Coefficients.SE(2));
        
        dCRH_dM(n,m)=mdl.Coefficients.Estimate(3);
        dCRH_dM_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,M_fraca,mdl.Coefficients.SE(3));
        
        dCRH_dH(n,m)=mdl.Coefficients.Estimate(4);
        dCRH_dH_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,H_fraca,mdl.Coefficients.SE(4));
        
        dCRH_dLM(n,m)=mdl.Coefficients.Estimate(5);
        dCRH_dLM_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,LM_fraca,mdl.Coefficients.SE(5));
        
        dCRH_dLH(n,m)=mdl.Coefficients.Estimate(6);
        dCRH_dLH_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,LH_fraca,mdl.Coefficients.SE(6));
        
        dCRH_dMH(n,m)=mdl.Coefficients.Estimate(7);
        dCRH_dMH_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,MH_fraca,mdl.Coefficients.SE(7));
        
        dCRH_dLMH(n,m)=mdl.Coefficients.Estimate(8);
        dCRH_dLMH_CI(n,m)=compute_regression_slope_95_CI(LW_CRHa,LMH_fraca,mdl.Coefficients.SE(8));       
        clearvars mdl X L_fraca M_fraca H_fraca LM_fraca LH_fraca MH_fraca LMH_fraca LW_CRHa              
    end
end
clearvars n m    





