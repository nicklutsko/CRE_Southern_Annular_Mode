% function to compute dCRH_dSAM from decomposition into contributions from
% different cloud regimes
%
% Inputs:
% (1) lat, time, z: vetors of latitude, time, and height
% (2) L_frac, M_frac, etc.: fractional occurrence of cloud regimes (units: %,
% dimensions: lat, time)
% (3) LW_CRH: LW cloud radiative heating (units: K d-1, dimensions: lat,
% z, time)
% (4) num_obs: number of observations (dimensions: lat, time)
% (5) SAM: SAM index (dimensions: time)
%
% Outputs:
% (1) dCRH_dSAM: change in CRH with respect to SAM estimated from
% decomposition by cloud regimes
% (units: K d-1 sigma_SAM-1; dimensions: lat, z)
% (2) dCRH_dSAM_CI: 95% confidence intervals for dCRH_dSAM
% (units: K d-1 sigma_SAM-1; dimensions: lat,z)

function [dCRH_dSAM,dCRH_dSAM_CI] = ...
    dCRH_dSAM_from_decomposition(lat,time,z,L_frac,M_frac,H_frac,LM_frac,...
    LH_frac,MH_frac,LMH_frac,LW_CRH,num_obs,SAM)

% initialize arrays
dCRH_dSAM=NaN*ones(numel(lat),numel(z));
dCRH_dSAM_CI=NaN*ones(numel(lat),numel(z));

%%%%%%%% regress cloud regime fraction on SAM index %%%%%%%%%
[dL_frac_dSAM, dL_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,L_frac,SAM,num_obs);
[dM_frac_dSAM, dM_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,M_frac,SAM,num_obs);
[dH_frac_dSAM, dH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,H_frac,SAM,num_obs);
[dLM_frac_dSAM, dLM_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LM_frac,SAM,num_obs);
[dLH_frac_dSAM, dLH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LH_frac,SAM,num_obs);
[dMH_frac_dSAM, dMH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,MH_frac,SAM,num_obs);
[dLMH_frac_dSAM, dLMH_frac_dSAM_CI] = regress_cloud_regime_frac_on_SAM(lat,time,LMH_frac,SAM,num_obs);

%%%%%%%%% regress CRH on cloud regime fraction %%%%%%%%%
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
        
        dL_frac_dSAMa=dL_frac_dSAM(n);
        dL_frac_dSAM_CIa=dL_frac_dSAM_CI(n);
        dM_frac_dSAMa=dM_frac_dSAM(n);
        dM_frac_dSAM_CIa=dM_frac_dSAM_CI(n);
        dH_frac_dSAMa=dH_frac_dSAM(n);
        dH_frac_dSAM_CIa=dH_frac_dSAM_CI(n);
        dLM_frac_dSAMa=dLM_frac_dSAM(n);
        dLM_frac_dSAM_CIa=dLM_frac_dSAM_CI(n);
        dLH_frac_dSAMa=dLH_frac_dSAM(n);
        dLH_frac_dSAM_CIa=dLH_frac_dSAM_CI(n);
        dMH_frac_dSAMa=dMH_frac_dSAM(n);
        dMH_frac_dSAM_CIa=dMH_frac_dSAM_CI(n);
        dLMH_frac_dSAMa=dLMH_frac_dSAM(n);
        dLMH_frac_dSAM_CIa=dLMH_frac_dSAM_CI(n);
        
        X=NaN*ones(numel(time),7);
        X(:,1)=L_fraca;
        X(:,2)=M_fraca;
        X(:,3)=H_fraca;
        X(:,4)=LM_fraca;
        X(:,5)=LH_fraca;
        X(:,6)=MH_fraca;
        X(:,7)=LMH_fraca;
        
        mdl=fitlm(X,LW_CRHa);
        dCRH_dLa=mdl.Coefficients.Estimate(2);        
        dCRH_dMa=mdl.Coefficients.Estimate(3);       
        dCRH_dHa=mdl.Coefficients.Estimate(4);        
        dCRH_dLMa=mdl.Coefficients.Estimate(5);       
        dCRH_dLHa=mdl.Coefficients.Estimate(6);       
        dCRH_dMHa=mdl.Coefficients.Estimate(7);       
        dCRH_dLMHa=mdl.Coefficients.Estimate(8);
        clearvars X L_fraca M_fraca H_fraca LM_fraca LH_fraca MH_fraca LMH_fraca 
        
        % uncertainty from dCRH_df
        X_pred=[dL_frac_dSAMa dM_frac_dSAMa dH_frac_dSAMa dLM_frac_dSAMa dLH_frac_dSAMa dMH_frac_dSAMa dLMH_frac_dSAMa];
        [y_pred,y_CI]=predict(mdl,X_pred,'alpha',0.05);
        y_CI=(max(y_CI)-min(y_CI))./2;
        N_nom=numel(LW_CRHa(~isnan(LW_CRHa)));
        r1=corrcoef(LW_CRHa(1:(end-1)),LW_CRHa(2:end));
        r1=r1(1,2);
        N_eff=N_nom*(1-r1)/(1+r1);        
        dCRH_dSAM_CI1=y_CI*tinv(0.975,N_eff)/tinv(0.975,N_nom).*sqrt(N_nom/N_eff);
        clearvars X_pred y_CI N_nom N_eff r1
        
        % uncertinaty from dL_frac_dSAM, etc.
        dCRH_dSAM_CI2=sqrt((dCRH_dLa.*dL_frac_dSAM_CIa).^2 + (dCRH_dMa.*dM_frac_dSAM_CIa).^2 + ...
            (dCRH_dHa.*dH_frac_dSAM_CIa).^2 + (dCRH_dLMa.*dLM_frac_dSAM_CIa).^2 + ...
            (dCRH_dLHa.*dLH_frac_dSAM_CIa).^2 + (dCRH_dMHa.*dMH_frac_dSAM_CIa).^2 + ...
            (dCRH_dLMHa.*dLMH_frac_dSAM_CIa).^2);
        
        % combine uncertinaties in quadrature
        dCRH_dSAM_CI(n,m)=sqrt(dCRH_dSAM_CI1.^2 + dCRH_dSAM_CI2.^2);
        dCRH_dSAM(n,m)=y_pred;
        clearvars y_pred dCRH_dSAM_CI1 dCRH_dSAM_CI2 dCRH_dLa dCRH_dMa dCRH_dHa dCRH_dLMa dCRH_dLHa ...
            dCRH_dMHa dCRH_dLMHa dL_frac_dSAM_CIa dM_frac_dSAM_CIa dH_frac_dSAM_CIa dLM_frac_dSAM_CIa ...
            dLH_frac_dSAM_CIa dMH_frac_dSAM_CIa dLMH_frac_dSAM_CIa
    end
end
clearvars n m    

