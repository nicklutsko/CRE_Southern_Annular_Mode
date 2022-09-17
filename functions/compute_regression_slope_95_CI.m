% function to determine determine 95% confidence interval for regression 
% coefficient dy/dx  The function is intended to be used with y representing cloud fraction or
% cloud radaitive heating and x representing the standardized SAM index.
%
% Inputs:
% (1) y: vector cloud cloud radiative heating or cloud fraction (1, time)
% (2) x: vector of standardized SAM index (1, time)
% (3) dy_dx_SE: standard error of regression slope dy_dx (scalar)
%
% Outputs:
% (1) dy_dx_CI: 95% confidence interval for dy_dx regression slope

function dy_dx_CI = compute_regression_slope_95_CI(y,x,dy_dx_SE)

% compute DOFs in the time dimension
N_nom=numel(y);

r1_x=corrcoef(x(1:(end-1)),x(2:end));
r1_x=r1_x(1,2);
r1_y=corrcoef(y(1:(end-1)),y(2:end));
r1_y=r1_y(1,2);
N_eff=N_nom*(1-r1_x*r1_y)/(1+r1_x*r1_y); % equation 31 of Bretherton et al. 1999
clearvars r1_x r1_y

% compute critical value from students t-distribution
t=tinv(0.975,N_eff-2);

% compute confidence interval
dy_dx_CI=t*dy_dx_SE*sqrt(N_nom/N_eff);








