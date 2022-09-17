% function to compute mean frequency of occurrence for cloud categories
% Inputs:
% (1) freq - frequency of occurrence of a cloud regime (lat,time)
% (2) lat - latitude vector
% (3) num_obs - number of observations (lat,time)
% 
% Outputs:
% (1) mean_freq - mean frequency of occurrence weighted by area and number
% of observations

function freq_mean = compute_mean_freq_of_occurrence(freq,lat,num_obs)

% average in time
freq_time_mean=NaN*ones(size(lat));

for n=1:numel(lat)
    num_obsa=num_obs(n,:);
    Wa=num_obsa./max(num_obsa);
    freqa=freq(n,:);

    freq_time_mean(n)=nansum(freqa(:).*Wa(:))./nansum(Wa(:));
end
clearvars n freqa Wa num_obsa

% avearge over latitude
W=cos(lat.*pi./180);
freq_mean=nansum(freq_time_mean(:).*W(:))./nansum(W(:));


