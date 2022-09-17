% function to calculate overcast heating rates and ACRE
% Inputs:
% (1) ind - index over which averaging is performed (can be empty, and if so NaN is returned)
% (2) data_in - data to be averaged
% (3) dp - matrix for pressure weighting of average
%
% Outputs:
% data_out - pressure-weighted average of "data_in" over the indices "ind"


function [data_out]=compute_mean_heating(ind,data_in,dp)

if isempty(ind)
    data_out=NaN;
else
    data_out=nansum(data_in(ind).*dp(ind))./nansum(dp(ind));
end
    
    
