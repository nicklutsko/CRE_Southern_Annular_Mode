% function to categorize clouds scenes based on cloud vertical structure
% Inputs:
% (1) CloudLayers - number of cloud layers in scene (num_granules)
% (2) LayerTop/LayerBase - altitude of cloud layer top and base (m)
%
% Outputs
% (1) cloud_labels - flag to indicate if high/mid/low cloud is present.
% first entry is H, second entry is M, thrid entry is L

function [cloud_labels] = cloud_categorization(CloudLayers,LayerTop,LayerBase)

low_threshold=3000; % height threshold for difference between L/M categories (m)
high_threshold=7000; % height threshold for difference between M/H categories (m)

% set entries with no cloud layers to NaN
for n=1:numel(CloudLayers)
    if CloudLayers(n)<5
        LayerTop(n,(CloudLayers(n)+1):end)=NaN;
        LayerBase(n,(CloudLayers(n)+1):end)=NaN;
    end
end

cloud_labels=zeros(numel(CloudLayers),3);

for n=1:numel(CloudLayers)
    if CloudLayers(n)>0
        % check for high cloud
        if any(LayerTop(n,:)>=high_threshold)
            cloud_labels(n,1)=1;
        end
        
        % check for low cloud
        if any(LayerBase(n,:)<low_threshold)
            cloud_labels(n,3)=1;
        end
        
        % check for mid-level cloud
        if any(LayerTop(n,:)>=low_threshold & LayerTop(n,:)<high_threshold) | ...
                any(LayerBase(n,:)<high_threshold & LayerBase(n,:)>=low_threshold) | ...
                any(LayerTop(n,:)>=high_threshold & LayerBase(n,:)<low_threshold)
            cloud_labels(n,2)=1;
        end
    end
end
        
end
