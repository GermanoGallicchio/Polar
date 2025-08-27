function [theta_centered, rho_centered ] = po_recenter(theta, rho)

% DESCRIPTION: Center to the phase-only mean along first dimension
%
% INPUT:        
%
%
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% implementation

data = rho.*exp(1i*theta);

data_mean = mean(data,1);

data_centered = nan(size(data));
for rIdx = 1:size(data,1)
    data_centered(rIdx,:) = data(rIdx,:)-data_mean(1,:);
end

rho_centered   = abs(data_centered);
theta_centered = angle(data_centered);

