function [theta_centered, rho_centered ] = po_recenter(theta, rho)

% DESCRIPTION: Center to the phase-only mean along first dimension
%
% INPUT:        
%
%
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%%


%% implementation

data_mean = mean(exp(1i*theta),1);  % phase clustering over dimension 1


% center data
data_centered = rho.*(exp(1i*theta)-data_mean); % automatic broadcasting in modern matlab


rho_centered   = abs(data_centered);
theta_centered = angle(data_centered);

