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

theta_mean = mean(exp(1i*theta),1);  % phase clustering over dimension 1
centered_complex = exp(1i*theta) - theta_mean;   % complex difference vectors


theta_centered = angle(centered_complex);
rho_centered = abs(centered_complex);



