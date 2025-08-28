function [theta_centered, rho_centered ] = po_demean(theta, rho)

% DESCRIPTION: Center to the phase-only mean along first dimension
%
% INPUT:        
%
%
% usage: feed all 1s as rho for dPAC, then element wise product of magnitude and complex vector created by combining rho_centered and theta_centered
% 
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%%


%% implementation

compl_mean = mean(rho.*exp(1i*theta),1); 

% phase clustering over dimension 1
centered_complex = rho.*exp(1i*theta) - compl_mean;   % complex difference vectors


theta_centered = angle(centered_complex);
rho_centered = abs(centered_complex);



