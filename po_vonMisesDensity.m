function [theta, r] = po_vonMisesDensity(kappa,mu,theta)

% DESCRIPTION: von Mises density determined for given theta values
%
%
% INPUT:
%       kappa   - concentration parameter kappa >= 0
%       mu      - mean direction
%       theta   - vector of angles (radians) in [-pi, pi)
%
% OUTPUT:
%       theta   - n x 1 vector of input angles in [-pi, pi) (passed through)
%       r       - n x 1 vector of density
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% input check

% kappa must be numeric and non-negative
if ~isnumeric(kappa)
    error('kappa must be numeric')
else
    if kappa<0
        error('kappa must be non negative')
    end
end

% mu must be a value in [-pi, pi)
if ~isnumeric(mu)
    error('mu must be numeric')
else
    if mu<-pi | mu>pi
        error('mu must be between -pi and pi')
    end
end


% number of distributions (length of kappa and mu)
if numel(kappa) ~= numel(mu)
    error('kappa and mu must have the same number of elements');
end
nDistributions = numel(mu);



% theta must be not empty and values in [-pi, pi)
if isempty(theta)
    error('theta must be non-empty.');
else
    if any(theta < -pi) || any(theta >= pi)
        error('All theta values must lie in [-pi, pi).');
    end
end
nTheta = numel(theta);


%% implementation

% compute per-distribution density using Bessel function (Berens book p.16, Pewsey book p.56)
vMdensity_perDistribution = nan(nDistributions, nTheta);
for dIdx = 1:nDistributions
    I0 = besseli(0, kappa(dIdx));   % normalization constant
    for tIdx = 1:nTheta
        d = theta(tIdx) - mu(dIdx);
        vMdensity_perDistribution(dIdx, tIdx) = (1 / (2*pi * I0)) * exp(kappa(dIdx) * cos(d));
    end
end

% Average with equal weighting across modes to create a mixture of 
% vonMises distributions. No mixtures if only one mode is provided.
r = mean(vMdensity_perDistribution, 1).';   % column vector

% Basic sanity check
if any(r < 0)
    error('Computed density contains negative values. Check inputs.');
end


end