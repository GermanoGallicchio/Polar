function [theta, r] = po_vonMisesMixSample(kappaVec,muVec,n)

% DESCRIPTION: sample from a mixture of von Mises distributions
%
% INPUT:        
%   kappaVec    - vector of concentration parameters (one per distribution)
%   muVec       - vector of mean directions (one per distribution)
%   n           - number of samples to draw
%
%
% OPTIONAL INPUT:
%   none
%
%
% OUTPUT:
%   theta  -    n x 1 vector angles (in radians)
%   r      -    n x 1 vector of lengths (all 1s)
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% input check

% kappaVec and muVec must have the same size
if size(kappaVec) ~= size(muVec)
    error('kappaVec and muVec must have same size');
end

%% shortcuts

nDistr = length(kappaVec); % number of distributions

%% implementation
% There are multiple samples to draw (randomly) from N distributions.
% Each one-sample draw is taken from one of the possible distributions,
% with equal probability.

theta = zeros(n,1);

for sIdx = 1:n  %  loop through each sample index
    
    dIdx = randi(nDistr);  % choose randomly the distribution index

    % draw one sample from the distribution dIdx
    theta(sIdx) = po_vonMisesSample(kappaVec(dIdx), muVec(dIdx), 1);
end

r = ones(n,1);