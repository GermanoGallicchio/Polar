function [RLent] = ...
    PolarS_transformation(RLen,transformation)

% transforms data with domain [0, 1] into ...various options
% 
% USAGE:
%
% Rt = PolarS_transformation(R,'logit')
%
% Rt = PolarS_transformation(R,'atanh')
%
%
%
%
% INPUT:        
%
% RLen          resultant length
%               any size
%               must be a real-valued number
%               must have [0, 1] domain
%
% tranformation       char
%                       logit transformation
%

% OUTPUT
%
% RLent         transformed resultant length
%               complex-valued
%               resultant length can be obtained as the modulus of R (i.e., abs(R))
%               resultant angle can be obtained as the argument of R (i.e., angle(R))
% 
%
% Author: Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sort input

p = inputParser;

% add required parameter
addRequired(p, 'RLen')
addRequired(p, 'transformation')


%% sanity checks

% RLen is real
if imag(RLen)~=0
    error('RLen must be real')
end

%% functions

logit = @(x) log((x./(1-x)));

%% output

switch transformation
    case 'logit'
        RLent = logit(RLen);
    case {'atanh' 'FisherZ'}
        RLent = atanh(RLen);
    otherwise
        error('this transformation has not been coded')
end