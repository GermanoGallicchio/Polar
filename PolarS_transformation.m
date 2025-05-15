function [MRLent] = ...
    PolarS_transformation(MRLen,transformation)

% transforms data with domain [0, 1] into ...various options
% 
% SYNTAX
%
% MRLent = PolarS_transformation(MR,'logit')
%
% MRLent = PolarS_transformation(MR,'atanh')
%
%
%
%
% INPUT:        
%
% MRLen         mean resultant length
%               any size
%               must be a real-valued number
%               must have [0, 1] domain
%
% tranformation       char
%                       'logit', 'atanh', 'FisherZ'
%

% OUTPUT
%
% MRLent         transformed mean resultant length
%               real-valued
% 
%
% Author: Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sort input

p = inputParser;

% add required parameter
addRequired(p, 'MRLen')
addRequired(p, 'transformation')


%% sanity checks

% RLen is real
if imag(MRLen)~=0
    error('RLen must be real')
end

%% functions

logit = @(x) log((x./(1-x)));

%% output

switch transformation
    case 'logit'
        MRLent = logit(MRLen);
    case {'atanh' 'FisherZ'}
        MRLent = atanh(MRLen);
    otherwise
        error('this transformation has not been coded')
end