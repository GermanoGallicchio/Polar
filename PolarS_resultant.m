function [MR, idx_toTrim] = ...
    PolarS_resultant(phaseVec, varargin)

% computes resultant vector 
% 
% USAGE:
%
% MR = PS_meanResultant(phaseVec)
%
% MR = PS_meanResultant(phaseVec, magnitudeVec)
%
% MR = PS_meanResultant(phaseVec, magnitudeVec, 'trim', 20)
%
%
%
% INPUT:        
%
% phaseVec      phase angles [radians]
%               size: N, 1 (column vector) 
%               
%
% Optional INPUT:
%
% magnitudeVec  magnitudes [any unit]
%               must be non negative
%               size: N, 1 (column vector) 
%
%
% Optional 'Name',Value pairs
%
% 'Trim'        a number expressing a percentage to trim
%

% OUTPUT
%
% MR            mean resultant
%               complex-valued
%               mean resultant length can be obtained as the modulus of R (i.e., abs(R))
%               mean resultant angle can be obtained as the argument of R (i.e., angle(R))
% 
%
% Author: Germano Gallicchio 
% germano.gallicchio@gmail.com

%% sort input

p = inputParser;

% add required parameter
addRequired(p, 'phaseVec')

% add optional parameters with default values
addParameter(p, 'magnitudeVec', ones(size(phaseVec)))
addParameter(p, 'trim', 0);

% parse the input arguments
parse(p, phaseVec, varargin{:});

% get optional parameters into variables
magnitudeVec = p.Results.magnitudeVec;
trim = p.Results.trim;

%% sanity checks

% phaseVec is column vector
if size(phaseVec,2) ~= 1
    error('phaseVec should have size N-by-1, that is a column vector')
end

%% get data

N = length(phaseVec);

%% trimming

if trim>0
    
    % compute angular distance matrix (two at a time: vector 'a' and in turn all other vectors)
    D = nan(N,N);
    for vIdx_a = 1:N
        for vIdx_b = 1:N
            D(vIdx_a,vIdx_b) = min([abs(phaseVec(vIdx_b)-phaseVec(vIdx_a)) 2*pi-abs(phaseVec(vIdx_b)-phaseVec(vIdx_a))]);
        end
    end

    % sum of angular distances of each vector from all other vectors 
    d = sum(D,2);

    % compute weight from magnitudeVec 
    %
    % the weight is essentially magnitudeVec reversed so that small values 
    % become large and viceversa. the boundaries of small and large are 
    % defined as the min and max of the sum of distances.
    if length(unique(magnitudeVec))~=1
        w = normalize(-abs(magnitudeVec), 'range', [min(d) max(d)]);
    else
        % if magnitudeVec is all 1s, then weight plays no role
        w = median(d);
    end

    
    % weighted sum of angular distances
    % large values are achieved as combination of 
    % - vector being far from the bulk of the other vectors => large sum of distances
    % - vector having a lower magnitude => small magnitude
    % vectors that are far from the bulk are likely to be removed unless their large magnitude counters the large sum of distances
    % vectors that are in the bulk are likely to be retained unless their small magnitude counters the small sum of distances

    wd = w.*d;
    
    % trim percentage (defined by 'trim') of outliers
    % find the largest p% and remove them
    idx_toTrim = wd > prctile(wd,100-trim);
    phaseVec(idx_toTrim) = [];
    magnitudeVec(idx_toTrim) = [];

    
end



%% implementation

MR = mean(magnitudeVec.* exp(1i.*phaseVec)) / mean(magnitudeVec);



