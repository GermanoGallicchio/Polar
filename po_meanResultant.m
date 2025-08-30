function metrics = po_meanResultant(po_cfg, theta, rho)

% Computes mean resultant metrics
% 
% SYNTAX:
%
% DESCRIPTION:
%   This function calculates the mean resultant vector from a set of phase
%   angles (and optionally, magnitudes) in circular data. The resultant
%   vector represents the mean direction and magnitude of the data points.
%
%
% INPUT:        
%
%   theta - A column vector (N-by-1) of phase angles in radians.
%               
%
% OPTIONAL INPUT:
%   rho   - A column vector (N-by-1) of magnitudes (any unit). 
%                  Must be non-negative. If no input, it assumes as default
%                  a vector of ones (equal weights for all angles).
%
% 
% EXAMPLES:
%
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% sort input

%% sanity checks

% theta and rho are same size
if size(theta)~=size(rho)
    error('theta and rho must be of the same size')
end

% because these metrics work on the rows (1st dim) a row vector is likely a user mistake
if length(size(theta))==2
    if size(theta,1)==1  &  size(theta,2)>1
        warning('metrics are computed on the first dimension. you might want to transpose theta and rho');
    end
end


% requests are possible
metricsList = [
    "meanResultantLength" 
    "meanResultantLengthNorm" 
    "meanResultantAngle" 
    "UmeanResultantSquaredLength" 
    "UmeanResultantSquaredLengthNorm"
    ];

metricsRequested = metricsList(ismember(metricsList,po_cfg.requests));
if isempty(metricsRequested)
    disp(metricsList)
    error('po_cfg.requests should be a string vector of at least one of the above')
end

%% get data

nSamples = size(theta,1);

% assemble theta and rho to make complex numbers
complVec = rho .* exp(1i.*theta);

% matrixify complex numbers (so that all computations will happen across rows)
nDims = size(complVec);
nRows = nDims(1);
nCols = prod(nDims(2:end));
complVec_matrix = reshape(complVec, nRows,nCols);

%% implementation

% meanResultantLength
if any(strcmp('meanResultantLength',metricsRequested))
    metric = nan(1,nCols);
    for colIdx = 1:nCols
        c = complVec_matrix(:,colIdx);
        metric(1,colIdx) = abs(sum(c,1)) / nSamples;
    end
    metric_reshaped = reshape(metric,[ 1 nDims(2:end)]);  % reshape into original size, except the row dimension is squeezed
    meanResultantLength = metric_reshaped;

    % sort output
    metrics.meanResultantLength             = meanResultantLength;
end


% meanResultantLengthNorm
if any(strcmp('meanResultantLengthNorm',metricsRequested))
    metric = nan(1,nCols);
    for colIdx = 1:nCols
        c = complVec_matrix(:,colIdx);
        numerator = abs(sum(c,1));
        denominator = sum(abs(c),1);
        if denominator~=0
            metric(1,colIdx) = numerator ./ denominator;
        else
            metric(1,colIdx) = 0;
        end
    end
    metric_reshaped = reshape(metric,[ 1 nDims(2:end)]);  % reshape into original size, except the row dimension is squeezed
    meanResultantLengthNorm = metric_reshaped;

    % sort output
    metrics.meanResultantLengthNorm         = meanResultantLengthNorm;
end


% meanResultantAngle
if any(strcmp('meanResultantAngle',metricsRequested))
    metric = nan(1,nCols);
    for colIdx = 1:nCols
        c = complVec_matrix(:,colIdx);
        metric(1,colIdx) = angle(sum(c,1));
    end
    metric_reshaped = reshape(metric,[ 1 nDims(2:end)]);  % reshape into original size, except the row dimension is squeezed
    meanResultantAngle = metric_reshaped;

    % sort output
    metrics.meanResultantAngle              = meanResultantAngle;
end


% UmeanResultantSquaredLength 
% U-statistic estimator of the squared length of mean resultant
if any(strcmp('UmeanResultantSquaredLength',metricsRequested))
    metric = nan(1,nCols);
    for colIdx = 1:nCols
        c = complVec_matrix(:,colIdx);
        total_sum = sum(sum(c*c'));  % sum all elements of the matrix
        diag_sum  = sum(diag(c*c')); % sum only the diagonal
        metric(1,colIdx) = (total_sum - diag_sum) / (nSamples*(nSamples-1));
    end
    metric = real(metric);    % due to computer rounding precision, there can be a tiny residual imaginary part that should not be there and that must be ignored
    metric_reshaped = reshape(metric,[ 1 nDims(2:end)]);  % reshape into original size, except the row dimension is squeezed
    UmeanResultantSquaredLength = metric_reshaped;
    
    % sort output
    metrics.UmeanResultantSquaredLength     = UmeanResultantSquaredLength;
end

% UmeanResultantSquaredLengthNorm
% normalized U-statistic estimator of the squared length of mean resultant
if any(strcmp('UmeanResultantSquaredLengthNorm',metricsRequested))
    metric = nan(1,nCols);
    for colIdx = 1:nCols
        c = complVec_matrix(:,colIdx);

        total_sum = sum(sum(c*c'));  % sum all elements of the matrix
        diag_sum  = sum(diag(c*c')); % sum only the diagonal
        numerator = (total_sum - diag_sum);

        total_sum = sum(sum(abs(c) * abs(c)'));
        diag_sum = sum(diag(abs(c)*abs(c)'));
        denominator = total_sum - diag_sum;

        if denominator ~= 0
            metric(1,colIdx) =  numerator ./ denominator;
        else
            metric(1,colIdx) = 0;
        end
    end
    metric = real(metric);    % due to computer rounding precision, there can be a tiny residual imaginary part that should not be there and that must be ignored
    metric_reshaped = reshape(metric,[ 1 nDims(2:end)]);  % reshape into original size, except the row dimension is squeezed
    UmeanResultantSquaredLengthNorm = metric_reshaped;

    % sort output
    metrics.UmeanResultantSquaredLengthNorm = UmeanResultantSquaredLengthNorm;
end




%disp('values are:')
%disp(metrics)
% disp(meanResultantLength)
% disp(meanResultantLengthNorm)
% disp(meanResultantAngle)
% disp(UmeanResultantSquaredLength)
% disp(UmeanResultantSquaredLengthNorm)
% disp(PPC)



