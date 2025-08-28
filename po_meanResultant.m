function metrics = po_meanResultant(theta, rho)

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
%% get data

nSamples = length(theta);

% assemble theta and rho to make complex numbers
complVec = rho .* exp(1i.*theta);

% matrixify complex numbers (so that all computations will happen across rows)
nDims = size(complVec);
nRows = nDims(1);
nCols = prod(nDims(2:end));
complVec_matrix = reshape(complVec, nRows,nCols);

%% implementation

% length of mean resultant
metric = nan(1,nCols);
for colIdx = 1:nCols
    c = complVec_matrix(:,colIdx);
    metric(1,colIdx) = abs(sum(c,1)) / nSamples;
end
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
meanResultantLength = metric_reshaped;


% normalized length of mean resultant
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
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
meanResultantLengthNorm = metric_reshaped;


% TO DO
% "debiased" (just centered) length of mean resultant

% TO DO
% "debiased" (just centered) normalized length of mean resultant


% angle of mean resultant
metric = nan(1,nCols);
for colIdx = 1:nCols
    c = complVec_matrix(:,colIdx);
    metric(1,colIdx) = angle(sum(c,1));
end
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
meanResultantAngle = metric_reshaped;



% U-statistic estimator of the squared length of mean resultant
metric = nan(1,nCols);
for colIdx = 1:nCols
    c = complVec_matrix(:,colIdx);
    total_sum = sum(sum(c*c'));  % sum all elements of the matrix
    diag_sum  = sum(diag(c*c')); % sum only the diagonal
    metric(1,colIdx) = (total_sum - diag_sum) / (nSamples*(nSamples-1));
end
metric = real(metric);    % due to computer rounding precision, there can be a tiny residual imaginary part that should not be there and that must be ignored
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
UmeanResultantSquaredLength = metric_reshaped;


% normalized U-statistic estimator of the squared length of mean resultant
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
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
UmeanResultantSquaredLengthNorm = metric_reshaped;

% --- temporary just to check that strategy applied outside this function works
% dPAC
% build unit-length vectors with phase information only
metric = nan(1,nCols);
for colIdx = 1:nCols
    
    c = complVec_matrix(:,colIdx);
    %theta2 = angle(c);
    
    %metric(1,colIdx) = abs(mean(  rho .* exp(1i*theta2)    )) ;

    
    
    phi = mean(1.*exp(1i*theta));
    metric(1,colIdx) = abs(mean(  rho .* (exp(1i*theta)-phi)    )) ;
    
end
metric_reshaped = reshape(metric,1,nDims(2:end));  % reshape into original size, except the row dimension is squeezed
metrics.dPAC = metric_reshaped;
% ---


% --- PPC --- temporary just to check that my code above works
Z = exp(1i * theta);
% Sum across trials: [freq × time × 1 × chan]
sumZ = sum(Z, 1);
% Squared magnitude of the vector sum
sqMag = abs(sumZ).^2;
% Apply PPC formula
numer = sqMag - nSamples;
denominator = nSamples * (nSamples - 1);
% Result: [freq × time × 1 × chan] → squeeze to [freq×time×chan]
PPC = squeeze(numer ./ denominator);
%----



%% sort output

metrics.meanResultantLength             = meanResultantLength;
metrics.meanResultantLengthNorm         = meanResultantLengthNorm;
metrics.meanResultantAngle              = meanResultantAngle;
metrics.UmeanResultantSquaredLength     = UmeanResultantSquaredLength;
metrics.UmeanResultantSquaredLengthNorm = UmeanResultantSquaredLengthNorm;


disp('values are:')
disp(metrics)
% disp(meanResultantLength)
% disp(meanResultantLengthNorm)
% disp(meanResultantAngle)
% disp(UmeanResultantSquaredLength)
% disp(UmeanResultantSquaredLengthNorm)
% disp(PPC)



