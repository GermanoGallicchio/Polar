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

%% get data

nSamples = length(theta);

complVec = rho .* exp(1i.*theta);

%% implementation

% length of mean resultant
meanResultantLength = abs(sum(complVec)) / nSamples;

% normalized length of mean resultant
denom = sum(abs(complVec));
if denom~=0
    meanResultantLengthNorm = abs(sum(complVec)) / denom;
else
    meanResultantLengthNorm = 0;
end

% TO DO
% "debiased" (just centered) length of mean resultant

% TO DO
% "debiased" (just centered) normalized length of mean resultant


% angle of mean resultant
meanResultantAngle = angle(sum(complVec));




% U-statistic estimator of the squared length of mean resultant
total_sum = sum(sum(complVec * complVec'));
diag_sum = sum(complVec .* conj(complVec));
UmeanResultantSquaredLength = (total_sum - diag_sum) / (nSamples*(nSamples-1));
UmeanResultantSquaredLength = real(UmeanResultantSquaredLength);    % due to computer rounding precision, there can be a tiny residual imaginary part that should not be there and that must be ignored

% normalized U-statistic estimator of the squared length of mean resultant
total_sum = sum(sum(complVec * complVec'));
diag_sum = sum(complVec .* conj(complVec));
numerator = (total_sum - diag_sum);
total_sum = sum(sum(abs(complVec) * abs(complVec)'));
diag_sum = sum(abs(complVec) .* abs(complVec));
denom = total_sum - diag_sum;
if denom ~= 0
    UmeanResultantSquaredLengthNorm =  numerator / denom;
    UmeanResultantSquaredLengthNorm = real(UmeanResultantSquaredLengthNorm);    % due to computer rounding precision, there can be a tiny residual imaginary part that should not be there and that must be ignored
else
    UmeanResultantSquaredLengthNorm = 0;
end



% --- PPC --- temporary just to check that my code above works

Z = exp(1i * theta);
% Sum across trials: [freq × time × 1 × chan]
sumZ = sum(Z, 1);
% Squared magnitude of the vector sum
sqMag = abs(sumZ).^2;
% Apply PPC formula
numer = sqMag - nSamples;
denom = nSamples * (nSamples - 1);
% Result: [freq × time × 1 × chan] → squeeze to [freq×time×chan]
PPC = squeeze(numer ./ denom);
%----

disp('values are:')
disp(meanResultantLength)
disp(meanResultantLengthNorm)
disp(meanResultantAngle)
disp(UmeanResultantSquaredLength)
disp(UmeanResultantSquaredLengthNorm)
disp(PPC)

%% sort output

metrics.meanResultantLength             = meanResultantLength;
metrics.meanResultantLengthNorm         = meanResultantLengthNorm;
metrics.meanResultantAngle              = meanResultantAngle;
metrics.UmeanResultantSquaredLength     = UmeanResultantSquaredLength;
metrics.UmeanResultantSquaredLengthNorm = UmeanResultantSquaredLengthNorm;





