function [theta, r] = po_vonMisesSample(kappa,mu,n)

% DESCRIPTION: sample from a von Mises distribution on [-pi, pi)
%              Handles both single and mixture distributions
%
%       [theta, r] = po_vonMisesSample(kappa,mu,n)
%
% INPUT:
%       kappa   - concentration parameter kappa >= 0 (scalar or vector for mixture)
%       mu      - mean direction (scalar or vector for mixture)
%       n       - number of samples
%
% OUTPUT:
%       theta   - n x 1 vector of sampled angles in [-pi, pi)
%       r       - n x 1 vector of amplitudes (all 1s)
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% input check

% kappa must be numeric and non-negative
if ~isnumeric(kappa)
    error('kappa must be numeric')
else
    if any(kappa<0)
        error('kappa must be non negative')
    end
end

% mu must be a value in [-pi, pi)
if ~isnumeric(mu)
    error('mu must be numeric')
else
    if any(mu<-pi) || any(mu>pi)
        error('mu must be between -pi and pi')
    end
end

% kappa and mu must have the same number of elements
if numel(kappa) ~= numel(mu)
    error('kappa and mu must have the same number of elements');
end

% n must be a number greater than 0
if ~isnumeric(n)
    error('n must be numeric')
else
    if n <= 0
        error('n must be greater than 0')
    end
end

%% implementation

% Handle mixture case: if multiple distributions, call po_vonMisesMixSample
if numel(mu) > 1
    [theta, r] = po_vonMisesMixSample(kappa, mu, n);
    return
end


if kappa == 0
    theta = (2*pi) * rand(n,1) - pi;   % uniform in [-pi, pi)
    r = ones(n,1);
    return
end

% use Best & Fisher (1979) rejection sampling algorithm
a = 1 + sqrt(1 + 4 * kappa^2);
b = (a - sqrt(2 * a)) / (2 * kappa);
r_const = (1 + b^2) / (2 * b);

theta = zeros(n,1);
for i = 1:n
    while true
        u1 = rand;
        z = cos(pi * u1);
        f = (1 + r_const * z) / (r_const + z);
        c = kappa * (r_const - f);

        u2 = rand;
        if (u2 <= c * (2 - c)) || (log(u2) <= (1 - c))
            break;
        end
    end

    f = max(-1, min(1, f));            % clamp for numerical safety
    sgn = 2 * (rand < 0.5) - 1;    % +1 or -1 with equal probability
    theta(i) = mu + sgn * acos(f);
end

% wrap sampler output to [-pi, pi)
theta = mod(theta + pi, 2*pi) - pi;

r = ones(n,1);

end % function end

