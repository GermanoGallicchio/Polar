function [theta, r] = po_simulate(po_cfg)

% DESCRIPTION:
%   Simulate polar data (angles and magnitudes) from circular distributions.
%
%   Three simulation modes:
%   1. Simulate theta by sampling stochastically from a theoretical distribution, with r = 1
%      (for phase consistency measures like ITPC, PLV)
%   2. Determine r using distribution density for uniform values of theta
%      (for distribution visualization and shape analysis)
%   3. Determine r for theta themselves sampled from another distribution
%      (for phase-amplitude coupling simulations with coupled distributions)
%
%   Stochasticity:
%   - case 1: theta is sampled stochastically from the theoretical distribution
%   - case 2: r is deterministic (density evaluations) for uniform theta, unless noise is added
%   - case 3: theta is sampled stochastically from its distribution like in case 1, r is deterministically like in case 2, unless noise is added
%
% INPUT:        
%   po_cfg - Configuration structure with the following fields:
%
%     po_cfg.sampling.theta  - [logical] true to sample theta from distribution (modes 1,3)
%                                        false for uniform theta (mode 2)
%     po_cfg.sampling.r      - [logical] true to compute r from density (modes 2,3)
%                                        false to set r = 1 (mode 1)
%
%     Mode selection:
%       Mode 1: po_cfg.sampling.theta = true,  po_cfg.sampling.r = false
%       Mode 2: po_cfg.sampling.theta = false, po_cfg.sampling.r = true
%       Mode 3: po_cfg.sampling.theta = true,  po_cfg.sampling.r = true
%
%     po_cfg.sampling.nSamples            - [numeric] number of samples to generate
%     po_cfg.sampling.distribution.family - [string] 'vonMises' or 'wrappedCauchy'
%     
%     For vonMises distributions:
%       po_cfg.sampling.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       po_cfg.sampling.distribution.theta.kappa - [numeric] concentration parameter(s) 
%                                                   for theta (kappa >= 0)
%       po_cfg.sampling.distribution.r.mu        - [numeric] mean direction(s) for r density
%       po_cfg.sampling.distribution.r.kappa     - [numeric] concentration parameter(s) 
%                                                   for r density (kappa >= 0)
%     
%     For wrappedCauchy distributions:
%       po_cfg.sampling.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       po_cfg.sampling.distribution.theta.rho   - [numeric] concentration parameter(s) 
%                                                   for theta (0 <= rho < 1)
%       po_cfg.sampling.distribution.r.mu        - [numeric] mean direction(s) for r density
%       po_cfg.sampling.distribution.r.rho       - [numeric] concentration parameter(s) 
%                                                   for r density (0 <= rho < 1)
%     
%     po_cfg.sampling.r_noise             - [numeric] (optional, default=0) multiplicative 
%                                            Gaussian noise level for r (std of noise)
%     po_cfg.sampling.randomSeed          - [numeric] (optional, default=42) random seed
%     po_cfg.verbose                      - [logical] (optional, default=true) display warnings
%
%
% OUTPUT:
%   theta  -    vector angles (in radians)
%   r      -    vector lengths 
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% playground section (try some parameters)

toggle = true;

if toggle

% decide sampling parameters
po_cfg.sampling.nSamples = 400;
po_cfg.sampling.distribution.family = 'wrappedCauchy';  % 'vonMises' | 'wrappedCauchy'
po_cfg.sampling.theta = true;   % choose if sampling theta
po_cfg.sampling.r     = true;  % choose if sampling r
%po_cfg.sampling.distribution.theta.kappa  = [1]; % concentration parameter for von Mises
po_cfg.sampling.distribution.theta.rho    = 0.8; % concentration parameter for wrapped Cauchy
po_cfg.sampling.distribution.theta.mu     = [pi/2]; % direction parameter for either von Mises or wrapped Cauchy
%po_cfg.sampling.distribution.r.kappa  = 3; % concentration parameter for von Mises
po_cfg.sampling.distribution.r.rho    = 0.5; % concentration parameter for wrapped Cauchy
po_cfg.sampling.distribution.r.mu     = pi; % direction parameter for either von Mises or wrapped Cauchy


end

%% defaults

% verbose = true
fieldLbl = fieldnames(po_cfg);
if any(strcmp(fieldLbl,'verbose')) ~= 1
    po_cfg.verbose = true;
    warning('po_cfg.verbose not defined by user. i am setting it to true, but you can specify it explicitly to avoid this warning')
end

% randomSeed
fieldLbl = fieldnames(po_cfg);
if any(strcmp(fieldLbl,'randomSeed')) ~= 1
    po_cfg.sampling.randomSeed = 42;
    warning(['po_cfg.sampling.randomSeed not defined by user. i am setting it to ' num2str(po_cfg.sampling.randomSeed) ' , but you can specify it explicitly to avoid this warning'])
end

% r_noise = 0 (no noise by default)
% Multiplicative Gaussian noise applied to r (density values) to make the process stochastic.
% This adds trial-to-trial variability in magnitudes, simulating realistic noise in data
% (e.g., for phase-amplitude coupling simulations with varying amplitude across trials).
fieldLbl = fieldnames(po_cfg.sampling);
if any(strcmp(fieldLbl,'r_noise')) ~= 1
    po_cfg.sampling.r_noise = 0;
    if po_cfg.verbose
        warning(['po_cfg.sampling.r_noise not defined by user. i am setting it to ' num2str(po_cfg.sampling.r_noise) ' (no noise), but you can specify it explicitly to avoid this warning'])
    end
end

    
%% input checks

% sanity check: a "sampling" field exists
fieldLbl = fieldnames(po_cfg);
if any(strcmp(fieldLbl,'sampling')) ~= 1
    error('po_cfg structure needs to have a "sampling" field')
end

% sanity check: sampling "distribution" field exists
fieldLbl = fieldnames(po_cfg.sampling);
if any(strcmp(fieldLbl,'distribution')) ~= 1
    error('po_cfg.sampling structure needs to have a "distribution" field')
end

% sanity check: sampling distribution "family" field exists
fieldLbl = fieldnames(po_cfg.sampling.distribution);
if any(strcmp(fieldLbl,'family')) ~= 1
    error('po_cfg.sampling.distribution structure needs to have a "family" field')
end

% sanity check: sampling distribution family is one of the allowed ones
familyDistributions = ["vonMises" "wrappedCauchy"];
if ~ismember(po_cfg.sampling.distribution.family, familyDistributions)
    disp(familyDistributions)
    error('po_cfg.sampling.distribution.family needs to be one of the above')
end

% sanity check: user has decided if wanting to sample theta
fieldLbl = fieldnames(po_cfg.sampling);
if any(strcmp(fieldLbl,'theta')) ~= 1
    error('po_cfg.sampling structure needs to have a "theta" field')
end
% sanity check: theta needs to be a logical value
if ~islogical(po_cfg.sampling.theta)
    error('theta needs to be a logical value (true = you want to sample it')
end

% sanity check: user has decided if wanting to sample r
fieldLbl = fieldnames(po_cfg.sampling);
if any(strcmp(fieldLbl,'r')) ~= 1
    error('po_cfg.sampling structure needs to have a "r" field')
end
% sanity check: r needs to be a logical value
if ~islogical(po_cfg.sampling.r)
    error('theta needs to be a logical value (true = you want to sample it')
end

% sanity check: if both po_cfg.sampling.theta and po_cfg.sampling.theta are false, there is nothing to sample
if po_cfg.sampling.theta==false  &&  po_cfg.sampling.r==false
    error('nothing to sample: both po_cfg.sampling.theta and po_cfg.sampling.r are set to false')
end

% sanity check: check the right parameters are given for the distribution
switch po_cfg.sampling.distribution.family
    case 'vonMises' % vonMises distribution requires parameters kappa and mu

        if po_cfg.sampling.theta
            fieldLbl = fieldnames(po_cfg.sampling.distribution.theta);
            if any(strcmp(fieldLbl,'kappa')) ~= 1
                error('po_cfg.sampling.distribution.theta.kappa is required')
            end
            if any(strcmp(fieldLbl,'mu')) ~= 1
                error('po_cfg.sampling.distribution.theta.mu is required')
            end
        end

        if po_cfg.sampling.r
            fieldLbl = fieldnames(po_cfg.sampling.distribution.r);
            if any(strcmp(fieldLbl,'kappa')) ~= 1
                error('po_cfg.sampling.distribution.r.kappa is required')
            end
            if any(strcmp(fieldLbl,'mu')) ~= 1
                error('po_cfg.sampling.distribution.r.mu is required')
            end
        end

    case 'wrappedCauchy' % wrappedCauchy distribution requires parameters rho and mu

        if po_cfg.sampling.theta
            fieldLbl = fieldnames(po_cfg.sampling.distribution.theta);
            if any(strcmp(fieldLbl,'rho')) ~= 1
                error('po_cfg.sampling.distribution.theta.rho is required')
            end
            if any(strcmp(fieldLbl,'mu')) ~= 1
                error('po_cfg.sampling.distribution.theta.mu is required')
            end
        end

        if po_cfg.sampling.r
            fieldLbl = fieldnames(po_cfg.sampling.distribution.r);
            if any(strcmp(fieldLbl,'rho')) ~= 1
                error('po_cfg.sampling.distribution.r.rho is required')
            end
            if any(strcmp(fieldLbl,'mu')) ~= 1
                error('po_cfg.sampling.distribution.r.mu is required')
            end
        end

end


% sanity check: specified how many samples to draw from this distribution
fieldLbl = fieldnames(po_cfg.sampling);
if any(strcmp(fieldLbl,'nSamples')) ~= 1
    error('po_cfg.sampling structure needs to have a "nSample" field')
end
% sanity check: nSamples is a number
if ~isnumeric(po_cfg.sampling.nSamples)
    error('po_cfg.sampling.nSamples must be numeric')
end



%% shortcuts 

% for distribution parameters
nSamples    = po_cfg.sampling.nSamples;
familyDistr = po_cfg.sampling.distribution.family;
randomSeed = po_cfg.sampling.randomSeed;
switch familyDistr
    case 'vonMises'

        if po_cfg.sampling.theta
            theta_concentration = po_cfg.sampling.distribution.theta.kappa;
            theta_direction = po_cfg.sampling.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end

        if po_cfg.sampling.r
            r_concentration = po_cfg.sampling.distribution.r.kappa;
            r_direction = po_cfg.sampling.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end

    case 'wrappedCauchy'

        if po_cfg.sampling.theta
            theta_concentration = po_cfg.sampling.distribution.theta.rho;
            theta_direction = po_cfg.sampling.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end
        
        if po_cfg.sampling.r
            r_concentration = po_cfg.sampling.distribution.r.rho;
            r_direction = po_cfg.sampling.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end
end

%% feedback to user on what they want to sample
if po_cfg.verbose
    fprintf('\n')
    if po_cfg.sampling.theta  &&  ~po_cfg.sampling.r % only theta
        textLbl = "you want to sample %u values of theta from a %s distribution with parameters %f (concentration) and %f (direction)";
        fprintf(textLbl, ...
            nSamples, familyDistr, theta_concentration, theta_direction)

    elseif ~po_cfg.sampling.theta  &&  po_cfg.sampling.r  % only r
        textLbl = "you want to sample %u values of r (for uniform theta) from a %s distribution with parameters %f (concentration) and %f (direction)";
        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction)

    elseif po_cfg.sampling.theta  &&  po_cfg.sampling.r  % both theta and r
        textLbl = "sampling %u values of r from a %s distribution with parameters %f (concentration) and %f (direction) and of theta from a %s distribution with parameters %f (concentration) and %f (direction)";
        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction, familyDistr, theta_concentration, theta_direction)

    end
    fprintf('\n')
end

%% implementation

% UNTIL HERE OK
rng(randomSeed)



if      po_cfg.sampling.theta  &&  ~po_cfg.sampling.r   % simulate theta

    % -- old section ---
    % drawing theta samples (step 3) from a probability mass
    % distribution obtained from an empirical dendisty distribution
    % (step 2), which is obtained from a theoretical distribution (sum
    % of von Mises).
    % % step 1. create an angular grid in (-pi, pi]
    % thetaCoverage = 10000; % hard-coded but not important. it needs to provide a good coverage of the angular space
    % thetaLim = [-pi, pi];  % currently: [-pi, pi]
    % theta = linspace(thetaLim(1), thetaLim(2), thetaCoverage+1);
    % theta(1) = [];         % now: (-pi, pi] 
    % 
    % % step 2. empirical density distribution
    % kappa = theta_concentration;
    % mu    = theta_direction;
    % vmDensity = po_vonMisesDensity(kappa,mu,theta);
    % 
    % % step 3. probability mass distribution
    % vmMass = vmDensity * (theta(2) - theta(1)); % approximate area per bin
    % vmMass = vmMass / sum(vmMass);              % make sure it sums to 1
    % 
    % % step 4. draw theta samples from probability mass distribution
    % sampleIdx = randsample(length(theta), nSamples, true, vmMass);
    % samples = theta(sampleIdx);
    % -- --


    switch familyDistr
        case 'vonMises'
            mu    = po_cfg.sampling.distribution.theta.mu;
            kappa = po_cfg.sampling.distribution.theta.kappa;
            n     = po_cfg.sampling.nSamples;
            % po_vonMisesSample handles both single and mixture distributions internally
            [theta, ~] = po_vonMisesSample(kappa, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1

        case 'wrappedCauchy'
            mu = po_cfg.sampling.distribution.theta.mu;
            rho = po_cfg.sampling.distribution.theta.rho;
            n = po_cfg.sampling.nSamples;
            % po_wrappedCauchySample handles both single and mixture distributions internally
            [theta, ~] = po_wrappedCauchySample(rho, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1
    end


    % ensure vectors are column
    theta = theta(:);
    r = r(:);


elseif ~po_cfg.sampling.theta  &&   po_cfg.sampling.r   % simulate r, for uniform theta

    switch familyDistr
        case 'vonMises'
            mu    = po_cfg.sampling.distribution.r.mu;
            kappa = po_cfg.sampling.distribution.r.kappa;
            n     = po_cfg.sampling.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = po_vonMisesDensity(kappa, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.sampling.r_noise > 0
                r = r .* (1 + po_cfg.sampling.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            mu    = po_cfg.sampling.distribution.r.mu;
            rho   = po_cfg.sampling.distribution.r.rho;
            n     = po_cfg.sampling.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = po_wrappedCauchyDensity(rho, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.sampling.r_noise > 0
                r = r .* (1 + po_cfg.sampling.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end


elseif  po_cfg.sampling.theta  &&   po_cfg.sampling.r   % simulate both theta and r

    switch familyDistr
        case 'vonMises'
            % Sample theta from its distribution
            mu_theta    = po_cfg.sampling.distribution.theta.mu;
            kappa_theta = po_cfg.sampling.distribution.theta.kappa;
            n           = po_cfg.sampling.nSamples;
            [theta, ~]  = po_vonMisesSample(kappa_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r    = po_cfg.sampling.distribution.r.mu;
            kappa_r = po_cfg.sampling.distribution.r.kappa;
            [~, r] = po_vonMisesDensity(kappa_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.sampling.r_noise > 0
                r = r .* (1 + po_cfg.sampling.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            % Sample theta from its distribution
            mu_theta  = po_cfg.sampling.distribution.theta.mu;
            rho_theta = po_cfg.sampling.distribution.theta.rho;
            n         = po_cfg.sampling.nSamples;
            [theta, ~] = po_wrappedCauchySample(rho_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r  = po_cfg.sampling.distribution.r.mu;
            rho_r = po_cfg.sampling.distribution.r.rho;
            [~, r] = po_wrappedCauchyDensity(rho_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.sampling.r_noise > 0
                r = r .* (1 + po_cfg.sampling.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end
    
    % ensure vectors are column
    theta = theta(:);
    r = r(:);


end
