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
%     po_cfg.simulate.theta  - [logical] true to sample theta from distribution (modes 1,3)
%                                        false for uniform theta (mode 2)
%     po_cfg.simulate.r      - [logical] true to compute r from density (modes 2,3)
%                                        false to set r = 1 (mode 1)
%
%     Mode selection:
%       Mode 1: po_cfg.simulate.theta = true,  po_cfg.simulate.r = false
%       Mode 2: po_cfg.simulate.theta = false, po_cfg.simulate.r = true
%       Mode 3: po_cfg.simulate.theta = true,  po_cfg.simulate.r = true
%
%     po_cfg.simulate.nSamples            - [numeric] number of samples to generate
%     po_cfg.simulate.distribution.family - [string] 'vonMises' or 'wrappedCauchy'
%     For vonMises distributions:
%       po_cfg.simulate.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       po_cfg.simulate.distribution.theta.kappa - [numeric] concentration parameter(s) for theta (kappa >= 0)
%       po_cfg.simulate.distribution.r.mu        - [numeric] mean direction(s) for r density
%       po_cfg.simulate.distribution.r.kappa     - [numeric] concentration parameter(s) for r density (kappa >= 0)
%     For wrappedCauchy distributions:
%       po_cfg.simulate.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       po_cfg.simulate.distribution.theta.rho   - [numeric] concentration parameter(s) for theta (0 <= rho < 1)
%       po_cfg.simulate.distribution.r.mu        - [numeric] mean direction(s) for r density
%       po_cfg.simulate.distribution.r.rho       - [numeric] concentration parameter(s) for r density (0 <= rho < 1)
%     po_cfg.simulate.r_noise             - [numeric] (optional, default=0) multiplicative Gaussian noise level for r (std of noise)
%     po_cfg.simulate.randomSeed          - [numeric] (optional, default=42) random seed
%     po_cfg.verbose                      - [logical] (optional, default=true) display warnings
%
% OUTPUT:
%   theta  -    vector angles (in radians)
%   r      -    vector lengths 
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% playground section (try some parameters)

% playground: set to true only for manual testing
toggle = false;
if toggle
    % decide sampling parameters
    po_cfg.simulation.nSamples = 400;
    po_cfg.simulation.distribution.family = 'wrappedCauchy';  % 'vonMises' | 'wrappedCauchy'
    po_cfg.simulation.theta = true;   % choose if simulating theta
    po_cfg.simulation.r     = true;   % choose if simulating r
    %po_cfg.simulate.distribution.theta.kappa  = [1]; % concentration parameter for von Mises
    po_cfg.simulation.distribution.theta.rho    = 0.8; % concentration parameter for wrapped Cauchy
    po_cfg.simulation.distribution.theta.mu     = [pi/2]; % direction parameter for either von Mises or wrapped Cauchy
    %po_cfg.simulate.distribution.r.kappa  = 3; % concentration parameter for von Mises
    po_cfg.simulation.distribution.r.rho    = 0.5; % concentration parameter for wrapped Cauchy
    po_cfg.simulation.distribution.r.mu     = pi; % direction parameter for either von Mises or wrapped Cauchy
end

%% defaults

% verbose = true
fieldLbl = fieldnames(po_cfg);
if ~any(strcmp(fieldLbl,'verbose'))
    po_cfg.verbose = true;
    warning('po_cfg.verbose not defined by user. i am setting it to true, but you can specify it explicitly to avoid this warning')
end

% randomSeed (used to make stochastic simulations replicable)
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'randomSeed'))
    po_cfg.simulation.randomSeed = 42;
    warning(['po_cfg.simulate.randomSeed not defined by user. i am setting it to ' num2str(po_cfg.simulation.randomSeed) ' , but you can specify it explicitly to avoid this warning'])
end

% r_noise = 0 (no noise by default)
% Multiplicative Gaussian noise applied to r (density values) to make the process stochastic.
% This adds trial-to-trial variability in magnitudes, simulating realistic noise in data
% (e.g., for phase-amplitude coupling simulations with varying amplitude across trials).
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'r_noise'))
    po_cfg.simulation.r_noise = 0;
    if po_cfg.verbose
        warning(['po_cfg.simulate.r_noise not defined by user. i am setting it to ' num2str(po_cfg.simulation.r_noise) ' (no noise), but you can specify it explicitly to avoid this warning'])
    end
end

    
%% input checks

% sanity check: a "simulation" field exists
fieldLbl = fieldnames(po_cfg);
if ~any(strcmp(fieldLbl,'simulation'))
    error('po_cfg structure needs to have a "simulation" field')
end

% sanity check: simulation "distribution" field exists
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'distribution'))
    error('po_cfg.simulate structure needs to have a "distribution" field')
end

% sanity check: simulation distribution "family" field exists
fieldLbl = fieldnames(po_cfg.simulation.distribution);
if ~any(strcmp(fieldLbl,'family'))
    error('po_cfg.simulate.distribution structure needs to have a "family" field')
end

% sanity check: simulation distribution family is one of the allowed ones
familyDistributions = ["vonMises" "wrappedCauchy"];
if ~ismember(po_cfg.simulation.distribution.family, familyDistributions)
    disp(familyDistributions)
    error('po_cfg.simulate.distribution.family needs to be one of the above')
end

% sanity check: user has decided if wanting to sample theta
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'theta'))
    error('po_cfg.simulate structure needs to have a "theta" field')
end
% sanity check: theta needs to be a logical value
if ~islogical(po_cfg.simulation.theta)
    error('theta needs to be a logical value (true = you want to sample it')
end

% sanity check: user has decided if wanting to sample r
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'r'))
    error('po_cfg.simulate structure needs to have a "r" field')
end
% sanity check: r needs to be a logical value
if ~islogical(po_cfg.simulation.r)
    error('r needs to be a logical value (true = you want to sample it)')
end

% sanity check: if both po_cfg.simulate.theta and po_cfg.simulate.theta are false, there is nothing to sample
if po_cfg.simulation.theta==false  &&  po_cfg.simulation.r==false
    error('nothing to sample: both po_cfg.simulate.theta and po_cfg.simulate.r are set to false')
end

% sanity check: check the right parameters are given for the distribution
switch po_cfg.simulation.distribution.family
    case 'vonMises' % vonMises distribution requires parameters kappa and mu

        if po_cfg.simulation.theta
            fieldLbl = fieldnames(po_cfg.simulation.distribution.theta);
            if ~any(strcmp(fieldLbl,'kappa'))
                error('po_cfg.simulate.distribution.theta.kappa is required')
            end
            if ~any(strcmp(fieldLbl,'mu'))
                error('po_cfg.simulate.distribution.theta.mu is required')
            end
        end

        if po_cfg.simulation.r
            fieldLbl = fieldnames(po_cfg.simulation.distribution.r);
            if ~any(strcmp(fieldLbl,'kappa'))
                error('po_cfg.simulate.distribution.r.kappa is required')
            end
            if ~any(strcmp(fieldLbl,'mu'))
                error('po_cfg.simulate.distribution.r.mu is required')
            end
        end

    case 'wrappedCauchy' % wrappedCauchy distribution requires parameters rho and mu

        if po_cfg.simulation.theta
            fieldLbl = fieldnames(po_cfg.simulation.distribution.theta);
            if ~any(strcmp(fieldLbl,'rho'))
                error('po_cfg.simulate.distribution.theta.rho is required')
            end
            if ~any(strcmp(fieldLbl,'mu'))
                error('po_cfg.simulate.distribution.theta.mu is required')
            end
        end

        if po_cfg.simulation.r
            fieldLbl = fieldnames(po_cfg.simulation.distribution.r);
            if ~any(strcmp(fieldLbl,'rho'))
                error('po_cfg.simulate.distribution.r.rho is required')
            end
            if ~any(strcmp(fieldLbl,'mu'))
                error('po_cfg.simulate.distribution.r.mu is required')
            end
        end

end


% sanity check: specified how many samples to draw from this distribution
fieldLbl = fieldnames(po_cfg.simulation);
if ~any(strcmp(fieldLbl,'nSamples'))
    error('po_cfg.simulate structure needs to have a "nSamples" field')
end
% sanity check: nSamples is a number
if ~isnumeric(po_cfg.simulation.nSamples)
    error('po_cfg.simulate.nSamples must be numeric')
end



%% shortcuts 

% for distribution parameters
nSamples    = po_cfg.simulation.nSamples;
familyDistr = po_cfg.simulation.distribution.family;
randomSeed = po_cfg.simulation.randomSeed;
switch familyDistr
    case 'vonMises'

        if po_cfg.simulation.theta
            theta_concentration = po_cfg.simulation.distribution.theta.kappa;
            theta_direction = po_cfg.simulation.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end

        if po_cfg.simulation.r
            r_concentration = po_cfg.simulation.distribution.r.kappa;
            r_direction = po_cfg.simulation.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end

    case 'wrappedCauchy'

        if po_cfg.simulation.theta
            theta_concentration = po_cfg.simulation.distribution.theta.rho;
            theta_direction = po_cfg.simulation.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end
        
        if po_cfg.simulation.r
            r_concentration = po_cfg.simulation.distribution.r.rho;
            r_direction = po_cfg.simulation.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end
end

%% feedback to user on what they want to sample
if po_cfg.verbose
    fprintf('\n')
    if po_cfg.simulation.theta  &&  ~po_cfg.simulation.r % only theta
        textLbl = "you want to simulate %u values of theta by sampling from a %s distribution with parameters %f (concentration) and %f (direction)";
        fprintf(textLbl, ...
            nSamples, familyDistr, theta_concentration, theta_direction)

    elseif ~po_cfg.simulation.theta  &&  po_cfg.simulation.r  % only r
        textLbl = "you want to simulate %u values of r by determining the density values --based on a %s distribution with parameters %f (concentration) and %f (direction)-- corresponding with some (uniform) theta values .";
        % TO DO: include feedback about r_noise making the simulation of r non-deterministic and stochastic

        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction)

    elseif po_cfg.simulation.theta  &&  po_cfg.simulation.r  % both theta and r
        textLbl = "you want to simulate %u values of r by determining the density values --based on a %s distribution with parameters %f (concentration) and %f (direction)-- corresponding with some theta values themselves sampled from a %s distribution with parameters %f (concentration) and %f (direction)";
        % TO DO: include feedback about r_noise making the simulation of r non-deterministic and stochastic
        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction, familyDistr, theta_concentration, theta_direction)

    end
    fprintf('\n')
end

%% implementation

rng(randomSeed)

if      po_cfg.simulation.theta  &&  ~po_cfg.simulation.r   % simulate theta

    switch familyDistr
        case 'vonMises'
            mu    = po_cfg.simulation.distribution.theta.mu;
            kappa = po_cfg.simulation.distribution.theta.kappa;
            n     = po_cfg.simulation.nSamples;
            % po_vonMisesSample handles both single and mixture distributions internally
            [theta, ~] = po_vonMisesSample(kappa, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1

        case 'wrappedCauchy'
            mu = po_cfg.simulation.distribution.theta.mu;
            rho = po_cfg.simulation.distribution.theta.rho;
            n = po_cfg.simulation.nSamples;
            % po_wrappedCauchySample handles both single and mixture distributions internally
            [theta, ~] = po_wrappedCauchySample(rho, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1
    end


    % ensure vectors are column
    theta = theta(:);
    r = r(:);


elseif ~po_cfg.simulation.theta  &&   po_cfg.simulation.r   % simulate r, for uniform theta

    switch familyDistr
        case 'vonMises'
            mu    = po_cfg.simulation.distribution.r.mu;
            kappa = po_cfg.simulation.distribution.r.kappa;
            n     = po_cfg.simulation.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = po_vonMisesDensity(kappa, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.simulation.r_noise > 0
                r = r .* (1 + po_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            mu    = po_cfg.simulation.distribution.r.mu;
            rho   = po_cfg.simulation.distribution.r.rho;
            n     = po_cfg.simulation.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = po_wrappedCauchyDensity(rho, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.simulation.r_noise > 0
                r = r .* (1 + po_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end


elseif  po_cfg.simulation.theta  &&   po_cfg.simulation.r   % simulate both theta and r

    switch familyDistr
        case 'vonMises'
            % Sample theta from its distribution
            mu_theta    = po_cfg.simulation.distribution.theta.mu;
            kappa_theta = po_cfg.simulation.distribution.theta.kappa;
            n           = po_cfg.simulation.nSamples;
            [theta, ~]  = po_vonMisesSample(kappa_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r    = po_cfg.simulation.distribution.r.mu;
            kappa_r = po_cfg.simulation.distribution.r.kappa;
            [~, r] = po_vonMisesDensity(kappa_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.simulation.r_noise > 0
                r = r .* (1 + po_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            % Sample theta from its distribution
            mu_theta  = po_cfg.simulation.distribution.theta.mu;
            rho_theta = po_cfg.simulation.distribution.theta.rho;
            n         = po_cfg.simulation.nSamples;
            [theta, ~] = po_wrappedCauchySample(rho_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r  = po_cfg.simulation.distribution.r.mu;
            rho_r = po_cfg.simulation.distribution.r.rho;
            [~, r] = po_wrappedCauchyDensity(rho_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if po_cfg.simulation.r_noise > 0
                r = r .* (1 + po_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end
    
    % ensure vectors are column
    theta = theta(:);
    r = r(:);


end
