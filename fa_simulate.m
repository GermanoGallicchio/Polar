function [theta, r] = fa_simulate(fa_cfg)
% SYNTAX:
%        [theta, r] = fa_simulate(fa_cfg)
%
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
%   fa_cfg - Configuration structure with the following fields:
%
%     fa_cfg.simulate.theta  - [logical] true to sample theta from distribution (modes 1,3)
%                                        false for uniform theta (mode 2)
%     fa_cfg.simulate.r      - [logical] true to compute r from density (modes 2,3)
%                                        false to set r = 1 (mode 1)
%
%     Mode selection:
%       Mode 1: fa_cfg.simulate.theta = true,  fa_cfg.simulate.r = false
%       Mode 2: fa_cfg.simulate.theta = false, fa_cfg.simulate.r = true
%       Mode 3: fa_cfg.simulate.theta = true,  fa_cfg.simulate.r = true
%
%     fa_cfg.simulate.nSamples            - [numeric] number of samples to generate
%     fa_cfg.simulate.distribution.family - [string] 'vonMises' or 'wrappedCauchy'
%     For vonMises distributions:
%       fa_cfg.simulate.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       fa_cfg.simulate.distribution.theta.kappa - [numeric] concentration parameter(s) for theta (kappa >= 0)
%       fa_cfg.simulate.distribution.r.mu        - [numeric] mean direction(s) for r density
%       fa_cfg.simulate.distribution.r.kappa     - [numeric] concentration parameter(s) for r density (kappa >= 0)
%     For wrappedCauchy distributions:
%       fa_cfg.simulate.distribution.theta.mu    - [numeric] mean direction(s) for theta
%       fa_cfg.simulate.distribution.theta.rho   - [numeric] concentration parameter(s) for theta (0 <= rho < 1)
%       fa_cfg.simulate.distribution.r.mu        - [numeric] mean direction(s) for r density
%       fa_cfg.simulate.distribution.r.rho       - [numeric] concentration parameter(s) for r density (0 <= rho < 1)
%     fa_cfg.simulate.r_noise             - [numeric] (optional, default=0) multiplicative Gaussian noise level for r (std of noise)
%     fa_cfg.simulate.randomSeed          - [numeric] (optional, default=42) random seed
%     fa_cfg.verbose                      - [logical] (optional, default=true) display warnings
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
    fa_cfg.simulation.nSamples = 400;
    fa_cfg.simulation.distribution.family = 'wrappedCauchy';  % 'vonMises' | 'wrappedCauchy'
    fa_cfg.simulation.theta = true;   % choose if simulating theta
    fa_cfg.simulation.r     = true;   % choose if simulating r
    %fa_cfg.simulation.distribution.theta.kappa  = [1]; % concentration parameter for von Mises
    fa_cfg.simulation.distribution.theta.rho    = 0.8; % concentration parameter for wrapped Cauchy
    fa_cfg.simulation.distribution.theta.mu     = [pi/2]; % direction parameter for either von Mises or wrapped Cauchy
    %fa_cfg.simulation.distribution.r.kappa  = 3; % concentration parameter for von Mises
    fa_cfg.simulation.distribution.r.rho    = 0.5; % concentration parameter for wrapped Cauchy
    fa_cfg.simulation.distribution.r.mu     = pi; % direction parameter for either von Mises or wrapped Cauchy
end

%% input check and defaults

% if verbose not set by user, set it to true
if ~isfield(fa_cfg,'verbose')
    fa_cfg.verbose = true;
    warning('fa_cfg.verbose not defined by user. i am setting it to true, but you can specify it explicitly to avoid this warning. This is not a bug. See fa_documentation() for details')
end

% sanity check: fa_cfg.simulation field exists
if ~isfield(fa_cfg,'simulation')
    error('fa_cfg structure needs to have a "simulation" field. see fa_documentation() for details')
end

% sanity check: fa_cfg.simulation.distribution field exists
if ~isfield(fa_cfg.simulation,'distribution')
    error('fa_cfg.simulation structure needs to have a "distribution" field. see fa_documentation() for details')
end

% sanity check: fa_cfg.simulation.distribution.family field exists
if ~isfield(fa_cfg.simulation.distribution,'family')
    error('fa_cfg.simulation.distribution structure needs to have a "family" field. see fa_documentation() for details')
end

% sanity check: simulation distribution family is one of the allowed ones
familyDistributions = ["vonMises" "wrappedCauchy"];
if ~ismember(fa_cfg.simulation.distribution.family, familyDistributions)
    disp(familyDistributions)
    error('fa_cfg.simulate.distribution.family needs to be one of the families above. see fa_documentation() for details')
end



% sanity check: fa_cfg.simulation.theta field exists 
% i.e., do you want to simulate theta?
if ~isfield(fa_cfg.simulation,'theta')
    error('fa_cfg.simulation structure needs to have a "theta" field. see fa_documentation() for details')
end
% sanity check: fa_cfg.simulation.theta is logical
if ~islogical(fa_cfg.simulation.theta)
    error('fa_cfg.simulation.theta needs to be a logical value (true = you want to simulate it). see fa_documentation() for details')
end

% sanity check: if simulating theta, theta distribution parameters exist
if fa_cfg.simulation.theta
    if ~isfield(fa_cfg.simulation.distribution,'theta')
        error('fa_cfg.simulation.distribution structure needs to have a "theta" field. see fa_documentation() for details')
    end
end



% sanity check: fa_cfg.simulation.r field exists 
% i.e., do you want to simulate r?
if ~isfield(fa_cfg.simulation,'r')
    error('fa_cfg.simulation structure needs to have a "r" field. see fa_documentation() for details')
end

% sanity check: fa_cfg.simulation.r is logical
if ~islogical(fa_cfg.simulation.r)
    error('fa_cfg.simulation.r needs to be a logical value (true = you want to simulate it). see fa_documentation() for details')
end

% sanity check: if simulating r, r distribution parameters exist
if fa_cfg.simulation.r
    if ~isfield(fa_cfg.simulation.distribution,'r')
        error('fa_cfg.simulation.distribution structure needs to have a "r" field. see fa_documentation() for details')
    end
end

% r_noise = 0 (no noise by default)
% Multiplicative Gaussian noise applied to r (density values) to make the process stochastic.
% This adds trial-to-trial variability in magnitudes, simulating realistic noise in data
% (e.g., for phase-amplitude coupling simulations with varying amplitude across trials).
if fa_cfg.simulation.r
    if ~isfield(fa_cfg.simulation,'r_noise')
        fa_cfg.simulation.r_noise = 0;
        if fa_cfg.verbose
            warning(['fa_cfg.simulation.r_noise was not defined. i am setting it to ' num2str(fa_cfg.simulation.r_noise) ' (0 = no noise), but you can specify it explicitly to avoid this warning. This is not a bug. See fa_documentation() for details'])
        end
    end
end




% sanity check: if both fa_cfg.simulate.theta and fa_cfg.simulate.theta are false, there is nothing to sample
if fa_cfg.simulation.theta==false  &&  fa_cfg.simulation.r==false
    disp('you set both fa_cfg.simulation.theta and fa_cfg.simulation.r to false, so there is nothing to simulate. see fa_documentation() for details')
    return
end

% sanity check: the right parameters are given for the distribution
switch fa_cfg.simulation.distribution.family
    case 'vonMises' % vonMises distribution requires parameters kappa and mu

        if fa_cfg.simulation.theta
            if ~isfield(fa_cfg.simulation.distribution.theta,'kappa')
                error('fa_cfg.simulation.distribution.theta structure needs to have a "kappa" field (concentration). see fa_documentation() for details')
            end
            if ~isfield(fa_cfg.simulation.distribution.theta,'mu')
                error('fa_cfg.simulation.distribution.theta structure needs to have a "mu" field (direction). see fa_documentation() for details')
            end
        end

        if fa_cfg.simulation.r
            if ~isfield(fa_cfg.simulation.distribution.r,'kappa')
                error('fa_cfg.simulation.distribution.r structure needs to have a "kappa" field (concentration). see fa_documentation() for details')
            end
            if ~isfield(fa_cfg.simulation.distribution.r,'mu')
                error('fa_cfg.simulation.distribution.r structure needs to have a "mu" field (direction). see fa_documentation() for details')
            end
        end

    case 'wrappedCauchy' % wrappedCauchy distribution requires parameters rho and mu

        if fa_cfg.simulation.theta
            if ~isfield(fa_cfg.simulation.distribution.theta,'rho')
                error('fa_cfg.simulation.distribution.theta structure needs to have a "rho" field (concentration). see fa_documentation() for details')
            end
            if ~isfield(fa_cfg.simulation.distribution.theta,'mu')
                error('fa_cfg.simulation.distribution.theta structure needs to have a "mu" field (direction). see fa_documentation() for details')
            end
        end

        if fa_cfg.simulation.r
            if ~isfield(fa_cfg.simulation.distribution.r,'rho')
                error('fa_cfg.simulation.distribution.r structure needs to have a "rho" field (concentration). see fa_documentation() for details')
            end
            if ~isfield(fa_cfg.simulation.distribution.r,'mu')
                error('fa_cfg.simulation.distribution.r structure needs to have a "mu" field (direction). see fa_documentation() for details')
            end
        end

end




    



% sanity check: fa_cfg.simulation.nSamples field exists
% i.e., specified how many samples to draw from this distribution
if ~isfield(fa_cfg.simulation,'nSamples')
    error('fa_cfg.simulation structure needs to have a "nSamples" field. see fa_documentation() for details')
end
% sanity check: nSamples is a number
if ~isnumeric(fa_cfg.simulation.nSamples)
    error('fa_cfg.simulate.nSamples must be numeric. see fa_documentation() for details')
end


% if randomSeed not set by user, set it to 42
% (randomSeed used to make stochastic simulations replicable)
if ~isfield(fa_cfg.simulation,'randomSeed')
    fa_cfg.simulation.randomSeed = 42;
    warning(['fa_cfg.simulation.randomSeed not defined by user. i am setting it to ' num2str(fa_cfg.simulation.randomSeed) ' , but you can specify it explicitly to avoid this warning. This is not a bug. See fa_documentation() for details'])
end



%% shortcuts 

% for distribution parameters
nSamples    = fa_cfg.simulation.nSamples;
familyDistr = fa_cfg.simulation.distribution.family;
switch familyDistr
    case 'vonMises'

        if fa_cfg.simulation.theta
            theta_concentration = fa_cfg.simulation.distribution.theta.kappa;
            theta_direction = fa_cfg.simulation.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end

        if fa_cfg.simulation.r
            r_concentration = fa_cfg.simulation.distribution.r.kappa;
            r_direction = fa_cfg.simulation.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end

    case 'wrappedCauchy'

        if fa_cfg.simulation.theta
            theta_concentration = fa_cfg.simulation.distribution.theta.rho;
            theta_direction = fa_cfg.simulation.distribution.theta.mu;
        else
            theta_concentration = 'none';
            theta_direction = 'none';
        end
        
        if fa_cfg.simulation.r
            r_concentration = fa_cfg.simulation.distribution.r.rho;
            r_direction = fa_cfg.simulation.distribution.r.mu;
        else
            r_concentration = 'none';
            r_direction = 'none';
        end
end

%% feedback to user on what they want to sample

% TO DO: improve feedback by creating a table to allow distribution
% mixtures

if fa_cfg.verbose
    fprintf('\n')
    if fa_cfg.simulation.theta  &&  ~fa_cfg.simulation.r % only theta
        textLbl = "you want to simulate %u values of theta by sampling from a %s distribution with parameters %f (concentration) and %f (direction)";
        fprintf(textLbl, ...
            nSamples, familyDistr, theta_concentration, theta_direction)

    elseif ~fa_cfg.simulation.theta  &&  fa_cfg.simulation.r  % only r
        textLbl = "you want to simulate %u values of r by determining the density values --based on a %s distribution with parameters %f (concentration) and %f (direction)-- corresponding with some (uniform) theta values .";
        % TO DO: include feedback about r_noise making the simulation of r non-deterministic and stochastic

        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction)

    elseif fa_cfg.simulation.theta  &&  fa_cfg.simulation.r  % both theta and r
        textLbl = "you want to simulate %u values of r by determining the density values --based on a %s distribution with parameters %f (concentration) and %f (direction)-- corresponding with some theta values themselves sampled from a %s distribution with parameters %f (concentration) and %f (direction)";
        % TO DO: include feedback about r_noise making the simulation of r non-deterministic and stochastic
        fprintf(textLbl, ...
            nSamples, familyDistr, r_concentration, r_direction, familyDistr, theta_concentration, theta_direction)

    end
    fprintf('\n')
end


%% implementation

rng(fa_cfg.simulation.randomSeed)


if      fa_cfg.simulation.theta  &&  ~fa_cfg.simulation.r   % simulate theta

    switch familyDistr
        case 'vonMises'
            mu    = fa_cfg.simulation.distribution.theta.mu;
            kappa = fa_cfg.simulation.distribution.theta.kappa;
            n     = fa_cfg.simulation.nSamples;
            % fa_vonMisesSample handles both single and mixture distributions internally
            [theta, ~] = fa_vonMisesSample(kappa, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1

        case 'wrappedCauchy'
            mu = fa_cfg.simulation.distribution.theta.mu;
            rho = fa_cfg.simulation.distribution.theta.rho;
            n = fa_cfg.simulation.nSamples;
            % fa_wrappedCauchySample handles both single and mixture distributions internally
            [theta, ~] = fa_wrappedCauchySample(rho, mu, n);
            r = ones(n, 1);  % unit magnitude for Case 1
    end


    % ensure vectors are column
    theta = theta(:);
    r = r(:);


elseif ~fa_cfg.simulation.theta  &&   fa_cfg.simulation.r   % simulate r, for uniform theta

    switch familyDistr
        case 'vonMises'
            mu    = fa_cfg.simulation.distribution.r.mu;
            kappa = fa_cfg.simulation.distribution.r.kappa;
            n     = fa_cfg.simulation.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = fa_vonMisesDensity(kappa, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if fa_cfg.simulation.r_noise > 0

                r = r .* (1 + fa_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            mu    = fa_cfg.simulation.distribution.r.mu;
            rho   = fa_cfg.simulation.distribution.r.rho;
            n     = fa_cfg.simulation.nSamples;
            theta = (-pi) + (0:(n-1))' * (2*pi / n);
            [~, r] = fa_wrappedCauchyDensity(rho, mu, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if fa_cfg.simulation.r_noise > 0

                r = r .* (1 + fa_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end


elseif  fa_cfg.simulation.theta  &&   fa_cfg.simulation.r   % simulate both theta and r

    switch familyDistr
        case 'vonMises'
            % Sample theta from its distribution
            mu_theta    = fa_cfg.simulation.distribution.theta.mu;
            kappa_theta = fa_cfg.simulation.distribution.theta.kappa;
            n           = fa_cfg.simulation.nSamples;
            [theta, ~]  = fa_vonMisesSample(kappa_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r    = fa_cfg.simulation.distribution.r.mu;
            kappa_r = fa_cfg.simulation.distribution.r.kappa;
            [~, r] = fa_vonMisesDensity(kappa_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if fa_cfg.simulation.r_noise > 0
                r = r .* (1 + fa_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
            
        case 'wrappedCauchy'
            % Sample theta from its distribution
            mu_theta  = fa_cfg.simulation.distribution.theta.mu;
            rho_theta = fa_cfg.simulation.distribution.theta.rho;
            n         = fa_cfg.simulation.nSamples;
            [theta, ~] = fa_wrappedCauchySample(rho_theta, mu_theta, n);
            
            % Evaluate r density at the sampled theta values
            mu_r  = fa_cfg.simulation.distribution.r.mu;
            rho_r = fa_cfg.simulation.distribution.r.rho;
            [~, r] = fa_wrappedCauchyDensity(rho_r, mu_r, theta);
            
            % Apply multiplicative Gaussian noise to r to make the process stochastic
            % (adds trial-to-trial variability, simulating realistic noisy data)
            if fa_cfg.simulation.r_noise > 0
                r = r .* (1 + fa_cfg.simulation.r_noise * randn(size(r)));
                r = max(r, 0);  % ensure non-negative
            end
    end
    
    % ensure vectors are column
    theta = theta(:);
    r = r(:);


end
