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
%   - Simulations of theta are stochastic when sampled from theoretical distributions
%     (modes 1 and 3) due to random sampling.
%   - Simulations of r are deterministic (density evaluations) unless multiplicative
%     Gaussian noise is added via the controllable parameter po_cfg.sampling.r_noise.
%
% INPUT:        
%   po_cfg
%
%
%
% OPTIONAL INPUT:
%   none
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
po_cfg.sampling.distribution.family = 'vonMises';  % 'vonMises' | 'wrappedCauchy'
po_cfg.sampling.theta = true;   % choose if sampling theta
po_cfg.sampling.r     = true;  % choose if sampling r
po_cfg.sampling.distribution.theta.kappa  = [1]; % concentration parameter for von Mises
%po_cfg.sampling.distribution.theta.rho    = 0.5; % concentration parameter for wrapped Cauchy
po_cfg.sampling.distribution.theta.mu     = [pi/2]; % direction parameter for either von Mises or wrapped Cauchy
po_cfg.sampling.distribution.r.kappa  = 3; % concentration parameter for von Mises
%po_cfg.sampling.distribution.r.rho    = 0.5; % concentration parameter for wrapped Cauchy
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
            % po_wrappedCauchy handles both single and mixture distributions internally
            [theta, ~] = po_wrappedCauchy(rho, mu, n);
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
            error('not coded yet')
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
            error('not coded yet')
    end
    
    % ensure vectors are column
    theta = theta(:);
    r = r(:);


end

% temporary section
figure(101)
viewParams=struct();
po_view(theta, r, viewParams)


%% older section
switch num2str(theta_rho)
    case num2str([0 1])  % simulate rho with uniform or given theta
        % build an empirical density distribution by drawing from a
        % theoretical distribution (sum of von Mises).
    
        fieldLbl = fieldnames(po_cfg);
        if ~any(strcmp(fieldLbl,'theta'))
            thetaLim = [-pi, pi];  % create a range of phase values in (-pi, pi]
            theta = linspace(thetaLim(1), thetaLim(2), nSamples+1);
            theta(1) = [];
        else
            theta = po_cfg.theta;
        end

        % empirical density distribution 
        r = po_vonMisesDensity(kappa,mu,theta);

        % ensure vectors are column
        theta = theta(:);
        r = r(:);

    case num2str([1 0])  % simulate theta, with unit rho
        % drawing theta samples (step 3) from a probability mass
        % distribution obtained from an empirical dendisty distribution
        % (step 2), which is obtained from a theoretical distribution (sum
        % of von Mises).

        % step 1. empirical density distribution 
        thetaLim = [-pi, pi];  % create a range of phase values in (-pi, pi]
        theta = linspace(thetaLim(1), thetaLim(2), 10000+1);
        theta(1) = [];
        r = po_vonMisesDensity(kappa,mu,theta);
        % step 2. probability mass distribution
        vmMass = r * (theta(2) - theta(1)); % approximate area per bin
        vmMass = vmMass / sum(vmMass);              % make sure it sums to 1
        % step 3. draw theta samples from probability mass distribution
        sampleIdx = randsample(length(theta), nSamples, true, vmMass);
        samples = theta(sampleIdx);

        % ensure vectors are column
        theta = samples(:);
        r = ones(size(theta));

    

    otherwise  % simulate both theta and rho
        if strcmp(num2str(po_cfg.theta_rho),num2str([1 1]))  % simulate both theta and rho
            warning(['if you want to simulate both theta and rho at once,' ...
                ' do it in two steps' ...
                ' 1) simulate theta with its parameters' ...
                ' 2) simulate rho with its parameters, but also use "simulationParams.theta = theta" to inherit the simulated theta'])
        end
        error('po_cfg.simulationParams.theta_rho must be either [0 1] or [1 0]')
end













