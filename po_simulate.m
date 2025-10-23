function [theta, rho] = po_simulate(simulationParams)

% DESCRIPTION:
%
% INPUT:        
%
% OPTIONAL:
%
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% shortcuts

nSamples  = simulationParams.nSamples;
kappa     = simulationParams.kappa;
mu        = simulationParams.mu;
theta_rho = simulationParams.theta_rho;
nModes = length(mu);

%% random seed




fieldLbl = fieldnames(simulationParams);
if ~any(strcmp(fieldLbl,'randomSeed'))
    randomSeed = 42;
    warning(['setting random seed for you: ' num2str(randomSeed)])
else
    randomSeed = simulationParams.randomSeed;
end

rng(randomSeed)
%% implementation
switch num2str(theta_rho)
    case num2str([0 1])  % simulate rho with uniform or given theta
        % build an empirical density distribution by drawing from a
        % theoretical distribution (sum of von Mises).
    
        fieldLbl = fieldnames(simulationParams);
        if ~any(strcmp(fieldLbl,'theta'))
            thetaLim = [-pi, pi];  % create a range of phase values in (-pi, pi]
            theta = linspace(thetaLim(1), thetaLim(2), nSamples+1);
            theta(1) = [];
        else
            theta = simulationParams.theta;
        end

        % empirical density distribution 
        rho = po_vonMisesDensity(kappa,mu,theta);

        % ensure vectors are column
        theta = theta(:);
        rho = rho(:);

    case num2str([1 0])  % simulate theta, with unit rho
        % drawing theta samples (step 3) from a probability mass
        % distribution obtained from an empirical dendisty distribution
        % (step 2), which is obtained from a theoretical distribution (sum
        % of von Mises).

        % step 1. empirical density distribution 
        thetaLim = [-pi, pi];  % create a range of phase values in (-pi, pi]
        theta = linspace(thetaLim(1), thetaLim(2), 10000+1);
        theta(1) = [];
        rho = po_vonMisesDensity(kappa,mu,theta);
        % step 2. probability mass distribution
        pmf = rho * (theta(2) - theta(1)); % approximate area per bin
        pmf = pmf / sum(pmf);              % make sure it sums to 1
        % step 3. draw theta samples from probability mass distribution
        sampleIdx = randsample(length(theta), nSamples, true, pmf);
        samples = theta(sampleIdx);

        % ensure vectors are column
        theta = samples(:);
        rho = ones(size(theta));

    

    otherwise  % simulate both theta and rho
        if strcmp(num2str(simulationParams.theta_rho),num2str([1 1]))  % simulate both theta and rho
            warning(['if you want to simulate both theta and rho at once,' ...
                ' do it in two steps' ...
                ' 1) simulate theta with its parameters' ...
                ' 2) simulate rho with its parameters, but also use "simulationParams.theta = theta" to inherit the simulated theta'])
        end
        error('po_cfg.simulationParams.theta_rho must be either [0 1] or [1 0]')
end













