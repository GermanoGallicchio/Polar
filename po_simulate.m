function [theta, rho] = po_simulate(po_cfg)

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

nSamples  = po_cfg.simulationParams.nSamples;
kappa     = po_cfg.simulationParams.kappa;
mu        = po_cfg.simulationParams.mu;
theta_rho = po_cfg.simulationParams.theta_rho;
nModes = length(mu);

switch num2str(theta_rho)
    case num2str([0 1])  % simulate rho with uniform theta
        % build an empirical density distribution by drawing from a
        % theoretical distribution (sum of von Mises).
    
        % empirical density distribution 
        thetaLim = [-pi, pi];  % create a range of phase values in (-pi, pi]
        theta = linspace(thetaLim(1), thetaLim(2), nSamples+1);
        theta(1) = [];
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
        if strcmp(num2str(po_cfg.simulationParams.theta_rho),num2str([1 1]))  % simulate both theta and rho
            warning('simulate theta and rho one at the time with same or different distributional parameters')
        end
        error('po_cfg.simulationParams.theta_rho must be either [0 1] or [1 0]')
end













