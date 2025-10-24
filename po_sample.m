function [theta, r] = po_sample(po_cfg)

% DESCRIPTION:
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

po_cfg.sampling.distribution = 'vonMises';  % 'vonMises' | 'wrappedCauchy'
po_cfg.sampling.theta_mu = 
po_cfg.sampling.theta_kappa

end

%% input check and defaults

% sanity check: vonMises distribution requires parameters mu and kappa
% if simulating theta, we need mu and kappa for theta
% if simulating r, we need mu and kappa for r
% if simulating both, we need mu and kappa for each
% TO DO

% sanity check: wrappedCauchy distribution requires parameters mu and rho
% if simulating theta, we need mu and kappa for theta
% if simulating r, we need mu and kappa for r
% if simulating both, we need mu and kappa for each
% TO DO


%% shortcuts

nSamples  = po_cfg.nSamples;
kappa     = po_cfg.kappa;
mu        = po_cfg.mu;
theta_rho = po_cfg.theta_rho;
nModes = length(mu);

%% random seed


fieldLbl = fieldnames(po_cfg);
if ~any(strcmp(fieldLbl,'randomSeed'))
    randomSeed = 42;
    warning(['setting random seed for you: ' num2str(randomSeed)])
else
    randomSeed = po_cfg.randomSeed;
end

rng(randomSeed)
%% implementation
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
        pmf = r * (theta(2) - theta(1)); % approximate area per bin
        pmf = pmf / sum(pmf);              % make sure it sums to 1
        % step 3. draw theta samples from probability mass distribution
        sampleIdx = randsample(length(theta), nSamples, true, pmf);
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













