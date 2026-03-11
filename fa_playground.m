%% script to play (not formally validate) with Fase
%
% DESCRIPTIOM
%
%   A simple script for plug and play values to
%   - simulate values
%   - test whether they differ from a null distribution (Monte Carlo)
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% initial setup

close all; clearvars; clear global; clc;

%% add folder: Fase (toolbox)

% relative path to toolbox
targetFolder = fullfile('work','functions','Fase');
found = false;

if ispc % --- Windows, loop over drive letters A: to Z: ---
    letterVector = 65:90;   % ASCII codes for A–Z
    for diskLbl = string(char(letterVector) + ":")'
        addOnPath = fullfile(diskLbl, targetFolder);
        if exist(addOnPath,'dir')
            addpath(addOnPath)
            disp(['added to path: ' addOnPath])
            found = true;
            break
        end
    end

elseif isunix % --- Linux, look under /media/<user>
    mountFolder = fullfile('/media', getenv('USER'));
    % check one level down for drives, then into work/functions/PhysioExplorer
    matches = dir(fullfile(mountFolder, '*', targetFolder));
    if ~isempty(matches)
        addOnPath = matches(1).folder;
        addpath(addOnPath)
        disp(['added to path: ' addOnPath])
        found = true;
    end
        
end

if ~found
    warning([ targetFolder ' not added to path: folder not found'])
    error('Fase folder not found. Cannot run fa_playground.m without toolbox on path.')
end

cd(addOnPath)

%% simulate values
% mode 1: simulate theta with fixed r

% simulation parameters
clear fa_cfg
fa_cfg.verbose = true;
fa_cfg.simulation.randomSeed = 95;
fa_cfg.simulation.theta = true;
fa_cfg.simulation.r     = false;
fa_cfg.simulation.nSamples = 200;
fa_cfg.simulation.distribution.family = 'wrappedCauchy';
fa_cfg.simulation.distribution.theta.mu    = [pi 0];
fa_cfg.simulation.distribution.theta.rho = [0.999 0.999];

% simulate data
[theta, r] = fa_simulate(fa_cfg);

% compute polar metrics (ITPC and its associated angle)
fa_cfg.metrics.requests = ["meanResultantLength" "meanResultantAngle"];
metrics = fa_meanResultant(theta,r,fa_cfg);
disp(metrics)

% view simulation
figure(101); clf
fa_cfg.viewParams.type = 'line';
fa_view(theta,r,fa_cfg)

%% simulate values
% mode 2: simulate r with uniform theta

% simulation parameters
clear fa_cfg
fa_cfg.verbose = true;
fa_cfg.simulation.randomSeed = 42;
fa_cfg.simulation.theta = false;
fa_cfg.simulation.r     = true;
fa_cfg.simulation.r_noise = 0.7;
fa_cfg.simulation.nSamples = 300;
fa_cfg.simulation.distribution.family = 'wrappedCauchy';
fa_cfg.simulation.distribution.r.mu    = [-pi/2];
fa_cfg.simulation.distribution.r.rho = [0.05];


% simulate data
[theta, r] = fa_simulate(fa_cfg);

% this allows to study phase-amplitude avoidance
% normMinMax = @(x,y) ( (max(y)-min(y))*((x-min(x))/(max(x)-min(x)))+min(y) );  % normalize distribution x to have same min and max as distribution y
%r = normMinMax(-r,r); 

% view simulation
figure(101); clf
fa_cfg.viewParams.type = 'point';
fa_view(theta,r,fa_cfg)


% compute polar metrics (PAC)
% TO DO add dPAC and uPAC
fa_cfg.metrics.requests = ["meanResultantLengthNorm"];
metrics = fa_meanResultant(theta,r,fa_cfg);
disp(metrics)

%% simulate values
% mode 3: simulate both r and theta

% simulation parameters
clear fa_cfg
fa_cfg.verbose = true;
fa_cfg.simulation.randomSeed = 42;
fa_cfg.simulation.theta = true;
fa_cfg.simulation.r     = true;
fa_cfg.simulation.r_noise = 0;
fa_cfg.simulation.nSamples = 200;
fa_cfg.simulation.distribution.family = 'wrappedCauchy';
fa_cfg.simulation.distribution.theta.mu    = [0];
fa_cfg.simulation.distribution.theta.rho = [0];
fa_cfg.simulation.distribution.r.mu    = [0];
fa_cfg.simulation.distribution.r.rho = [0.6];

% simulate data
[theta, r] = fa_simulate(fa_cfg);

% decouple theta and r
fa_cfg.decoupling.offsetRange = [10 90];
fa_cfg.decoupling.randomSeed = 42;
fa_cfg.decoupling.perColumn = true;
[theta_decoupled, r_decoupled ] = fa_decouple(theta, r, fa_cfg);
theta = theta_decoupled;
%r = r_decoupled;

figure(103); clf
plot(theta)

%%
% PAC
fa_cfg.metrics.requests = ["meanResultantLength" "meanResultantLengthNorm" "UmeanResultantSquaredLengthNorm"];
metrics = fa_meanResultant(theta, r, fa_cfg);
disp(metrics)



% view simulation
figure(101); clf
fa_cfg.viewParams.type = 'line';
fa_view(theta,r,fa_cfg)



%%
% debias the phase clustering
[theta_demeaned, r_demeaned ] = fa_demean(theta, ones(size(theta)));

% dPAC
fa_cfg.metrics.requests = ["meanResultantLength" "meanResultantLengthNorm" "UmeanResultantSquaredLengthNorm"];
metrics = fa_meanResultant(theta_demeaned, r.*r_demeaned, fa_cfg);
disp(metrics)

% view debiased data
figure(102); clf
fa_cfg.viewParams.type = 'line';
fa_view(theta_demeaned,r.*r_demeaned,fa_cfg)

