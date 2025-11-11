%% script to play (not formally validate) with Polar
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

%% add folder: Polar (toolbox)

% relative path to toolbox
targetFolder = fullfile('work','functions','Polar');
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
end

cd(addOnPath)

%% simulate values
% mode 1: simulate theta with fixed r

% simulation parameters
clear po_cfg
po_cfg.verbose = 1;
po_cfg.simulation.randomSeed = randi(100);
po_cfg.simulation.theta = true;
po_cfg.simulation.r     = false;
po_cfg.simulation.nSamples = 50;
po_cfg.simulation.distribution.family = 'wrappedCauchy';
po_cfg.simulation.distribution.theta.mu    = [pi];
po_cfg.simulation.distribution.theta.rho = [0.8];

% simulate data
[theta, r] = po_simulate(po_cfg);

% compute polar metrics (ITPC and its associated angle)
po_cfg.metrics.requests = ["meanResultantLength" "meanResultantAngle"];
metrics = po_meanResultant(theta,r,po_cfg);
disp(metrics)

% view simulation
figure(101); clf
po_cfg.viewParams.type = 'line';
po_view(theta,r,po_cfg)

%% simulate values
% mode 2: simulate r with uniform theta

% simulation parameters
clear po_cfg
po_cfg.verbose = 1;
po_cfg.simulation.randomSeed = 42;
po_cfg.simulation.theta = false;
po_cfg.simulation.r     = true;
po_cfg.simulation.r_noise = 0.1;
po_cfg.simulation.nSamples = 300;
po_cfg.simulation.distribution.family = 'wrappedCauchy';
po_cfg.simulation.distribution.r.mu    = [pi/2];
po_cfg.simulation.distribution.r.rho = [0.175];


% simulate data
[theta, r] = po_simulate(po_cfg);

% view simulation
figure(101); clf
po_cfg.viewParams.type = 'line';
po_view(theta,r,po_cfg)


% compute polar metrics (PAC)
% TO DO add dPAC and uPAC
po_cfg.metrics.requests = ["meanResultantLengthNorm"];
metrics = po_meanResultant(theta,r,po_cfg);
disp(metrics)

%% simulate values
% mode 3: simulate both r and theta

% simulation parameters
clear po_cfg
po_cfg.verbose = 1;
po_cfg.simulation.randomSeed = 42;
po_cfg.simulation.theta = true;
po_cfg.simulation.r     = true;
po_cfg.simulation.r_noise = 0;
po_cfg.simulation.nSamples = 200;
po_cfg.simulation.distribution.family = 'wrappedCauchy';
po_cfg.simulation.distribution.theta.mu    = [0];
po_cfg.simulation.distribution.theta.rho = [0];
po_cfg.simulation.distribution.r.mu    = [0];
po_cfg.simulation.distribution.r.rho = [0.6];

% simulate data
[theta, r] = po_simulate(po_cfg);

% decouple theta and r
po_cfg.decoupling.offsetRange = [10 90];
[theta_decoupled, r_decoupled ] = po_decouple(theta, r, po_cfg);
theta = theta_decoupled;
%r = r_decoupled;

figure(103); clf
plot(theta)

%%
% PAC
po_cfg.metrics.requests = ["meanResultantLength" "meanResultantLengthNorm" "UmeanResultantSquaredLengthNorm"];
metrics = po_meanResultant(theta, r, po_cfg);
disp(metrics)



% view simulation
figure(101); clf
po_cfg.viewParams.type = 'line';
po_view(theta,r,po_cfg)



%%
% debias the phase clustering
[theta_demeaned, r_demeaned ] = po_demean(theta, ones(size(theta)));

% dPAC
po_cfg.metrics.requests = ["meanResultantLength" "meanResultantLengthNorm" "UmeanResultantSquaredLengthNorm"];
metrics = po_meanResultant(theta_demeaned, r.*r_demeaned, po_cfg);
disp(metrics)

% view debiased data
figure(102); clf
po_cfg.viewParams.type = 'line';
po_view(theta_demeaned,r.*r_demeaned,po_cfg)

