function [rho] = po_vonMisesDensity(kappa,mu,theta)

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
nModes = length(mu);
nSamples = length(theta);

%% implementation
% von Mises density distribution (Berens book p.16, Pewsey book p.56)
% mean of (potentially multiple) von Mises distribution(s)

vMdensity_perMode = nan(nModes,nSamples);
for mIdx = 1:nModes
    
    vMdensity_perMode(mIdx,:) = (1/(2*pi*besseli(0,kappa(mIdx)))) * exp(kappa(mIdx) * cos(theta - mu(mIdx)));
end
vMdensity = mean(vMdensity_perMode,1);

%% sanity checks
% sanity check: valid density distribution
% non negative everywhere
if any(vMdensity<0)
    error('likely coding bug: there aer some negative values in the density distribution')
end
% its integral is 1
distrIntegral = trapz(theta,vMdensity);
if abs(distrIntegral-1) > 1e-4
    %             warning([' the integral of the function is ' num2str(distrIntegral) ', just short of 1. Increase the number of observations (currently: '  num2str(nSamples) ') to improve the distribution'])
end

%% add noise
% --- Add noise to density values
%         noise_level = 0; % noise amplitude, adjust as needed
%         vMdensity_noisy = vMdensity + noise_level*randn(size(vMdensity));
%         % Ensure non-negativity (optional, for density)
%         vMdensity_noisy(vMdensity_noisy < 0) = 0;
%         % (Optional) Renormalize so area = 1
%         vMdensity_noisy = vMdensity_noisy / trapz(theta, vMdensity_noisy);
%         vMdensity = vMdensity_noisy;
        % ---
        
%%
rho = vMdensity(:);













