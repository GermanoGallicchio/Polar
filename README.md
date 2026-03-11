# _Fase_ for MATLAB

Fase is a toolbox for circular statistics. It includes functions for calculating angular and radial metrics used, for example, for Inter-Trial Phase Clustering, Phase Locking, Phase Amplitude Coupling. "Fase" (pronounced /ˈfaːze/ or "FAH-zeh") is the Italian for "phase"  

Author: <br>
Germano Gallicchio, [Bangor University](https://www.bangor.ac.uk/)

## How to cite

Gallicchio, G. (2025). Fase. GitHub. https://github.com/GermanoGallicchio/Fase/



## Overview / What Fase can do

- Compute mean resultant vector metrics (e.g., used for Inter-Trial Phase Clustering, Phase Locking, Phase Amplitude Coupling)
- Optional trimming of outliers
- Perform transformations (e.g., logit, FisherZ)
- Monte Carlo testing of metrics against null models via fa_test + fa_simulate

Assume we have a set of complex numbers c, each describing:
- a direction (or angle or phase) arg(c)
- a length (or magnitude) |c|
  
In formula (with [Euler's notation](https://en.wikipedia.org/wiki/Euler%27s_formula)):

$$
c = |c| \ e^{i \ arg(c)}
$$

The complex number could be made of angle and magnitude from the same signal (e.g., coefficients of continuous wavelent transform). Or the complex number could be artificially assembled with the angle of a signal and the magnitude of another signal (e.g., to study the phase-magnitude coupling). Sometimes we are only interested in magnitude (you won't need Fase for this). Other times we are only interested in angles (e.g., preference for a certain angle across replicates). Below is a series of formulas describing important polar metrics. By assembling angle and magnitude appropriately, **these polar metrics correspond with popular measures of angle-based consistency, such as Inter-Trial Phase Clustering (ITPC), phase-amplitude coupling (PAC), and others**. A non-exhaustive list will be provided below on how to use Fase metrics to obtain these measures.

### Length of mean resultant
$$
n^{-1} \left| \sum_{t=1}^n c_t \right|
$$

### Normalized length of mean resultant
$$
\frac{n^{-1} \left| \sum_{t=1}^n c_t \right|}{n^{-1} \sum_{t=1}^n \left| c_t \right|}
$$

### Angle of mean resultant
$$
\arg\left( \sum_{t=1}^n c_t \right)
$$

### U-statistic estimateor of the squared length of the mean resultant
upcoming...

### Normalized U-statistic estimateor of the squared length of the mean resultant
upcoming...






## Documentation, installation, and tutorials

Documentation for Fase will be available at this (link)[https://germanogallicchio.github.io/Fase_documentation/index.html]

Code developed on MATLAB R2025b on a Linux OS (Kubuntu).

## Tutorials

## How to...

### ...compute Inter-Trial Phase Clustering (aka, PLV) or Pairwise Phase Consistency

theta contains the phases of signal A

rho contains the magnitudes of signal A ... we don't care about it
```matlab
fa_meanResultant(theta, ones(size(theta)))
```
ask for "meanResultantLength" for ITPC and "UmeanResultantSquaredLength" for PPC

### ...compute Phase Amplitude Coupling

theta contains the phases of signal A

rho contains the magnitudes of signal B

we just need the (normalized, if you prefer) mean resultant length of rho exp(1i*theta)
```matlab
fa_meanResultant(theta, rho)
```

### ...compute debiased Phase Amplitude Coupling

theta contains the phases of signal A

rho contains the magnitudes of signal B
```matlab
[theta_demeaned, rho_demeaned] = fa_demean(theta,ones(size(theta)));
fa_meanResultant(theta_demeaned, rho.*rho_demeaned);
```

### ...test a metric against a null model (Monte Carlo)

Test whether an observed mean resultant length exceeds what you'd expect under a uniform-phase null (here using von Mises with kappa=0):

```matlab
% observed data
theta_obs = yourAngles(:);           % radians, N x 1
rho_obs   = ones(size(theta_obs));   % unit lengths for phase consistency

% which metrics to compute
metrics_cfg.requests = ["meanResultantLength"];  % or add others

% null model (uniform phases)
sim_cfg.sampling.theta = true;    % sample theta from null
sim_cfg.sampling.r     = false;   % r = 1
sim_cfg.sampling.distribution.family = 'vonMises';
sim_cfg.sampling.distribution.theta.kappa = 0;   % uniform
sim_cfg.sampling.distribution.theta.mu    = 0;   % irrelevant when kappa=0

% test options
test_cfg.nIter      = 2000;
test_cfg.alpha      = 0.05;
test_cfg.tail       = "greater"; % one-sided
test_cfg.randomSeed = 12345;

% run test
res = fa_test(theta_obs, rho_obs, metrics_cfg, sim_cfg, test_cfg);
disp(res.p.meanResultantLength)
```





