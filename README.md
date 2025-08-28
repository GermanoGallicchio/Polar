# _Polar_ for MATLAB
A circular statistics toolbox for computing polar metrics. It includes functions for calculating angular and radial metrics.

Author: <br>
Germano Gallicchio, [Bangor University](https://www.bangor.ac.uk/)

Reference: <br>
Gallicchio, G. (2025). Polar. GitHub. https://github.com/GermanoGallicchio/Polar/

Code developed on MATLAB R2022b on a Windows 10 device.

## Features
- Compute mean resultant vector metrics (e.g., used for Inter-Trial Phase Clustering, Phase Locking, Phase Amplitude Coupling)
- Optional trimming of outliers
- Perform transformations (e.g., logit, FisherZ)


## What Polar can do
Assume we have a set of complex numbers c, each describing:
- a direction (or angle or phase) arg(c)
- a length (or magnitude) |c|
  
In formula (with [Euler's notation](https://en.wikipedia.org/wiki/Euler%27s_formula)):

$$
c = |c| \ e^{i \ arg(c)}
$$

The complex number could be made of angle and magnitude from the same signal (e.g., coefficients of continuous wavelent transform). Or the complex number could be artificially assembled with the angle of a signal and the magnitude of another signal (e.g., to study the phase-magnitude coupling). Sometimes we are only interested in magnitude (you won't need Polar for this). Other times we are only interested in angles (e.g., preference for a certain angle across replicates). Below is a series of formulas describing important polar metrics. By assembling angle and magnitude appropriately, **these polar metrics correspond with popular measures of angle-based consistency, such as Inter-Trial Phase Clustering (ITPC), phase-amplitude coupling (PAC), and others**. A non-exhaustive list will be provided below on how to use Polar metrics to obtain these measures.

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


## How to...

### ...compute Inter-Trial Phase Clustering (aka, PLV)
theta contains the phases of signal A
rho contains the magnitudes of signal A ... we don't care about it
```matlab
po_meanResultant(theta, ones(size(theta)))
```

### ...compute Phase Amplitude Coupling
theta contains the phases of signal A
rho contains the magnitudes of signal B
we just need the (normalized in you prefer) mean resultant length of rho exp(1i*theta)
```matlab
po_meanResultant(theta, rho)
```

### ...compute debiased Phase Amplitude Coupling
theta contains the phases of signal A
rho contains the magnitudes of signal B
```matlab
[theta_demeaned, rho_demeaned] = po_demean(theta,ones(size(theta)));
po_meanResultant(theta_demeaned, rho.*rho_demeaned);
```
