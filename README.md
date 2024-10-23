# snapWavyForces 

A set of `python` scripts for assesing the stability of a droplet. The details of the contents of this repository is based on the paper by Bisquert, _et al._ (2024) [doi:...]().

The `python` scripts solve a system of differential equations that describe the shape of the interface of a 2D droplet sitting on a smooth sinusoidal surface of amplitude $a$. The interface is subject to Laplace's pressure jump at the interface which set the balance of forces. Two forces act on the interface:
 - The capillary force $\gamma \kappa$, where $\gamma$ is the surface tension and $\kappa$ is the curvature of the interface
 - The body force $\rho g_x x$, where $\rho$ is the density of the droplet and $g_x$ is the lateral projection of the gravitational acceleration. 

The system of differential equations is, in dimensionless units,
$ \frac{\mathrm{d}\phi}{\mathrm{d}s} = k0 + k1 x, $
$ \frac{\mathrm{d}x}{\mathrm{d}s} = \cos \phi, $
$ \frac{\mathrm{d}y}{\mathrm{d}s} = \sin \phi. $
The system is solved using the `solve_bvp` from the `scipy.integrate` package. As the variable $s$ ranges from 0 to $S$, and $S$ is uknown, a change of variables is necessary, thus using $\zeta = s/S$ and now $\zeta \in [0, 1]$ keeping $S$ as one of the parameters to be found by `solve_bvp`.

## Getting started

### Requisits

 - `scipy` a version that contains `solve_bvp`.
 - `numpy`.
 - `matplotlib` (optionally) for plotting results from the example.

### Installation

No need for installation. The `snapFixedRadius` class can be imported and run from another python script.

## Documentation

The script it documented by comments within.


> [!IMPORTANT] 
> Although it has been tested, this software comes with no warranty.
> Tested for `python-3.12.4` and the modules `numpy-1.26.4`, `scipy-1.13.1`, and `matplotlib-3.9.2` that come in [Anaconda3](https://www.anaconda.com/download/success).
