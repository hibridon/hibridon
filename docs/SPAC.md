In scattering or photodissociation calculations (flag [BOUNDC](BOUNDC) = .false.) SPAC is the interval (sector) width in the log-derivative integration. The range from [RSTART](RSTART,-RENDLD,-and-RENDAI) to [RENDLD](RSTART,-RENDLD,-and-RENDAI) is divided into N steps, where
```
N = INT[(RENDLD - RSTART)/SPAC] + 1
````

Propagation is done using the Manolopoulos's extension [D. E. Manolopoulos, J. Chem. Phys., **85**, 6425 (1986)] of the original log-derivative propagator of Johnson. Typically SPAC should be on the order of 0.1 - 0.25 times the deBroglie wavelength for the system, which can be determined by the [command](Commands) [DEBROGLI](DEBROGLI) The sector width SPAC scales as the square root of the collision energy.

In [bound-state](Bound-States) calculations (flag [BOUNDC](BOUNDC) = .true.), SPAC is the spacing between the equally distributed Gaussian basis
 
$\chi_m(R) = e^{-\alpha(R-R_m)^2}$

The scale factor  is determined by means of the additional parameter [C](BOUNDC) as follows:

$\alpha = \frac{c^2}{SPAC^2}$