If WAVEFL = .TRUE., then the transformation matrices are computed for the subsequent determination of the scattering wavefunction. These are written (in binary) to the file {jobname}.wfu

Subsequent calculation of the wavefunction is invoked by the command [PSI](PSI). Subsequent calculation of the incoming and outgoing fluxes is controlled by the command [FLUX](FLUX). 

⚠️ The total storage required to store the transformation matrix is

$(N_{log} + 2 N_{airy}) \times N_{ch}^2 + 3 N_{airy} \times N_{ch}$,

where N_{log} and N_{airy} are the number of LODG and AIRY integration steps and Nch is the number of channels. Thus the file {jobname}.wfu can often be very large.