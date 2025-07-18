In a typical scattering calculation, a large channel basis is needed only at small values of Jtot. As Jtot increases, the centrifugal barrier prevents the collision partners from penetrating into the strong coupling region, so that the effective distortion of the asymptotic Hamiltonian is much less. Considerable time can be saved by dropping channels as Jtot increases. 

This is done in the Hibridon code by means of the parameter RCUT, as follows:

In the subroutine basis, the centrifugal barrier in each channel, evaluated at R = RCUT,

$h^2 l_i(l_i+1) / (2 \mu {\rm RCUT}^2)$


where $l_i$ is the orbital angular momentum of the ith channel, is compared to the internal energy of the channel, i. If the centrifugal barrier at $R$ = RCUT is greater than i, then the channel is still closed at $R$ = RCUT. If any of the channels which are open asymptotically are found to be closed at R = RCUT, then these channels as well as all other closed channels are dropped from the channel basis.

We have found this to be an effective and simple way of dropping channels at large Jtot. The rapidity by which the channels are dropped is controlled by the magnitude of RCUT. Obviously, if RCUT is taken equal to [RENDAI](RSTART,-RENDLD,-and-RENDAI) - the ending point of the integration - then no channels will ever be dropped. A good rule of thumb is to chose initially a small value for RCUT in the range RENDAI/3 to RENDAI/2 and then increase RCUT only if necessary.

⚠️  In cases where the collision energy lies near (either above or below) the threshold for one (or more) channels (in other words near the point where a particular asymptotic state just becomes energetically allowed), accuracy may be lost be dropping channels. In this case, set RCUT negative, which will block the calculation from dropping any channels. **By default it is probably best to set RCUT negative**.