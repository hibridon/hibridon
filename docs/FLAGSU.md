To treat collisions of a molecule with a surface set FLAGSU = .TRUE.

In this case the squares of both the off-diagonal as well as the diagonal elements of the T matrix are equal to the squares of the corresponding S matrix elements. `FLAGSU` should be set .FALSE. for all gas-phase problems.

For collisions of molecules with surfaces, the orbital angular momentum is zero, and does not couple with the molecular rotational angular momentum. Thus, if FLAGSU = .TRUE., you must also set [CSFLAG](CSFLAG) = .TRUE., [JTOT1](JTOT) = 0, and [JTOT2](JTOT) = 0. Otherwise the program will abort with a diagnostic message.

For a good reference to the theory of inelastic collisions of a molecule with a flat, rigid surface see:

M. H. Alexander, J. Chem. Phys. **80**, 3485 (1984); G. C. Corey, J. E. Smedley, M. H. Alexander, and W.-K. Liu, Surf. Sci. **191**, 203 (1987), and references contained therein.