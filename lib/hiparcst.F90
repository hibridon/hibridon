! hibridon parameter constants
module mod_hiparcst

integer, parameter :: IPAR_COUNT = 10  ! number of system independent parameters of type integer
integer, parameter :: RPAR_COUNT = 9   ! number of system independent parameters of type real
integer, parameter :: LPAR_COUNT = 28  ! number of system independent parameters of type logical

end module mod_hiparcst

! fcod parameters
! note: we use a specific module for each enum because it's not practical for the user to use the recommanded keyword only for these enum (the user would have to list them all in some cases, which makes code maintaining tedious)
module fcod_enum
enum, bind( C )
   enumerator ::  &
      FCOD_AIRYFL = 1, &
      FCOD_BASTST = 2, &
      FCOD_BATCH = 3, &
      FCOD_CHLIST = 4, &
      FCOD_CSFLAG = 5, &
      FCOD_FLAGHF = 6, &
      FCOD_FLAGSU = 7, &
      FCOD_IHOMO = 8, &
      FCOD_IPOS = 9, &
      FCOD_LOGDFL = 10, &
      FCOD_NOPRIN = 11, &
      FCOD_NUCROS = 12, &
      FCOD_PHOTOF = 13, &
      FCOD_PRAIRY = 14, &
      FCOD_PRLOGD = 15, &
      FCOD_PRPART = 16, &
      FCOD_PRSMAT = 17, &
      FCOD_PRT2 = 18, &
      FCOD_PRXSEC = 19, &
      FCOD_READPT = 20, &
      FCOD_RSFLAG = 21, &
      FCOD_T2TEST = 22, &
      FCOD_TWOMOL = 23, &
      FCOD_WAVEFL = 24, &
      FCOD_WRPART = 25, &
      FCOD_WRSMAT = 26, &
      FCOD_WRXSEC = 27, &
      FCOD_BOUNDC = 28
end enum
end module fcod_enum

! lpar parameters
! note: we use a specific module for each enum because it's not practical for the user to use the recommanded keyword only for these enum (the user would have to list them all in some cases, which makes code maintaining tedious)
module lpar_enum
enum, bind( C )
   enumerator ::  &
      LPAR_AIRYFL = 1, &
      LPAR_PRAIRY = 2, &
      LPAR_BASTST = 3, &
      LPAR_BATCH = 4, &
      LPAR_CHLIST = 5, &
      LPAR_CSFLAG = 6, &
      LPAR_FLAGHF = 7, &
      LPAR_FLAGSU = 8, &
      LPAR_IHOMO = 9, &
      LPAR_IPOS = 10, &
      LPAR_LOGDFL = 11, &
      LPAR_PRLOGD = 12, &
      LPAR_NOPRIN = 13, &
      LPAR_PRPART = 14, &
      LPAR_READPT = 15, &
      LPAR_RSFLAG = 16, &
      LPAR_PRSMAT = 17, &
      LPAR_T2TEST = 18, &
      LPAR_PRT2 = 19, &
      LPAR_TWOMOL = 20, &
      LPAR_WRSMAT = 21, &
      LPAR_WRPART = 22, &
      LPAR_WRXSEC = 23, &
      LPAR_PRXSEC = 24, &
      LPAR_NUCROS = 25, &
      LPAR_PHOTOF = 26, &
      LPAR_WAVEFL = 27, &
      LPAR_BOUNDC = 28
end enum
end module lpar_enum

! ipar parameters
! note: we use a specific module for each enum because it's not practical for the user to use the recommanded keyword only for these enum (the user would have to list them all in some cases, which makes code maintaining tedious)
module ipar_enum
   enum, bind( C )
   enumerator ::  &
      IPAR_JTOT1   = 1, &
      IPAR_JTOT2   = 2, &
      IPAR_JTOTD   = 3, &
      IPAR_JLPAR   = 4, &
      IPAR_NERG    = 5, &
      IPAR_NUMAX   = 6, &
      IPAR_NUMIN   = 7, &
      IPAR_NUD     = 8, &
      IPAR_LSCREEN = 9, &
      IPAR_IPRINT  = 10
   end enum
end module ipar_enum

! rpar parameters
! note: we use a specific module for each enum because it's not practical for the user to use the recommanded keyword only for these enum (the user would have to list them all in some cases, which makes code maintaining tedious)
module rpar_enum
   enum, bind( C )
   enumerator ::  &

      ! scattering mode parameters
      RPAR_SCAT_FSTFAC  = 1, &
      RPAR_SCAT_RINCR   = 2, &
      RPAR_SCAT_RCUT    = 3, &
      RPAR_SCAT_RENDAI  = 4, &
      RPAR_SCAT_RENDLD  = 5, &
      RPAR_SCAT_RSTART  = 6, &
      RPAR_SCAT_SPAC    = 7, &
      RPAR_SCAT_TOLAI   = 8, &

      ! bound state mode parameters
      RPAR_BOUND_R1      = 1, &
      RPAR_BOUND_R2      = 2, &
      RPAR_BOUND_C       = 3, &
      RPAR_BOUND_SPAC    = 4, &
      RPAR_BOUND_DELR    = 5, &
      RPAR_BOUND_HSIMP   = 6, &
      RPAR_BOUND_EIGMIN  = 7, &
      RPAR_BOUND_TOLAI   = 8, &

      ! parameters common to scattering and bound state mode
      RPAR_XMU     = 9

   end enum
end module rpar_enum
