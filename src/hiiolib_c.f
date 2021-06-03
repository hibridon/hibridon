*********************************************************************
*                                                                   *
*                        i/o routines library                       *
*                                                                   *
*********************************************************************
*                         routines included:                        *
*                                                                   *
*   3. fimove   moves a sequential file to its end                  *
*   4. fimovs   moves s-matrix file to a specific record            *
*   6. gendat (sysdat) reads (saves) input for hibridon             *
*   7. io       block data, contains machine dependent settings     *
*   8. openf    utility to open files                               *
*   9. openfi   selects and opens all files for a run               *
*  10. restrt (rsave)   restart utility (useful when system crash)  *
*  11. sread (readrc/rdhead/inbfi/saddr) buffered s-matrix i/o      *
*  12. swrite (nexrec/wrhead/outbfi/clear) buffered s-matrix i/o    *
*  13. assgn    assigns and opens direct access files               *
*  14. rdabsf (wrabsf/fwait) read/write absolute routines           *
*  14. tmpnm    generates a unique filename                         *
*                                                                   *
*********************************************************************
*  NB cstart ultrix-dec for i/o with fortran instead of c routines
* ---------------------------------------------------------------
      subroutine fimove (nxfile)
* ---------------------------------------------------------------
*  subroutine to move sequential file on unit nxfile to its end
*  author:  millard alexander
*  current revision date: 23-sept-87
* ---------------------------------------------------------------
*    variables in call list:
*      nxfile:     logical unit
* ---------------------------------------------------------------
      character*10 line
      integer nxfile
10    read (nxfile, 15, err=20, end= 20) line
15    format (a)
      go to 10
20    return
      end
* -------------------------------------------------------------------
       subroutine fimovs(nfile,jtot,jlpar,nu,ien,ierr)
       implicit double precision (a-h,o-z)
*
* subroutine to move s-matrix file to jtot,jlpar partial wave
*
* note ! the file pointer is set to jtot-1, jlpar, which means that
*        the next "call writes" or "call nexrec" statement will over-
*        write the jtot record and not the jtot+1 record.
*
* author: b.follmeg
* current revision date: 4-oct-87
*
* -------------------------------------------------------------------
       common/cofil/ nfl,iofbuf,maxrec(60),iofrec(60),nwrec
       common /coener/ energ(1)
       ierr=0
       iaddr=0
       maxrec(nfile)=0
5      format(/' CHECKING S-MATRIX FILE FOR ENERGY',f10.3/)
*  read header block
       call  readrc(iaddr,nfile,lrec,jtot2,jlpar2,nu2,nopen,
     :              length,nnout)
       write(6,5) energ(ien)
10     call  readrc(iaddr,nfile,lrec,jtot2,jlpar2,nu2,nopen,
     :              length,nnout)
       if(lrec.lt.0) goto 20
       maxrec(nfile)=maxrec(nfile)+1
       write(6,15) jtot2,jlpar2,nu2
15     format(' READ JTOT=',i3,'  JLPAR=',i2,'  NU=',i3)
       iaddr=iaddr+lrec
       if(jtot2.ne.jtot .or. jlpar2.ne.jlpar.or. nu2.ne.nu) goto 10
* here if match
       iofrec(nfile)=iaddr
       return
* here if no match and eof
20     ierr=1
       return
       end
* ---------------------------------------------------------------
       subroutine gendat
*  subroutine to read system independent input parameters for hibridon code
*
*  author:  millard alexander
*  modifications by b. follmeg, h-j werner
*  current revision date:  10-jun-2006 by mha
* ------------------------------------------------------------------
*  line 1:
*    label:       label to uniquely specify your calculation (up to
*                 40 characters in length)
*  line 2:        basistype (integer 1-99)
*  line 3:
*    logdfl:      if .true., then logd propagation will take place
*    airyfl:      if .true., then airy propagation will take place
*    readpt:      if .true., then potential parameters are expected
*    bastst:      if .true., then execution terminates after the first call
*                 to basis
*  line 4:
*    rstart:      starting point (bohr) for integration
*                 at subsequent total angular momentum values rstart is
*                 adjusted to remain a constant distance inside the innermost
*                 classical turning point
*    rendld:      ending point (bohr) of log-derivative integration
*                 at subsequent total angular momentum values rendld is
*                 adjusted so that the range covered by the log-derivative
*                 integration remains constant
*    spac:        step size in the logd integration *    if logdfl=.false., 1.
*                          equal to rstart by this subroutine
*                       2. then spac will be the value of the first interval
*                          width in the airy propagation
*  line 5:
*      prairy:   if .true., then step-by-step information is printed out in
*                airy propagation
*  line 6:
*      tolai:    error parameter to determine subesequent step sizes
*      rincr:    power at which step sizes can increase
*      rendai:   ending point of airy integration
*      fstfac:   factor by which initial step size in logd integration is
*                multiplied to get initial step size in airy integration
*                if airyfl is .false., then the data in lines 4 and 5 will be
*                ignored!
*  line 7:
*    nerg:       the number of total energies at which calculation is to be
*                performed
*  line 8:
*    energ(i):   these energies in wavenumbers
*                if nerg.gt.1, then:
*                logd portion of integration involves writing quadrature
*                matrices onto file 11 at first energy and reading them back
*                in at subsequent energies
*                airy portion of integration involves writing transformation
*                matrices onto file 10 at first energy and reading them back
*                in at subsequent energies
*  line 9:
*    xmu:        collision reduced mass in carbon-12 amu
*                for a collision between atom a and molecule bc the collision
*                reduced mass is
*                                xmu = m  m    / (m  + m  )
*                                       a  bc      a    bc
*  line 11:
*    rcut:       cut-off point for keeping higher energy channels
*                if any open channel is still closed at r=rcut, then
*                all closed channels as well any open channels which
*                are still closed at r=rcut are dropped from basis
*  line 11:
*    jtot1:     starting value of total angular momentum (l-bar in coupled
*               states calculation
*    jtot2:     ending value of total angular momentum
*    jtotd:     step size for total angular momentum
*               this should be 1 for collisions with an uncorrugated surface
*    jlpar:     a switch to restrict which parity channels are included in a c
*               calculation.  only those channels for which
*                     eps * (-1)**(j + l -jtot) = jlpar
*               are included
*               in a cs calculation, this parameter can have any value on inpu
*               it is always set equal to +1
*    numin:     minimum coupled-states projection index
*    numax:     maximum coupled-states projection index
*               for molecules with even multiplicity (half-integer spin) the
*               true half-integer values of the total angular momentum and
*               the cs projection index are equal, respectively, to jtot + 1/2
*               and nu + 1/2
*               in the case of a cc calculation, numin and numax can be given
*               any value
*    nud:       step size for nu (nu=numin:nud:numax)
*  line 12:
*    nnout:     the values of the channel rotational quantum numbers for which
*               the s-matrix elements will be printed to file 91 (if writs = t
*    niout:     the values of the additional channel index for which
*               the s-matrix elements will be printed to file 91 (if writs = t
*  line 13:
*    jout:      an array containing these values of the rotational angular
*               momenta
*  line 14:
*    indout:    an array containing these values of the additional
*               channel index
*  line 15:
*   logical variables (input format l3) to control output:
*     prlogd:  if .true., then the lower triangle of log-derivative matrix
*              is printed at the end of the logd and the end of the airy
*              integrations
*     prsmat:  if .true., then the upper triangle of real and imaginary parts
*              of s-matrix are printed
*     prt2:    if .true., then the upper triangle of square modulus of t-matri
*              is printed
*     t2test:  if .true., then the first two columns of the square modulus of
*              the t-matrix are printed
*     wrsmat:  if .true., and nnout is > 0, then those s-matrix elements
*              for which both the initial and final rotational quantum
*              numbers are in the array jout (input line 12) are written to
*              files smat1, smat2, ...
*              if nnout < 0, then each column of the s-matrix whose initial
*              index is in the array jout is written to files smat1, smat2, ..
*  line 16:
*   logical variables (input format l3) to control output:
*     wrpart:  if .true., then input data and the matrix of partial cross
*              sections (summed over m-states) is written to file pxsec
*     prpart:  if .true, then the full matrix of partial cross sections (summe
*              over m-states) is printed
*     prxsec:  if .true., then the full matrix of integral cross sections
*              ((summed over m-states and summed from jtot1 to jtot2) is print
*     wrxsec:  if .true., then some input data and the full matrix of integral
*              cross sections ((summed over m-states and summed from
*              jtot1 to jtot2) is written to file xsec1, xsec2, ....
*     wavefl:  if .true. then information is written to calculate,
*              subsequently, wavefunctions, fluxes, and adiabatic energies
*     boundc:  if .true. then susan gregurick's bound state calculation
*              is implemented
*  line 17:
*   logical variables (input format l3) to control output:
*     noprin:  if .true., then most printing is suppressed
*     chlist:  if .true., then the channel quantum numbers and energies are
*              printed out at at each total-j
*              if .false., then  this is done only at first total-j
*     ipos:    if .true., then printout is suited for a 132-position printer
*              if .false., then printout is suited for a  80 -position printer
*     nucros:  parameter to control how CS integral cross sections are
*              computed
*     photof:  if .true. then photodissociation calculation
*  line 18:
*    flaghf:   if .true., then the system has even multiplicity (half-integer
*              total angular momentum)
*    csflag:   if .true., then coupled-states calculation is desired
*              if .false., then close-coupled calculation is desired
*    flagsu:   if .true., then the problem is assumed to a molecule scattering
*              of a surface, in which case the diagonal elements of the
*              transition probabilities are equal to the modulus squared of
*              the s-matrix (not t-matrix elements)
*    ihomo:    if .true., then the molecule is assumed to be homonuclear
*    twomol:   if .true., then molecule-molecule collision is assumed
*  line 19:
*    rsflag:      if .true., then calculation is to be restarted
*                 a check will be made to see if all required files
*                 are present:  these may include
*                    trstrt, tmp10, tmp11, xsecn (or tmpxn), smatn,
*                    psecn, tmp35, ...
*  variable in common block /coener/
*    energ:     array containing total energies (in cm-1) at which scattering
*               calculations are to be performed
*  variables in common block /cosout/
*    nnout:     number of different rotational levels for which s-matrix
*               elements are to be saved in files smat1, smat2, ...
*    jout(i):   values of rotational angular momentum for these lvels
*  variables in common block /cophot/
*     photof        true if photodissociation calculation
*                   false if scattering calculation
*     wavefn        true if g(a,b) transformation matrices are saved
*                   to be used later in computing the wavefunction
*  variable in common block /coselb/
*     ibasty    basistype
*
*
* subroutines called: open
*
* ----------------------------------------------------------------
      use mod_cosout, only: nnout, jout
      use mod_coiout, only: niout, indout
      implicit double precision (a-h,o-z)
      integer i, jlpar, jtot1, jtot2, jtotd, length, nerg,
     :        numax, numin, ibasty
      logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart,
     :        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf,
     :        csflag, flagsu, rsflag, t2test, existf, logdfl, batch,
     :        readpt, ihomo, bastst, twomol, nucros, photof, wavefl,
     :        boundc
      character*40 input, jobnam, output, savfil
      character*(*) filnam
      include "common/parpot"
      common /coselb/ ibasty
      common /coipar/ jtot1, jtot2, jtotd,jlpar,nerg,numax, numin, nud,
     :                lscreen, iprint
      common /corpar/ fstfac, rincr, rcut, rendai, rendld, rstart, spac,
     :                tolai, xmu
* NB if boundc = .true. then these variables are:
*      common /corpar/ r1,r2,c,spac,delr,hsimp,eigmin,tolai,xmu
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo, ipos, logdfl, logwr,
     :                noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr, nucros, photof, wavefl, boundc
      common /coskip/ nskip,iskip
      common /cofile/ input, output, jobnam, savfil
      common /coener/ energ(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv

* ----------------------------------------------------------------
*  open unit 8 for standard input
      length = index(input,' ') - 1
      inquire (file=input, exist=existf)
      if (.not. existf) then
        write (6, 20) input(1:length)
20      format
     :   (/'   *** INPUT FILE ',(a),' DOES NOT EXIST; ABORT ***')
        if( batch) call exit
        return
      end if
* open sequential/formatted (mode='sf') file
      call openf(8, input, 'sf', 0)
* ----------------------------------------------------------------
      rewind 8
      iline=1
*  read in input data
*  line 1
      read (8, 40, err=195) label
40    format(a40)
      iline = iline + 1
*  line 1a
      read (8, *, err=195) ibasty
      iline=iline+1
*  line 2
      read (8, 50, err=195) logdfl, airyfl, readpt, bastst
50    format (8l3)
      iline = iline + 1
*  line 3
      read (8, *, err=195) rstart, rendld, spac
      iline = iline + 1
*  set rendld = rstart if logdfl=.false.
      if (.not. logdfl) then
        rendld = rstart
      end if
*  line 4
        read (8, 50, err=195) airypr
      iline = iline + 1
*  line 5
        read (8, *, err=195) tolai, rincr, rendai, fstfac
      iline = iline + 1
      if (.not. airyfl) then
*  give dummy values to airy parameters if airyfl = .false.
        airypr = .false.
        tolai = 1.
        rincr = 1.
        rendai = rendld
      end if
*  line 6
      read (8, *, err=195) nerg
      iline = iline + 1
*  line 7
      read (8, *, err=195) (energ(i), i = 1, nerg)
      iline = iline + 1
*  line 8
      read (8, *, err=195) xmu
* convert to atomic units of mass and store in /coered/
      rmu=xmu/xmconv
      iline = iline + 1
*  line 9
      read (8, *, err=195) rcut
      iline = iline + 1
*  line 10
      read (8, *, err=195) jtot1,jtot2,jtotd,jlpar,numin,numax,nud
      iline = iline + 1
      read (8, *, err=195) lscreen, iprint
      iline = iline + 1
*  line 11
      read (8, *, err=195) nnout,niout
      iline = iline + 1
*  line 12
      read (8, *, err=195) (jout(i), i=1, iabs(nnout))
      iline = iline + 1
      if(niout.gt.0) then
*  line 13
        read (8, *, err=195) (indout(i), i=1, niout)
        iline = iline + 1
      end if
*  line 14
      read (8, 50, err=195) logwr, swrit, t2writ, t2test, writs
      iline = iline + 1
*  line 15
      read (8, 50, err=195) wrpart, partw, xsecwr, wrxsec, wavefl
      iline = iline + 1
*  line 16
      read (8, 50 ,err=195) noprin, chlist, ipos, nucros, photof
      iline = iline + 1
*  line 17
      read (8, 50, err=195) flaghf, csflag, flagsu, ihomo, twomol
      nskip=1
      if(ihomo) nskip=2
      iline = iline + 1
*  line 18
      read (8, 50, err=195) rsflag, boundc
*  open unit 9 for standard output under filename outpt
      call openf(9, output, 'sf', 0)
* ----------------------------------------------------------------
      entry genchk
* check that basistype is allowed
      call baschk(ibasty)
      if (.not. logdfl) then
*  force rendld to rstart if no logderivative propagation
         if (rendld.ne.rstart) then
           rendld=rstart
           write (6, 105) rendld
105        format ('*** RENDLD SET EQUAL TO RSTART =',f9.4,' ***')
         endif
      endif
*  if collision of a molecule with a surface, then abort unless
*  jtot1 = 0, jtot2 = jtot1, and csflag = .true.
      if (flagsu) then
        if (.not. csflag) then
          write (6, 110)
110       format
     :     (' *** CSFLAG = .FALSE. FOR SURFACE SCATTERING; ABORT ***')
          if (batch) call exit
        end if
        if (jtot1 .ne. 0) then
          write (6, 120) jtot1
120       format (' *** JTOT1 =',i4,
     :            ' .NE. 0 FOR SURFACE SCATTERING; ABORT ***')
          if (batch) call exit
        end if
        if (jtot2 .ne. 0) then
          write (6, 130) jtot2
130       format (' *** JTOT2 =',i4,
     :            ' .NE. 0 FOR SURFACE SCATTERING; ABORT ***')
          if (batch) call exit
        end if
*  set jtotd = 1 for surface collisions.  this guarantees that the
*  quantities accumulated as cross sections will correspond to the
*  degeneracy averaged transition probabilities
        if (jtotd .ne. 1) then
          write (6, 135) jtotd
135       format (' *** JTOTD =', i3,
     :            '; RESET TO 1 FOR SURFACE SCATTERING')
          jtotd = 1
        end if
      end if
*  if cc calculation set numin and numax to zero
* NB this is disabled for 2P+homonuclear
      if (.not. csflag.and.ibasty.ne.12) then
        numin = 0
        numax = 0
      end if
*  if cs calculation, check that numin and numax are both positive
*  and that numax .ge. numin
      if (csflag) then
        if (numin .lt. 0) then
          write (6, 140) numin
140       format (' *** NUMIN =',i3, ' .LT. 0; ABORT ***')
          if (batch) call exit
          return
        end if
        if (numax .lt. 0) then
          write (6, 145) numax
145       format (' *** NUMAX =',i3, ' .LT. 0; ABORT ***')
          if (batch) call exit
          return
        end if
        if (numin .gt. numax) then
          if (batch) then
             write (6, 150) numin, numax
150          format (' *** NUMIN =',i3,' .GT. NUMAX =',i3,'; ABORT ***')
             call exit
          else
             write (6, 151) numin, numax
151          format
     :        (' *** NUMIN =',i3,' .GT. NUMAX =',i3,'; NUMAX RESET ***')
             numax=numin
             return
          endif
        end if
      end if
*   if cs calculation, set jlpar to +1
      if (csflag) then
        jlpar = 1
      end if
       if (iabs(jlpar) .gt. 1) then
         write (6, 155) jlpar
 155     format (' *** JLPAR =',i3,' .NE. +/- 1; ABORT ***')
         if (batch) call exit
         return
       end if
*   if writs = .true., then warning unless jtotd = 1
      if (writs .and. jtotd .ne. 1) then
        write (6, 170) jtotd
170     format
     :   (' WARNING *** WRSMAT = .TRUE. BUT JTOTD=',i3,' .NE. 1 ')
      end if
* here if photodissociation calculation or wavefunction desired
* reset all flags accordingly
      if (photof .or. wavefl) then
         if (t2writ) then
          write (6, 171)
171       format (' PRT2 SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          t2writ = .false.
        endif
        if (t2test) then
          write (6, 172)
172       format (' T2TEST .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          t2test = .false.
        endif
       if (wrpart) then
          write (6, 175)
175       format (' WRPART SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          wrpart = .false.
        endif
        if (partw) then
          write (6, 176)
176       format (' PRPART SET .FALSE., SINCE PHOTFL OR WAVEFN .TRUE.')
          partw = .false.
        endif
        if (writs) then
          write (6, 178)
178       format (' WRSMAT SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          writs = .false.
        endif
        if (xsecwr) then
          write (6, 179)
179       format (' PRXSEC SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          xsecwr = .false.
        endif
        if (wrxsec) then
          write (6, 180)
180       format (' WRXSEC SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          wrxsec = .false.
        endif
        if (wrxsec) then
          write (6, 181)
181       format (' WRXSEC SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          wrxsec = .false.
        endif
        if (jtot2 .gt.jtot1) then
          write (6, 182)
182       format (' JTOT2 SET EQUAL TO JTOT1,',
     :              ' SINCE PHOTOF OR WAVEFN .TRUE.')
          jtot2=jtot1
        endif
        if (numax .gt.numin) then
          write (6, 183)
183       format (' NUMAX SET EQUAL TO NUMIN,',
     :              ' SINCE PHOTOF OR WAVEFN .TRUE.')
          numax=numin
        endif
        if (jlpar .eq. 0) then
          write (6, 184)
184       format (' JLPAR=0, RESET TO 1, SINCE PHOTOF OR WAVEFN .TRUE.')
          jlpar=1
        endif
        if (nucros) then
          write (6, 185)
185       format (' NUCROS SET .FALSE., SINCE PHOTOF OR WAVEFN .TRUE.')
          nucros=.false.
        endif
        if (nerg .gt. 1) then
          write (6, 186)
186       format(' NERG RESET TO 1, SINCE PHOTOF OR WAVEFN .TRUE.')
          nerg=1
        endif
      endif
      if (noprin .and. iprint .ne. -1) then
         write (6, 190)
190      format (' NOPRIN = .TRUE., SO IPRINT SET TO -1')
      endif
      if (jtotd .eq. 0) then
        write (6, 191)
191     format(' JTOTD = 0; SET EQUAL TO 1')
        jtotd=1
      endif
      if (nucros) then
        if (.not.wrxsec .and. .not.xsecwr .and. .not.partw
     :                  .and. .not.wrpart) then
          write (6, 192)
192       format(' NUCROSS SET .FALSE. BECAUSE PRPART, PRXSEC,',
     :           ' WRPART, AND WRPART ARE ALL FALSE')
          nucros = .false.
        endif
      endif
      return
* here if read error
195   write(6,200) iline, input(1:length)
200   format(/'   *** ERROR IN GENDAT READING LINE',i3,
     :        ' OF FILE ',(a),', ABORT ***')
      if(batch) call exit
      return
* ----------------------------------------------------------------
      entry savdat(inew, filnam)
*  entry here to save new system independent input parameters on unit 8
*  on input filnam -> file name of input file, if inew=1 this is a new
*  file, if inew=0, then existing input file is overwritten
      inquire(file=input,exist=existf)
      length=index(input,' ')
      if (inew .eq. 0 .and. existf) then
* here if input file already exists
          write (6, 202) input(1:length)
202        format(' *** INPUT FILE ',(a),'OVERWRITTEN')
      else
* create new input file
        if (inew .eq. 0) then
           write(6,200) input(1:length)
        else
           length=index(filnam,' ')
           inquire(file=filnam,exist=existf)
           if (existf) then
             write (6, 202) filnam(1:length)
           else
             write(6,205) filnam(1:length)
205          format('  *** NEW INPUT FILE ',(a),'CREATED')
           endif
        endif
      endif
      if (inew .eq. 0) then
         call openf(8, input, 'sf', 0)
      else
         close (8)
         call openf(8, filnam, 'sf', 0)
      endif
      rewind (8)
      nline=0
*  line 1
      nline=nline+1
      write (8, 210, err=999) label
210   format ((a))
*  line 1a
      nline=nline+1
      write (8, 215, err=999) ibasty
215   format(i4,25x,'    ibasty')
*  line 2
      nline=nline+1
      write (8, 220, err=999) logdfl, airyfl, readpt, bastst
220   format (4l3,18x, '   logdfl, airyfl, readpt, bastst')
*  line 3
      nline=nline+1
      if (.not.boundc) then
        write (8, 230, err=999) rstart, rendld, spac
230     format (3f10.4, '   rstart, rendld, spac')
      else
        write (8, 232, err=999) rstart, rendld, spac
232     format (2f10.4,g11.4, '  hsimp, delr, eigmin')
      endif
*  line 4
      nline=nline+1
      write (8, 240, err=999) airypr
240   format (l3, 27x,'   prairy')
*  line 5
      nline=nline+1
      if (.not.boundc) then
        write (8, 250, err=999) tolai, rincr, rendai, fstfac
250     format (g11.4, f6.2, f8.1, f6.2,
     :        '  tolai, rincr, rendai, fstfac')
      else
        write (8, 252, err=999) tolai, rincr, rendai, fstfac
252     format (g11.4, f6.2, f8.4, f6.3,
     :        '  tolai, r2, spac, r1')
      endif
*  line 6
      nline=nline+1
      write (8, 260, err=999) nerg
260   format (i4,26x, '   nerg')
*  line 7
      nline=nline+1
      write (8, 270, err=999) (energ(i), i = 1, nerg)
270   format (10f11.4)
*  line 8
      nline=nline+1
      write (8, 280, err=999) xmu
280   format (f11.5, 19x, '   xmu')
*  line 9
      nline=nline+1
      if (.not. boundc) then
        write (8, 290, err=999) rcut
290     format (f10.4, 20x,'   rcut')
      else
        write (8, 292, err=999) rcut
292     format (f10.4, 20x,'   c')
      endif
*  line 10
      nline=nline+1
      write (8, 300, err=999) jtot1, jtot2, jtotd, jlpar, numin,
     :                        numax, nud
300   format (3i4,4i4,2x,'   jtot1,jtot2,jtotd,jlpar,numin,numax,',
     :                       'nud')
      nline=nline+1
      write (8, 301, err=999) lscreen, iprint
301   format(2i4,23x,'  lscreen, iprint')
*  line 11
      nline=nline+1
      write (8, 305, err=999) nnout, niout
305   format (2i5,20x,'   nnout,niout')
*  line 12
      nline=nline+1
      if (twomol) then
       write (8,310,err=999)(jout(i)/10,mod(jout(i),10),i=1,iabs(nnout))
310    format (20('  ',2i1))
      else
        write (8, 315, err=999) (jout(i), i=1, iabs(nnout) )
315     format (20i4)
      end if
*  line 13
      if(niout.gt.0) then
        nline=nline+1
        write (8, 315, err=999) (indout(i), i=1, niout )
      endif
*  line 14
        nline=nline+1
      write (8, 330, err=999) logwr, swrit, t2writ, t2test, writs
330   format (5l3,15x,'   prlogd, prsmat, prt2, t2test, wrsmat')
*  line 15
        nline=nline+1
      write (8, 340, err=999) wrpart, partw, xsecwr, wrxsec, wavefl
340   format (5l3,15x,'   wrpart, prpart, prxsec, wrxsec, wavefl')
*  line 16
        nline=nline+1
      write (8, 350, err=999) noprin, chlist, ipos, nucros, photof
350   format (5l3,15x,'   noprin, chlist, ipos, nucros, photof')
*  line 17
        nline=nline+1
      write (8, 360, err=999) flaghf, csflag, flagsu, ihomo, twomol
360   format (5l3,15x,'   flaghf, csflag, flagsu, ihomo, twomol')
*  line 18
        nline=nline+1
      write (8, 370, err=999) rsflag, boundc
370   format (2l3,24x, '   rsflag, boundc')
      return
* here if write error
999   if (inew .eq. 0) then
        length = index(input,' ') - 1
        write(6, 380) input(1:length), nline
      else
        length = index(filnam,' ') - 1
        write(6, 380) filnam(1:length), nline
      endif
380   format(' *** ERROR WRITING ',(a),':  LINE ',i3)
      if(batch) call exit
      end
* ---------------------------------------------------------------
      block data io
* ---------------------------------------------------------------
* current revision date: 26-may-1991
*
*  variables in common block /clseg/
*    lseg:      number of integer words per disc sector
*    intrel:    number of integer words per real words
*    lchar:     number of characters per integer word
      common/clseg/ lseg, intrel, lchar
*  variables in common block /cobuf/
*    lbuf:    length of i/o buffer
*    ibuf:      i/o buffer
*  variable in common block /cosize/
*    isize:     size of files (only needed for univac, optional on vax)
*    isizes:    size of s-matrix file (only needed for univac, optional
*               on vax
* -----------------------------------------------------
      common /cosize/ isize,isizes
cstart unix
      common/cobuf/ lbuf,ibuf(1024)
      data lseg, intrel, lchar /1024,2,4/
      data lbuf/1024/
cend
cstart cray
c;      common/cobuf/ lbuf,ibuf(512)
c;      data lbuf/512/
c;      data lseg, intrel, lchar /512,1,8/
cend
      end
* ---------------------------------------------------------------
      subroutine openf(lunit,filnam,lmode,isize)
* --------------------------------------------------------
*
* subroutine to open files
* author: b. follmeg
* current revision date: 29-dec-2003 by mha
*
* on input: lunit  -> logical unit number, if lunit < 0 scratch file
*           filnam -> file name
*           lmode  -> 'SF' sequential/formatted
*                     'SU' sequential/unformatted
*                     'DU' direct/unformatted using assgn, wrabsf, rdabsf, fwa
*           isize  -> number of tracks (only for univac and vax)
*
*   lseg:  number of integer words per disc sector
      logical exstfl, openfl, tmpfil
      logical od
      character*12 fmt, stat, accs
      character*(*) filnam
      character*(*) lmode
      character*2  mode
      common /clseg/ lseg
      ierr=0
      iunit = lunit
      mode = lmode(1:2)
      call upper(mode)
* scratch file requested ?
      if(iunit.lt.0) then
         tmpfil=.true.
         iunit=iabs(iunit)
      else
         tmpfil=.false.
      end if
* all direct access i/o done via rdabsf, wrabsf, assgn, closf
      if(mode(1:1).eq.'D') then
         isiz=isize
         if(tmpfil) then
            call assgn(iunit,filnam,isiz,0)
         else
            call assgn(iunit,filnam,isiz,1)
         end if
         return
      end if
* select access
      if(mode(1:1).ne.'S') then
         write(6,10) mode
         stop
      end if
* select format
      if(mode(2:2).eq.'F') then
         fmt='formatted'
      else if(mode(2:2).eq.'U') then
         fmt='unformatted'
      else
         write(6,10) mode
         stop
      end if
* inquire file specifications
      inquire(file=filnam, exist=exstfl, opened=openfl)
      accs='sequential'
      if (exstfl) then
        if (openfl) return
        stat='old'
* make sure sequential formatted files are appended not overwritten
cstart unix-hp unix-dec unix-iris unix-sun unix-ifort
c;        accs='append'
cend
      else
        stat='new'
      end if
      if (tmpfil) then
cstart unix  mac
        stat = 'scratch'
cend
cstart cray
c;      stat = 'unknown'
cend
cstart unix-darwin
c;        inquire(unit=iunit,opened=od)
c;* if temporary file is already opened, close it
c;        if (od) close(unit=iunit)
cend
cstart unix mac cray
        open(unit=iunit,  access='sequential',
     :          form=fmt, status=stat, err=999, iostat=ierr)
cend
      else
*       print *,stat
        open(unit=iunit, file=filnam, access=accs,
     :          form=fmt, status=stat, err=999, iostat=ierr)
      end if
      return
10    format(' *** ERROR IN OPEN, UNKNOWN MODE=',a,' VALID MODES ARE',
     :       ' SF,SU,DU *** ABORT')
* here if open error
999   write(6,20) ierr,filnam,iunit,mode
20    format(' *** ERROR IN OPEN, IOSTAT=',i5,/,
     :       ' *** FILENAME=',a,/,
     :       ' *** IUNIT=',i4,' MODE=',a,/,
     :       ' *** PROGRAM WILL BE STOPPED !')
      stop
      end
*  --------------------------------------------------------------------
      subroutine openfi (nerg)
*  subroutine to open required i/o files for hibridon program
*  author:  millard alexander
*  modifications: bernd follmeg, g. v. s.
*  current revision date:  17-jan-2008 by mha
*  --------------------------------------------------------------------
*  variables in call list:
*    nerg:       number of different total energies at which scattering
*                calculation is to be done
*                if nerg.gt.1, then three files are opened:
*                              unit=10 (filename tmp10) for storage of
*                                      transformation matrices in airy
*                                      propagation
*                              unit=11 (filename tmp11) for storage of
*                                      quadrature matrices in logd propagation
*                              unit=12 (filename tmp12) for storage of
*                                      rotational angular momenta, orbital
*                                      angular momenta, extra quantum index,
*                                      and internal energies for all channels
*  variables in common block /colpar/  (see further description in subroutine
*                                       flow)
*    airyfl:     note, unit=10 is opened only if airyfl = .true.
*    writs:      if .true., then unit=91 to unit=(90+nerg) are opened as files
*                           smat1, smat2, ... smatnerg for
*                           storage of real and imaginary parts of
*                           selected elements of s-matrix
*    wrpart:     if .true., then unit=41 to unit=(40+nerg) are opened as files
*                           psec1, psec2,... psecnerg for
*                           storage of of some input date and degeneracy
*                           averaged partial cross sections
*    csflag:     if .true., and partw or or wrpart or xsecwr or wrxsec = .true
*                           then unit=66 to unit=(65+nerg)
*                           are opened as files tmp66, tmp67, ... etc.
*                           for accumulation of partial cross sections
*                           at each cs projection index
*    wrxsec, xsecwr:
*                if either of these variables is .true., then unit=16 to
*                           unit=(15+nerg) are opened for storage of some
*                           input data and degeneracy averaged integral
*                           cross sections
*                if wrxsec = .true., then unit=16 to unit=(15+nerg) are opened
*                                    as permanent files with filenames
*                                    xsec1, xsec2, xsec3, ... , xsecn
*                                    where n = nerg
*                          = .false., then unit=16 to unit=(15+nerg) are opene
*                                     as files tmpx1, tmpx2, ... tmpxn
*    rsflag:     if .true., then calculation is being restarted
*                abort will occur unless all requested i/o files already exist
*                if .false, then initial calculation
*                in any case unit=14 (filename trstrt) is opened to hold
*                temporary information in case of restart
*  variable in common block /cosize/
*    isize:     size of files (only needed for univac, optional on vax)
*    isizes:    size of s-matrix file (only needed for univac, optional
*               on vax
*
*  subroutines called: open
*  --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ifile, nerg, nfile, lenx, isize, isizes
      logical airyfl, airypr, bastst, batch, chlist, csflag, existf,
     :        flaghf, flagsu, ihomo, ipos, logdfl, logwr,
     :        noprin, partw, readpt, rsflag, swrit,
     :        t2test, t2writ, writs, wrpart, wrxsec,
     :        xsecwr, twomol, nucros, photof, wavefl, boundc
      character*40  oldlab,newlab
      include "common/parpot"
      character*40 xname,xnam1
      character*20 cdate
      character*40 input,output,jobnam,savfil
      common /cofile/ input,output,jobnam,savfil
      common /cosize/ isize, isizes
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos,
     :                logdfl, logwr, noprin, partw, readpt,
     :                rsflag, swrit, t2test, t2writ, twomol, writs,
     :                wrpart, wrxsec, xsecwr, nucros, photof, wavefl,
     :                boundc
      common /coisc1/ isc1(9)
      common /coisc2/ isc2(1)
      common /coisc3/ isc3(3)
      common /coisc4/  isc4(1)
      common /cosc1/  rsc1(2)
      common /cosc2/  rsc2(1)
      common /colsc1/  lsc1(5)
      common /coselb/ ibasty
      if (nerg .gt. 1) then
*  check to see if nerg .le. 25
        if (nerg .gt. 25) then
          write (9, 10) nerg
          write (6, 10) nerg
10        format (/' *** NERG =', i2,' > 25; ABORT ***')
          call exit
        end if
*  open units 10, 11, and 12 for storage of channel parameters, transformation
*  and quadrature matrices if more than one energy desired
cstart unix cray mac
        if (airyfl) then
          call tmpnm (10, xname)
* open scratch file (unit is therefore negativ here, see open)
* isize is only needed on a univac
          call openf(-10, xname, 'su', isize)
        endif
cend
        call tmpnm (11, xname)
        call openf(-11, xname, 'su', isize)
        call tmpnm (12, xname)
        call openf(-12, xname, 'sf', 0)
      end if
*   open files for storage of integral cross sections
      if (wrxsec .or. xsecwr) then
        do 60  ifile = 1, nerg
          nfile = 15 + ifile
          if (wrxsec) then
             call gennam (xname, jobnam, ifile, 'ics', lenx)
             if(ifile.eq.1) xnam1=xname
          else
            call tmpnm (nfile, xname)
          end if
          if (wrxsec) then
             call openf(nfile, xname, 'sf',0)
             rewind nfile
          else
             nnfile=-nfile
             call openf(nnfile, xname, 'sf',0)
          end if
60      continue
        if(wrxsec) then
          if (nerg .eq. 1) then
            write (9, 110) xnam1(1:lenx)
110         format (' ** INTEGRAL CROSS SECTIONS SAVED IN FILE ',(a))
          else
            write (9, 115) xnam1(1:lenx),xname(1:lenx)
115         format
     :    (' ** INTEGRAL CROSS SECTIONS SAVED IN FILES ',(a),
     :     ' THROUGH ',(a))
          end if
        end if
      end if
c.....open direct access file for interpolation and restart
*     only if partial or total cross sections desired, never if
*     photodissociation calculation or wavefunction calculation
      if (wrxsec .or. xsecwr .or. partw .or. wrpart) then
        if (.not. wavefl .and. .not. photof) then
           call dinit
           nfile = 14
           lenj=index(jobnam,' ')-1
           if (lenj .eq. 0) lenj=40
           if (lenj .gt. 8) then
             call gennam (xname, jobnam(1:8), 0, 'sav', lenx)
           else
             call gennam (xname, jobnam, 0, 'sav', lenx)
           endif
           if (rsflag) then
             inquire (file=xname, exist=existf)
             if (.not. existf) then
               write (9, 300) xname(1:lenx)
               write (6, 300) xname(1:lenx)
               call exit
             end if
           end if
           savfil=xname
           call dopen(1,nfile,savfil)
           write (9, 180) xname(1:lenx)
180        format (' ** RESTART INFORMATION SAVED IN FILE ',(a))
         endif
      endif
* open direct access file for storage of wavefunction
      if (wavefl.and. .not. boundc) then
        call dinit
        nfile = 15
        lenj=index(jobnam,' ')-1
        if (lenj .eq. 0) lenj=40
        if (lenj .gt. 8) then
          call gennam (xname, jobnam(1:8), 0, 'wfu', lenx)
        else
          call gennam (xname, jobnam, 0, 'wfu', lenx)
        endif
          if (rsflag) then
            inquire (file=xname, exist=existf)
            if (.not. existf) then
              write (9, 300) xname(1:lenx)
              write (6, 300) xname(1:lenx)
              call exit
            end if
          end if
        call dopen(2,nfile,xname)
        write (6, 210) xname(1:lenx)
        write (9, 210) xname(1:lenx)
210     format (' ** WAVEFUNCTION SAVED IN FILE ',(a))
      endif
*   open files for storage of partial cross sections
      if (wrpart) then
        do 230  ifile = 1, nerg
          nfile = 40 + ifile
          call gennam (xname, jobnam, ifile, 'pcs',lenx)
          if(ifile.eq.1) xnam1=xname
          inquire (file=xname, exist=existf)
          call openf(nfile, xname, 'sf', 0)
          rewind nfile
          if (rsflag.and.existf) then
            lenx=index(xname,' ')-1
            read (nfile, 220) oldlab
            read (nfile, 220) oldlab
220         format (1x, a)
            if (oldlab .ne. label) then
              write (9, 320) xname(1:lenx), xname, oldlab, label
              write (6, 320) xname(1:lenx), xname, oldlab, label
              call exit
            end if
            rewind nfile
            call fimove (nfile)
*            write (nfile, '('' ** RESTART **'')')
          end if
230     continue
      if (wrpart) then
        if (nerg .eq. 1) then
          write (9, 250) xnam1(1:lenx)
250       format (' ** PARTIAL CROSS SECTIONS SAVED IN FILE ',(a))
        else
          write (9, 260) xnam1(1:lenx),xname(1:lenx)
260       format
     :    (' ** PARTIAL CROSS SECTIONS SAVED IN FILES ',(a),
     :     ' THROUGH ',(a))
        end if
      end if
      end if
*   open files for accumulation of cs partial cross sections at each
*   projection index
      if (csflag) then
        do 280  ifile = 1, nerg
          nfile = 65 + ifile
          call tmpnm (nfile, xname)
          nnfile=-nfile
          call openf(nnfile, xname, 'su', 0)
280     continue
      end if
*  open files smatn for storage of selected s-matrix elements
      if (writs) then
        do 330  ifile = 1, nerg
        nfile = ifile + 90
          call gennam(xname, jobnam, ifile, 'smt', lenx)
          if(ifile.eq.1) xnam1=xname
          if (rsflag) then
            inquire (file=xname, exist=existf)
            if (.not. existf) then
              lenx=index(xname,' ')-1
              write (9, 300) xname(1:lenx)
              write (6, 300) xname(1:lenx)
300           format(/' RESTART ATTEMPTED, BUT FILE',a,
     :                ' DOES NOT EXIST')
              call exit
            end if
          end if
          call openf(nfile, xname, 'du', isizes)
          if (rsflag) then
* read s-matrix header
            newlab=label
            call rdhead(nfile,cdate,rsc1(1),rsc1(2),
     :                  lsc1(1),
     :                  lsc1(2),lsc1(3),lsc1(4),lsc1(5),isc1(1),isc1(2),
     :                  isc1(3),isc1(4),isc1(5),isc1(6),isc1(7),
     :                  isc1(8),isc1(9),isc2,isc3,rsc2,isc4)
            if (newlab .ne. label) then
              lenx=index(xname,' ')-1
              write (9, 320) xname(1:lenx), xname, label, newlab
              write (6, 320) xname(1:lenx), xname, label, newlab
320           format
     :     (/' *** LABEL IN FILE ',a, ' DOES NOT MATCH INPUT DATA',
     :         '; ABORT ***',/'     ', a5,':', a, /,
     :                        '     INPUT:', a)
              label=newlab
              call exit
            end if
          end if
330      continue
      if (writs) then
        if (nerg .eq. 1) then
          write (9, 340) xnam1(1:lenx)
340       format (' ** SELECTED S-MATRIX ELEMENTS SAVED IN FILE ',(a))
        else
          write (9, 350) xnam1(1:lenx),xname(1:lenx)
350       format
     :    (' ** SELECTED S-MATRIX ELEMENTS SAVED IN FILES ',(a),
     :     ' SMAT1 THROUGH ',(a))
        end if
      end if
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine restrt (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev,
     :               nchmax,rtmn1,rtmx1,dinsid, writs, csflag, nucros)
*  subroutine to read restart information
*  author:  h.-j. werner
*  latest revision:  17-jan-2008 by mha
* --------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      logical writs,csflag,nucros, llpar, wrpart, wrxsec
      character*40 oldlab
      integer jtot, nchmax
      character*40 input,output,jobnam,savfil
      include "common/parpot"
      common /cofile/ input,output,jobnam,savfil
      common /coener/ energ(1)
      common /cobmat/ q(1)
      common /cow/ q1(1)
      common /cozmat/ q2(1)
      common /coamat/ q3(1)
      common /colpar/ llpar(21),wrpart, wrxsec

      dimension word(30),iword(32),nlev(1)
      if (.not. wrpart .and. .not .wrxsec) then
        write (6, 5)
        write (9, 5)
5       format (' EITHER WRPART OR WRXSEC MUST HAVE BEEN TRUE IN',
     :    ' ORIGINAL DATA',/,'    FOR RESTART TO WORK; ABORT')
        call exit
      endif
      call drest(1)
      call dbri(iword,32,1,1)
      j=1
      do 10 i=1,10
      write(oldlab(j:j+3),'(A4)') iword(i)
10    j=j+4
      if (oldlab .ne. label) then
          write (9, 50) oldlab, label
          write (6, 50) oldlab, label
50        format
     : (/' *** LABEL IN FILE TRSTRT DOES NOT MATCH INPUT DATA',
     :     '; ABORT ***',/'     TRSTRT:', a, /,
     :                        '      INPUT:', a)
          call exit
      end if
      jtot=iword(11)
      if(iword(12).ne.jtotd) then
          write(6,60) iword(12),jtotd
60        format(/' RESTART: JTOTD=',i2,' RESET TO',i3)
      end if
      jtotd=iword(12)
      jlpar=iword(13)
      nu=iword(14)
      nutop=iword(15)
      if(iword(16).ne.nud) then
          write(6,70) iword(12),jtotd
70        format(/' RESTART: JTOTD=',i2,' RESET TO',i3)
      end if
      nud=iword(16)
      nerg=iword(17)
      nchmax=iword(18)
      csflag=iword(19).ne.0
      writs =iword(20).ne.0
      nucros=iword(21).ne.0
      jtop=iword(22)
      call dbrr(word,20,1,0)
      rtmn1=word(1)
      rtmx1=word(2)
      dinsid=word(3)
      write(6,75) jtot,jlpar,nu
      write(9,75) jtot,jlpar,nu
75    format(/' RESTART DATA FOUND FOR JTOT=',i3,'  JLPAR=',i2,
     1          '  NU=',i3)
      do 90 ien=1,nerg
      energ(ien)=word(10+ien)
      nlev(ien)=iword(22+ien)
      n=nlev(ien)
      irec=(ien-1)*5+2
      if(nucros) irec=irec+5*nerg
      write(6,80) energ(ien),nlev(ien)
      write(9,80) energ(ien),nlev(ien)
80    format(/' ENERGY=',f10.3,' NLEVOP=',i3)
      call intchk(irec,q,q1,q2,q3,jf,jp,n,n,nucros)
      if((.not.nucros.and.(jf.ne.jtot.or.jp.ne.jlpar))
     1     .or. (nucros.and.jf.ne.nu)) then
          write(6,85)
          write(9,85)
85        format(/' INCONSISTENT DATA DETECTED. ABORT')
          call exit
      end if
90    continue
      return
      entry rsave(jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev,
     :       nchmax,rtmn1,rtmx1,dinsid,writs, csflag, nucros)
      j=1
c.....reserve space for restart information
      call dres(64,1,1)
      do 95 i=1,10
      read(label(j:j+3),'(A4)') iword(i)
95    j=j+4
      iword(11)=jtot
      iword(12)=jtotd
      iword(13)=jlpar
      iword(14)=nu
      iword(15)=nutop
      iword(16)=nud
      iword(17)=nerg
      iword(18)=nchmax
      iword(19)=0
      iword(20)=0
      iword(21)=0
      iword(22)=jtop
      if(csflag) iword(19)=1
      if(writs)  iword(20)=1
      if(nucros) iword(21)=1
      word(1)=rtmn1
      word(2)=rtmx1
      word(3)=dinsid
      do 100 ien=1,nerg
      word(10+ien)=energ(ien)
      iword(22+ien)=nlev(ien)
100   word(8+ien)=energ(ien)
      call dbwi(iword,32,1,1)
      call dbwr(word,20,1,0)
      call dbwc(1,1)
      call dclos(1)
      call dopen(1,14,savfil)
      return
      end
      subroutine inbfi(ia,n)
*  buffer routine for s-matrices
*  author: h.j. werner
*  current revision date: 14-jun-1990 by mha
      implicit double precision (a-h,o-z)
      common/clseg/ lseg,intrel,lchar
      common/cobuf/ lbuf,ibuf(7)
      common/cofil/ nfl,iofbuf,maxrec(60),iofrec(60),nwrec
      dimension ia(1)
*
*  read n elements of array ia from buffer
*
*  lseg = sector length (as in open statement), provided in common clseg
*  intrel = number of integer words per real word
*  lchar = number of characters per integer word
*  lbuf = max. buffer length (should be multiple of lseg)
*  nfl = current file
*  iofbuf = current offset in buffer
*  iofrec(ifile) = offset in file
*
*  maxrec(ifile) = current number of s-matrices on file
*  nwrec = number of words for one logical record (each jtot one record)
*            (should be multiple of lseg)
*
      nr=n
      iofa=0
      goto 10
*
      entry inbfr(ia,n)
      nr=n*intrel
      iofa=0
10    if(nr.le.0) return
      if(iofbuf.eq.lbuf) then
          call rdabsf(nfl,ibuf,lbuf/intrel,iofrec(nfl)/intrel)
          call fwait(nfl)
          iofrec(nfl)=iofrec(nfl)+lbuf
          iofbuf=0
      end if
      m=min0(nr,lbuf-iofbuf)
      do 20 i=1,m
20    ia(iofa+i)=ibuf(iofbuf+i)
      iofa=iofa+m
      iofbuf=iofbuf+m
      nr=nr-m
      goto 10
*
* fill first buffer for given jtot
*
      entry readrc(iadr,nfile,lrec,jtot,jlpar,nu,nopen,length,nnout)
      nfl=nfile
*
*  if iadr.eq.0 read next record
*  if iadr.gt.0 read at address iadr
*  it is assumed here that lseg.ge.7
*
      if(iadr.gt.0) iofrec(nfl)=iadr
      call rdabsf(nfl,ibuf,lseg/intrel,iofrec(nfl)/intrel)
      call fwait(nfl)
* ibuf=-1 signals end of file
      if (ibuf(1).eq. -1) then
        lrec=-1
        return
      endif
      call rdabsf(nfl,ibuf,lbuf/intrel,iofrec(nfl)/intrel)
      call fwait(nfl)
      iadr=iofrec(nfl)
      iofrec(nfl)=iofrec(nfl)+lbuf
      lenbuf=lseg
      iofbuf=7
      nwrec=ibuf(1)-lseg
      lrec=ibuf(1)
      jtot=ibuf(2)
      jlpar=ibuf(3)
      nu=ibuf(4)
      nopen=ibuf(5)
      length=ibuf(6)
      nnout=ibuf(7)
      return
      end
*---------------------------------------------------------------------
      subroutine movrc8(ia,text)
      character*(*) text
      dimension ia(1)
      l=len(text)
      do 10 j=1,l
10    text(j:j)=char(ia(j))
      return
      end
      subroutine movcr8(text,ia)
      character*(*) text
      dimension ia(1)
      l=len(text)
      do 10 j=1,l
10    ia(j)=ichar(text(j:j))
      return
      end
* --------------------------------------------------------------------
      subroutine rdhead(nfile,cdate,ered,rmu,csflag,flaghf,
     :  flagsu,twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud,
     :  nlevel,nlevop,nnout,jlev,inlev,elev,jout)
*
*  subroutine to read header from file nfile
*  authors: h.j. werner and b. follmeg
*  current revision date: 27-oct-95
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag, flaghf, flagsu, twomol, nucros
      character*20 cdate
      include "common/parpot"
      common/clseg/ lseg,intrel,lchar
      common/cofil/ nfl,iofbuf,maxrec(60),iofrec(60),nwrec,lenbuf
      common/cobuf/ lbuf,ibuf(1)
      dimension idate(20),ilabel(48), ipotnm(48)
      dimension jlev(1), inlev(1), elev(1), jout(1)
*
      nfl=nfile
      print *,' nfl in rdhead:  ', nfl
      call rdabsf(nfl,ibuf,lbuf/intrel,0)
      print *, 'after rdabsf'
      call fwait(nfl)
      print *, 'after fwait'
      lenbuf=lbuf
      iofbuf=0
      nwrec=ibuf(1)-lseg
      iofrec(nfl)=lbuf
*
*     first word on file contains the length of the header
*
      ldate = 20
      llabel= 48
      lpotnm=48
      call inbfi(lenhd,1)
      call inbfi(idate,ldate)
      call inbfi(ilabel,llabel)
      call inbfi(ipotnm,lpotnm)
      call movrc8(idate,cdate)
      call movrc8(ilabel,label)
      call movrc8(ipotnm,potnam)
      call inbfr(ered,1)
      call inbfr(rmu,1)
      call inbfi(csflag,1)
      call inbfi(flaghf,1)
      call inbfi(flagsu,1)
      call inbfi(twomol,1)
      call inbfi(nucros,1)
      print *, 'nucros:  ', nucros
      call inbfi(jfirst,1)
      call inbfi(jfinal,1)
      call inbfi(jtotd,1)
      call inbfi(numin,1)
      call inbfi(numax,1)
      call inbfi(nud,1)
      call inbfi(nlevel,1)
      call inbfi(nlevop,1)
      print *, 'nlevop:  ', nlevop
      call inbfi(nnout,1)
      print *, 'nnout:  ', nnout
      call inbfi(jlev,nlevel)
      call inbfi(inlev,nlevel)
      print *, 'bfore inbfi inlev'
      call inbfr(elev,nlevel)
      call inbfi(jout,iabs(nnout))

      return
      end
* --------------------------------------------------------------------
      subroutine saddr(nfile,jfirst,jfinal,numin,numax,csflag,
     :                 iadr,nmax,jfsts,jlparf,jlpars,ierr)
* subroutine to built up a pointer table needed for direct access io
* author: h.j. werner and b. follmeg
* current revision date 6-sept-88
* variables in call list:
*  on input:
*  nfile:   logical unit number
*  jfirst:  the first jtot value for the first parity
*  jfinal:  the maximum value of jtot
*  numin:   min. value of coupled states proj. index
*  numax:   max. value of coupled states proj. index
*  csflag:  coupled states flag
*  nmax:    dimension of array iadr
*
*  on return:
*  iadr:    array, contains addresses of jtot records
*  jfsts:   the first jtot value for the second parity (if there is any,
*           else jfsts is set to -1)
*  jlparf:  the first parity
*  jlpars:  second parity, if not present jlpars = 0
*  ierr:    error flag, if .ne. 0 an error has occured
*
      logical csflag
      integer nfile, iadr, nmax, jfirst, jfsts, jfinal, numin,
     :        numax, ierr, jlpars, jlparf
      integer jln, jlp, nwaves, nrec, jlpold, nparit, jj, iaddr,
     :        lrec, jtot, jlpar, nu, nopen, length
      dimension iadr(1)
*
      ierr = 0
      jln=2
      if(csflag) then
        jln=1
      else
        numin=0
        numax=0
      end if
* nwaves is number of partial waves
      nwaves = jfinal - jfirst + 1
      nrec = nwaves * (numax-numin+1) * jln
      if(nrec.gt.nmax) then
        write(6,5) nrec,nmax
5       format(/' *** ADDRESS ARRAY TOO SMALL IN SADDR:',2i6,' ABORT')
        ierr=1
        return
      end if
* initialize address array
      do 10 i=1,nrec
10    iadr(i)=-1
*
      jlp = 1
      jlpold = -10
      nparit = 0
      jfsts=-1
      jlpars=0
      iaddr=0
*
20    call readrc(iaddr,nfile,lrec,jtot,jlpar,nu,nopen,length,nnout)
* end of file detected ?
      if(lrec.lt.0) goto 40
*
      if(jlpar.ne.jlpold) then
        nparit = nparit + 1
        jlpold = jlpar
        if(jlpar.eq.1) then
          jlp = 0
        else
          jlp = 1
        end if
        if(nparit .eq. 1) then
           jlparf = jlpar
        else if (nparit .eq. 2) then
           jfsts = jtot
           jlpars = jlpar
        else
           write(6,25)
25         format(/' ** PARITY ERROR IN SADDR, ABORT')
           ierr = 1
           return
        end if
      end if
      jj = (jlp + nu - numin) * nwaves + jtot + 1
      iadr(jj)=iaddr
*      write(6,30) jtot,jlpar,nu,jj,iaddr,jfsts
*30    format(' jtot=',i3,'  jlpar=',i2,'  nu=',i3,'  jj=',i3,
*     :       ' iad=',i6, ' jfsts=',i6)
      iaddr=iaddr+lrec
      goto 20
40    return
      end
* -----------------------------------------------------------------------
      subroutine sread (iadr, sreal, simag, jtot, jlpar, nu,
     :                  jq, lq, inq, inpack, jpack, lpack,
     :                  nfile, nmax, nopen, length, ierr)
* authors: h.j. werner and b. follmeg
* current revision date: 21-feb-2006 by mha
*     read real and imaginary parts of s-matrix together with
*     other information as written in soutpt, swrite
*     if iadr = 0 read sequential (next record)
*     if iadr > 0 read absolute
*     if nopen = -1, the lower triangle is filled
      implicit double precision (a-h,o-z)
      logical triang
      dimension sreal(nmax,1), simag(nmax,1),
     1    jpack(1), lpack(1),inpack(1),jq(1),lq(1),inq(1)
*  variable in common block /coselb/
*     ibasty    basistype
      common /coselb/ ibasty
      common  /coisc7/ isc7(9)
      common  /coisc8/ j12q(9)
      ierr=0
      triang =.false.
      if(nopen.lt.0) then
        triang = .true.
        nopen = iabs(nopen)
      end if
* on return:
* the vectors jpack, lpack, inpack will hold the column indices of the
*             packed basis (dimension length)
* the vectors jq, lq, inq will hold the row indices of the packed
*             basis (dimension nopen)
* read next s-matrix header
      iaddr = iadr
      print *, 'before call readrc, nfile = ', nfile
      call readrc(iaddr,nfile,lrec,jtot,jlpar,nu,nopen,length,nnout)
      if(lrec.lt.0) goto 400
*
      call inbfi(jpack, length)
      call inbfi(lpack, length)
      call inbfi(inpack,length)
      if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :   call inbfi(isc7,length)
      if (nnout .gt. 0) then
         do 50 i = 1, length
         jq(i) = jpack(i)
         lq(i) = lpack(i)
         inq(i)= inpack(i)
      if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :   j12q(i)=isc7(i)
50       continue
         nopen = length
         if (triang) then
            ioff =1
            do 70 irow = 1, length
               call inbfr(sreal(ioff,1), irow)
               call inbfr(simag(ioff,1), irow)
               ioff = ioff + irow
70          continue
            return
         end if
* read s-matrix
         do 80  icol = 1, length
            call inbfr(sreal(1,icol), icol)
            call inbfr(simag(1,icol), icol)
80       continue
* fill lower triangle
         do 90 icol=1,length
         do 90 irow=1,icol
         sreal(icol,irow)=sreal(irow,icol)
90       simag(icol,irow)=simag(irow,icol)
      else if (nnout .le. 0) then
*  here if you have written out columns of the s-matrix
         call inbfi(jq, nopen)
         call inbfi(lq, nopen)
         call inbfi(inq,nopen)
         if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :      call inbfi(j12q,nopen)
*  now read columns of the s-matrix
         do 140  icol = 1, length
            call inbfr(sreal(1,icol), nopen)
            call inbfr(simag(1,icol), nopen)
140      continue
      end if
      return
400   ierr=lrec
      return
      end
*--------------------------------------------------------------------
      subroutine swrite (sreal, simag, jtot, jlpar, nu,
     :                   jq, lq, inq, iorder, inpack, jpack, lpack,
     :                   epack, nfile, nmax, nopen)
*  subroutine to write selected elements of s-matrix to file nfile
*  author:  millard alexander
*  modified by  h.j. werner and b. follmeg
*  current revision date: 21-feb-2006 by mha
*  ------------------------------------------------------------------
*  variables in call list:
*    sreal:     on entry: contains real part of open-channel s-matrix
*               on return: contains real part of packed s-matrix
*    simag:     on entry: contains imaginary part of open-channel s-matrix
*               on return: contains imaginary part of packed s-matrix
*    jtot:      total angular momentum
*    csflag:    if .true. coupled-states calculation
*               if .false. close-coupled calculation
*    flaghf:    if .true., then system with half-integer spin
*                if .false., then system with integer spin
*    nu:        cs projection index (not used in cc calculation)
*    jq:        channel rotational angular momenta
*    lq:        channel orbital angular momenta
*    inq:       additional quantum index of each channel
*    note!!!   if flaghf = .true., then the true values
*    of the rotational quantum numbers, the total angular momentum,
*    and the coupled-states projection index are equal to the values
*    stored in jq, jtot, and nu plus 1/2
*    inpack,jpack,
*    lpack, epack
*    nfile:     logical unit for output of s-matrices
*    nmax:      maximum row dimension of matrices
*    nopen:     number of channels
*  variables in common block /cosout/
*    nnout:     number of different rotational levels for which s-matrix
*               elements are to be saved in file nfile
*    jout(i):   values of rotational angular momentum for these levels
*               if nnout is positive, than an s-matrix element will be saved
*               only if both the initial and final quantum numbers correspond
*               to one of the values of jout(i)
*               if nnout is negative, then every column of the s-matrix for
*               which the initial quantum numbers correspond to one of the
*               values of jout(i) will be printed
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units (mass of electron = 1)
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  ------------------------------------------------------------------
      use mod_cosout, only: nnout, jout
      use mod_coeint, only: eint
      implicit double precision (a-h,o-z)
      integer ic, icol, ii, ir, irow, jtot, jlpar, length, nmax,
     :        nopen, nfile, nu, mmout
      integer jq, jpack, lq, lpack, inq, inpack
*      real simag, sreal, eint, epack, ered, rmu
      common /clseg/ lseg,intrel,lchar
      common /coered/ ered, rmu
*  variable in common block /coselb/
*     ibasty    basistype
      common /coselb/ ibasty
      common  /coisc7/ isc7(9)
      common  /coisc8/ j12q(9)
      dimension sreal(nmax,nmax), simag(nmax,nmax),
     :          jq(1), lq(1), inq(1), jpack(1), lpack(1),
     :          epack(1), inpack(1), iorder(1)
*  the vector iorder will point to the position in the unpacked basis of each
*  state in the packed basis
*  the vector jpack will hold the rotational quantum numbers in the packed bas
*  the vector lpack will hold the orbital angular momenta of each channel in
*  the packed basis
*  the vector epack will hold the channel energies in the packed basis
*  the vector inpack will hold the symmetry indices in the packed basis
*  first sum over the unpacked states
        mmout=iabs(nnout)
        length = 0
        do 30  icol = 1, nopen
*  now sum over the packed states, find labels
          do 20  ii = 1, mmout
            if (jq(icol) .eq. jout(ii) ) then
*  here if match
              length = length + 1
              jpack(length) = jq(icol)
              lpack(length) = lq(icol)
              epack(length) = eint(icol)
              inpack(length) = inq(icol)
              if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :             isc7(length) = j12q(icol)
              iorder(length) = icol
              go to 30
            end if
20        continue
30      continue
*  calculate number of words that will be written
        nrecw = 7 + 3 * length
        if(nnout.gt.0) then
           nrecw = nrecw+length*(length+1)*intrel
        else
           nrecw = nrecw+3*nopen+2*length*nopen*intrel
        end if
*  write out general information on next record
        call nexrec(nfile,nrecw,jtot,jlpar,nu,nopen,length,nnout)
        call outbfi(jpack,length)
        call outbfi(lpack,length)
        call outbfi(inpack,length)
        if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :             call outbfi(isc7,length)
*  here we pack the s-matrix and print out just those elements for which
*  the initial and final rotational quantum numbers correspond to an element
*  in the array jout
      if (nnout .gt. 0) then
*  the dimension of the packed s-matrix is length x length
*  now pack the real part of the s-matrix
        do 45  icol = 1, length
          ic = iorder(icol)
          do 40  irow = 1, length
            ir = iorder(irow)
            sreal(irow,icol) = sreal(ir,ic)
            simag(irow,icol) = simag(ir,ic)
40        continue
45      continue
* write s-matrix into buffer
        do 80  icol = 1, length
           call outbfr(sreal(1,icol),icol)
           call outbfr(simag(1,icol),icol)
80      continue
*  here if you want to print out columns of the s-matrix
      else if (nnout .le. 0) then
         call outbfi(jq,nopen)
         call outbfi(lq,nopen)
         call outbfi(inq,nopen)
         if ((ibasty.eq.12).or.(ibasty.eq.13).or.(ibasty.eq.15))
     :             call outbfi(j12q,nopen)
*  now write out columns of the s-matrix into buffer
*  length is the number of columns of the s-matrix to save
        do 140  ii = 1, length
          icol = iorder(ii)
          call outbfr(sreal(1,icol),nopen)
          call outbfr(simag(1,icol),nopen)
140     continue
      end if
      call clear
      return
      end
*--------------------------------------------------------------------
      subroutine outbfi(ia,n)
*  buffer routine for s-matrices
*  author: h.j. werner
*  current revision date: 23-jan-90
      implicit double precision (a-h,o-z)
      dimension ia(1)
      logical ldum,rsflag,lpar
      common /colpar/ ldum(15),rsflag,lpar(11)
      common/clseg/ lseg,intrel,lchar
      common/cobuf/ lbuf,ibuf(7)
      common/cofil/ nfl,iofbuf,maxrec(60),iofrec(60),nwrec,lenbuf
*
*  buffer out n elements of array ia
*
*  lseg = sector length (as in open statement), provided in common clseg
*  intrel = number of integer words per real word
*  lchar = number of characters per integer word
*  lbuf = buffer length (should be multiple of lseg)
*  nfl = current file
*  iofbuf = offset in buffer
*  iofrec(ifile) = offset in file
*  maxrec(ifile) = current number of s-matrices on file
*  nwrec = number of words per logical record (should be multiple of lseg)
*  (n.b. the arrays for one jtot-value will be saved on one logical record,
*  where the number of logical records is not the same as the number of
*  physical records)
      nr=n
      iofa=0
      goto 10
*
      entry outbfr(ia,n)
      nr=n*intrel
      iofa=0
*
10    m=min0(nr,lbuf-iofbuf)
      do 20 i=1,m
20    ibuf(iofbuf+i)=ia(iofa+i)
      iofa=iofa+m
      iofbuf=iofbuf+m
      nr=nr-m
* write out buffer if necessary
      if(iofbuf.eq.lbuf) then
        call wrabsf(nfl,ibuf,lbuf/intrel,iofrec(nfl)/intrel)
        call fwait(nfl)
        iofrec(nfl)=iofrec(nfl)+lbuf
        iofbuf=0
      end if
*
      if(nr.gt.0) goto 10
      return
*
*  clear last buffer
*
      entry clear
      print *, 'nfile in clear:  ', nfl
*
      if(iofbuf.gt.0) then
        iadr=iofrec(nfl)
c       l=((iofbuf-1)/lseg+1)*lseg
c       call wrabsf(nfl,ibuf,l/intrel,iofrec(nfl)/intrel)
        call wrabsf(nfl,ibuf,lbuf/intrel,iofrec(nfl)/intrel)
        call fwait(nfl)
c       iofrec(nfl)=iofrec(nfl)+l
        iofrec(nfl)=iofrec(nfl)+lbuf
      end if
*  avoid writing end of file if header is replaced in restart
      if(rsflag.and.iadr.eq.0) return
*
*  write dummy record indicating eof
*
      ibuf(1)=-1
      call wrabsf(nfl,ibuf,lseg/intrel,iofrec(nfl)/intrel)
      call fwait(nfl)
      return
*
*  start new s-matrix record. it is assumed that buffer is cleared!
*
      entry nexrec(nfile,nrecw,jtot,jlpar,nu,nopen,length,nnout)
*
      nwrec=((nrecw-1)/lseg+1)*lseg
      lenbuf=0
      ibuf(1)=nwrec
      ibuf(2)=jtot
      ibuf(3)=jlpar
      ibuf(4)=nu
      ibuf(5)=nopen
      ibuf(6)=length
      ibuf(7)=nnout
      iofbuf=7
      nfl=nfile
      maxrec(nfl)=maxrec(nfl)+1
      print *, ' in nexrec, nfl, maxrec:  ', nfl, maxrec(nfl)
*  start address for each s-matrix is multiple of lseg
      return
      end
*---------------------------------------------------------------------
      subroutine wrhead(nfile,cdate,
     :  ered,rmu,csflag,flaghf,
     :  flagsu,twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud,
     :  nlevel,nlevop,nnout,jlev,inlev,elev,jout)
*
*  initialize buffering for file nfile and write general information
*  (header)
*
* author: h.j. werner
*  current revision date: 27-oct-1995 by mha
*---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag, flaghf, flagsu,twomol, nucros
      character*20 cdate
      include "common/parpot"
      common/clseg/ lseg,intrel,lchar
      common/cofil/ nfl,iofbuf,maxrec(60),iofrec(60),nwrec,lenbuf
*
*  equivalence to avoid trouble with string descriptors in buffering
*
      dimension jlev(1),inlev(1),elev(1)
      dimension idate(20),ilabel(48),ipotnm(48)
*
      call movcr8(label,ilabel)
      call movcr8(potnam,ipotnm)
      call movcr8(cdate,idate)
      iofrec(nfile)=0
      maxrec(nfile)=0
      iofbuf=0
      lenbuf=0
      nfl=nfile
*
      ldate = 20
      llabel= 48
      lpotnm = 48
      lenhd=ldate+llabel+lpotnm+
     :      2*intrel+15+(2+intrel)*nlevel+iabs(nnout)
      call outbfi(lenhd,1)
      call outbfi(idate,ldate)
      call outbfi(ilabel,llabel)
      call outbfi(ipotnm,lpotnm)
      call outbfr(ered,1)
      call outbfr(rmu,1)
      call outbfi(csflag,1)
      call outbfi(flaghf,1)
      call outbfi(flagsu,1)
      call outbfi(twomol,1)
      call outbfi(nucros,1)
      call outbfi(jfirst,1)
      call outbfi(jfinal,1)
      call outbfi(jtotd,1)
      call outbfi(numin,1)
      call outbfi(numax,1)
      call outbfi(nud,1)
      call outbfi(nlevel,1)
      call outbfi(nlevop,1)
      call outbfi(nnout,1)
      call outbfi(jlev,nlevel)
      call outbfi(inlev,nlevel)
      call outbfr(elev,nlevel)
      call outbfi(jout,iabs(nnout))
      call clear
      return
      end
* -------------------------------------------------------------------
      subroutine dread(ii,l,ifil,irec,iof)
* --------------------------------------------------------------
* subroutine to control direct access i/o of partial and integral cross
* sections, restart information, and wavefunction
* author:  h-j werner
* latest revision:  6-jun-1991 by hjw (maxrec extended 7-feb-1992)
* --------------------------------------------------------------
c
c.....read "l" integer words from record "irec" on file "ifil"
c.....with offset iof
c
*  variables in common block /cobuf/
*    lbuf:    length of i/o buffer
*    ibuf:      i/o buffer
*  variable in common block /cosize/
*    isize:     size of files (only needed for univac, optional on vax)
*    isizes:    size of s-matrix file (only needed for univac, optional
*               on vax
* --------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) name
      parameter (maxun=2, maxrec=5000)
      include "common/parllbf"
      common/cdbf/ ldbuf,libuf,ibfil,ibrec,ibof,ibstat,idbuf(llbuf)
      common/clseg/ lseg,intrel,lchar
      common /cobuf/ lbuf,ibuf(1)
      common /cdio/iadr(maxrec,maxun),len(maxrec,maxun),next(maxun),
     1             iun(maxun),iostat(maxun),last(maxun),lhea,junk
      common /cosize/ isize,isizes
      dimension ii(1)
      if(iun(ifil).eq.0) then
        write(6,10) ifil,irec
10      format(/' DIRECT READ ERROR ON FILE',i2,' RECORD',i3,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      if(irec.gt.maxrec) then
        write(6,15) irec,maxrec
15      format(/' DIRECT READ ERROR ON FILE',i2,' RECORDNUMBER',i3,
     1          ' OUT OF RANGE. ALLOWED',i4)
        call exit
      end if
      if(iadr(irec,ifil).lt.0) then
        write(6,20) ifil,irec
20      format(/' DIRECT READ ERROR ON FILE',i2,' RECORD',i3,
     1          ' UNDEFINED')
        call exit
      end if
      if(iof+l.gt.len(irec,ifil)) then
        write(6,30) ifil,irec,iof+l,len(irec,ifil)
30      format(/' DIRECT READ ERROR ON FILE',i2,' RECORD',i3,
     1          ' TRY TO READ',i6,' WORDS, RECORDLENGTH=',i6,' WORDS')
        call exit
      end if
      ir=iadr(irec,ifil)+iof
      if(mod(l,intrel).ne.0.or.mod(ir,intrel).ne.0) then
        write(6,31) ifil,irec,l,ir
31      format(/' DIRECT READ ERROR ON FILE',i2,' RECORD',i3,
     1    ' LENGTH=',i4,' OR ADRESS=',i6,' NOT MULTIPLE OF INTREL')
        call exit
      end if
      ll=l/intrel
      ir=ir/intrel
      call rdabsf(iun(ifil),ii,ll,ir)
      call fwait(iun(ifil))
      iostat(ifil)=0
      return
*
      entry dres(l,ifil,irec)
c
c.....reserve "l" real*8  words for record "irec" on file "ifil"
c
      ll=l*intrel
      if(iun(ifil).eq.0) then
        write(6,40) ifil,irec
40      format(/' DIRECT RESERVE ERROR ON FILE',i2,' RECORD',i3,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      if(irec.gt.maxrec) then
        write(6,45) irec,maxrec
45      format(/' DIRECT RESERVE ERROR ON FILE',i2,' RECORDNUMBER',i3,
     1          ' OUT OF RANGE. ALLOWED',i4)
        call exit
      end if
      if(iadr(irec,ifil).lt.0) then
        iadr(irec,ifil)=next(ifil)
        next(ifil)=next(ifil)+((ll-1)/lbuf+1)*lbuf
        len(irec,ifil)=next(ifil)-iadr(irec,ifil)
      end if
      if(ll.gt.len(irec,ifil)) then
        write(6,50) ifil,irec,ll,len(irec,ifil)
50      format(/' DIRECT RESERVE ERROR ON FILE',i2,' RECORD',i3,' TRY',
     1          ' TO RESERVE',i6,' WORDS, RECORDLENGTH=',i6,' WORDS')
        call exit
      end if
      return
*
      entry dwrite(ii,l,ifil,irec,iof)
c
c.....write "l" integer words to record "irec" on file "ifil"
c.....with offset iof
c
      if(iun(ifil).eq.0) then
        write(6,55) ifil,irec
55      format(/' DIRECT WRITE ERROR ON FILE',i2,' RECORD',i3,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      if(irec.gt.maxrec) then
        write(6,60) irec,maxrec
60      format(/' DIRECT WRITE ERROR ON FILE',i2,' RECORDNUMBER',i3,
     1          ' OUT OF RANGE. ALLOWED',i4)
        call exit
      end if
      if(iadr(irec,ifil).lt.0) then
        iadr(irec,ifil)=next(ifil)
        next(ifil)=next(ifil)+((l-1)/lbuf+1)*lbuf
        len(irec,ifil)=next(ifil)-iadr(irec,ifil)
      end if
      if(iof+l.gt.len(irec,ifil)) then
        if(irec.ge.last(ifil)) then
          next(ifil)=iadr(irec,ifil)+((iof+l-1)/lbuf+1)*lbuf
          len(irec,ifil)=next(ifil)-iadr(irec,ifil)
        else
      write (6,*) 'ii,l,ifil,irec,iof', ii,l,ifil,irec,iof
          write(6,65) ifil,irec,iof+l,len(irec,ifil)
65        format(/' DIRECT WRITE ERROR ON FILE',i2,' RECORD',i3,
     1          ' TRY TO WRITE',i6,' WORDS, RECORDLENGTH=',i6,' WORDS')
          call exit
        end if
      end if
      ir=iadr(irec,ifil)+iof
      if(mod(l,intrel).ne.0.or.mod(ir,intrel).ne.0) then
        write(6,66) ifil,irec,l,ir
66      format(/' DIRECT WRITE ERROR ON FILE',i2,' RECORD',i3,
     1    ' LENGTH=',i4,' OR ADRESS=',i6,' NOT MULTIPLE OF INTREL')
        call exit
      end if
      ll=l/intrel
      ir=ir/intrel
      call wrabsf(iun(ifil),ii,ll,ir)
      call fwait(iun(ifil))
      last(ifil)=max(last(ifil),irec)
      iostat(ifil)=0
      return
*
      entry dinit
c
c.....initialize direct i/o. dinit must be called before any call to
c.....another of these i/o routines
c
      if (llbuf .ne. lbuf) then
        write (6, 67) lbuf, llbuf
67      format ('** I/O BUFFER LENGTH OF', i4, ' IN BLOCK DATA IO .NE.',
     :        /,'   I/O BUFFER LENGTH OF', i4, ' IN DREAD; ABORT')
        call exit
      endif
      lhea=maxun*(2*maxrec+4)+2
      lnex=((lhea-1)/lbuf+1)*lbuf
      ibfil=0
      ibrec=0
      ibof=0
      ibstat=0
      libuf=lbuf
      ldbuf=lbuf/intrel
      do 70 if=1,maxun
      next(if)=lnex
      iun(if)=0
      last(if)=0
      iostat(if)=0
      do 70 ir=1,maxrec
      len(ir,if)=0
70    iadr(ir,if)=-1
      return
*
      entry dopen(ifil,iunit,name)
c
c.....open file "ifil" as unit "iunit" with filename "name"
c
      if(iun(ifil).eq.iunit) return
      if(iun(ifil).ne.0) then
        write(6,75) ifil,iunit,iun(ifil)
75      format(/' DIRECT OPEN ERROR FOR FILE',i2,' UNIT',i3,
     1          ' ALREADY OPEN AS UNIT',i3)
        call exit
      end if
      if(iunit.le.0) then
        write(6,80) ifil,iunit
80      format(/' DIRECT OPEN ERROR FOR FILE',i2,
     1          ' ILLEGAL UNIT NUMBER:',i4)
        call exit
      end if
      isize=0
      call openf(iunit,name,'DU',isize)
      iun(ifil)=iunit
      iostat(ifil)=0
      return
*
      entry dsave(ifil)
c
c.....save file information for restart
c
      if(iun(ifil).eq.0) then
        write(6,85) ifil
85      format(/' DIRECT SAVE ERROR ON FILE',i2,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      iofh=0
      do 86 ib=1,lhea,lbuf
      lh=min(lhea-iofh,lbuf)
      call imove(iadr(iofh+1,ifil),idbuf,lh)
      call wrabsf(iun(ifil),idbuf,ldbuf,iofh/intrel)
      call fwait(iun(ifil))
86    iofh=iofh+lbuf
      iostat(ifil)=0
      return
*
      entry drest(ifil)
c
c.....restart file "ifil"
c
      if(iun(ifil).eq.0) then
        write(6,90) ifil
90      format(/' DIRECT RESTART ERROR ON FILE',i2,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      iofh=0
      do 91 ib=1,lhea,lbuf
      ll=min(lhea-iofh,lbuf)
      call rdabsf(iun(ifil),idbuf,ldbuf,iofh/intrel)
      call fwait(iun(ifil))
      call imove(idbuf,iadr(iofh+1,ifil),ll)
91    iofh=iofh+lbuf
      iostat(ifil)=0
      return
*
      entry dclos(ifil)
c
c.....save file information and close file "ifil"
c
      if(iun(ifil).eq.0) then
        write(6,95) ifil
95      format(/' DIRECT CLOSE ERROR ON FILE',i2,
     1          ' UNIT NOT OPEN')
        call exit
      end if
      iofh=0
      do 96 ib=1,lhea,lbuf
      ll=min(lhea-iofh,lbuf)
      call imove(iadr(iofh+1,ifil),idbuf,ll)
      call wrabsf(iun(ifil),idbuf,ldbuf,iofh/intrel)
      call fwait(iun(ifil))
96    iofh=iofh+lbuf
      call closf(iun(ifil))
      iostat(ifil)=0
      iun(ifil)=0
      return
      end
* -------------------------------------------------------------------
      subroutine imove(ia,ib,n)
      dimension ia(1),ib(1)
      do 10 i=1,n
10    ib(i)=ia(i)
      return
      end
* -------------------------------------------------------------------
      subroutine dbri(ii,l,ifil,irec)
c
c.....sequentially buffer in "l" integer words from record "irec" on file "ifi
c.....positive irec indicates start at beginning
c.....zero or negative record means continue at present position
c
      implicit double precision (a-h,o-z)
      common/clseg/ lseg,intrel,lchar
      common/cdbf/ ldbuf,libuf,ibfil,ibrec,ibof,ibstat,idbuf(1)
      dimension ii(l)
      lre=l
      goto 5
c
      entry dbrr(ii,l,ifil,irec)
c
c.....sequentially buffer in "l" real words from record "irec" on file "ifil"
c.....positive irec indicates start at beginning
c.....zero or negative record means continue at present position
c
      lre=l*intrel
c
5     if(irec.gt.0) then
        if(ibstat.eq.1) call dwrite(idbuf,libuf,ibfil,ibrec,ibadr)
        ibstat=0
        ibof=0
        ibadr=0
        ibrec=irec
        ibfil=ifil
        call dread(idbuf,libuf,ifil,ibrec,ibadr)
      end if
      iof=0
10    len=min(lre,libuf-ibof)
      do 30 i=1,len
30    ii(iof+i)=idbuf(ibof+i)
      iof=iof+len
      lre=lre-len
      ibof=ibof+len
      if(lre.gt.0) then
        ibof=0
        ibadr=ibadr+libuf
        call dread(idbuf,libuf,ifil,ibrec,ibadr)
        goto 10
      end if
      return
c
      entry dbwi(ii,l,ifil,irec)
c
c.....sequentially buffer out "l" integer words to record "irec" on file "ifil
c.....positive irec indicates start at beginning or record
c.....zero or negative record means start at present position
c
      lre=l
      goto 50
c
      entry dbwr(ii,l,ifil,irec)
c
c.....sequentially buffer out "l" real words to record "irec" on file "ifil"
c.....positive irec indicates start at beginning of record
c.....zero or negative record means start at present position
c
      lre=l*intrel
c
50    if(irec.gt.0) then
        if(ibstat.gt.0) call dwrite(idbuf,libuf,ibfil,ibrec,ibadr)
        ibof=0
        ibadr=0
        ibrec=irec
        ibfil=ifil
        ibstat=1
      end if
      iof=0
70    len=min(lre,libuf-ibof)
      do 80 i=1,len
80    idbuf(ibof+i)=ii(iof+i)
      iof=iof+len
      lre=lre-len
      ibof=ibof+len
      if(lre.gt.0) then
        ibof=0
        call dwrite(idbuf,libuf,ifil,ibrec,ibadr)
        ibadr=ibadr+libuf
        goto 70
      end if
      return
c
      entry dbwc(ifil,irec)
c
c.....write out pending buffer
c
      call dwrite(idbuf,libuf,ibfil,ibrec,ibadr)
      ibstat=0
      return
      end
*********************************************************************
*                                                                   *
*                        i/o routines library                       *
*  these routines employ standart fortran i/o and may be quite      *
*  inefficient. If possible, use system i/o for direct access       *
*  lseg must be set in main program (file driver.f) appropriately   *
*  lseg should typically be 512. Note that the record length in the *
*  open statements are measured in bytes.                           *
*                                                                   *
*********************************************************************
*                         routines included:                        *
*                                                                   *
*   1. assgn    assigns and opens direct access files               *
*   2. rdabsf (wrabsf/fwait) read/write absolute routines           *
*   3. tmpnm    generates a unique filename                         *
*********************************************************************
cstart none
c;      subroutine assgn(luni,filnam,isize,icat)
c;*
c;*  routine to assign and open files
c;*
c;*  luni: logical unit number
c;*  luni .le.10:  formatted sequential files
c;*  luni .gt.10:  unformatted direct access files
c;*  filnam: filename
c;*  isize: size in tracks (output, not working with fortran i/o)
c;*  icat:  if icat=0, scratch
c;*         if icat>0, permanent
c;      character*(*) filnam
c;      character*14 nam
c;      character*12 stat
c;      character*(*) name
c;      common/disc/ ipos(98),iun(98),iostat(98),icatf(98),nam(98)
c;      logical openfl,exstfl
c;      common/clseg/ lseg,intrel,lchar
c;*
c;      isize=0
c;      if(filnam.eq.' ') call tmpnm(luni,filnam)
c;* here for standart fortran i/o
c;      inquire (file=filnam, opened=openfl, exist=exstfl)
c;* here if file exists
c;      if(exstfl) then
c;        if(openfl) return
c;        stat='OLD'
c;      else
c;        stat='NEW'
c;      end if
c;      iun(luni)=luni
c;      nam(luni)=filnam
c;      ipos(luni)=0
c;      iostat(luni)=0
c;      icatf(luni)=icat
c;      if(icat.ne.0) then
c;c.....permanent file
c;          open (luni, file=filnam, status=stat, access='DIRECT',
c;     :           recl=8*lseg, form='UNFORMATTED',
c;     :           iostat=ierr, err=999)
c;      else
c;c.....scratch file
c;          open (luni, status='SCRATCH',access='DIRECT',recl=8*lseg,
c;     :           form='UNFORMATTED',iostat=ierr, err=999)
c;      end if
c;      return
c;* here if open error occurs
c;999    lenfil=index(filnam,' ')-1
c;       write(6,1000) filnam(1:lenfil),ierr
c;1000   format(' %FATAL-ERROR IN ASSGN: cannot open ',(a),', iostat: ',
c;     :        i4,/)
c;       stop
c;c
c;      entry finit
c;c.....closes all files
c;      la=1
c;      le=30
c;      goto 43
c;c
c;      entry closf(luni)
c;c.....closes luni
c;      la=luni
c;      le=luni
c;43    do 44 l=la,le
c;      if(iun(l).eq.0) goto 44
c;      close(l)
c;      icatf(l)=0
c;      nam(l)=' '
c;      iostat(l)=0
c;      ipos(l)=1
c;      iun(l)=0
c;44    continue
c;      return
c;      entry chklun(luni,iflag,name)
c;c.....checks status of unit
c;      name=' '
c;      iflag=0
c;      if(iun(luni).eq.0) return
c;      iflag=icatf(luni)
c;      name=nam(luni)
c;      return
c;      end
c;* ---------------------------------------------------------------
c;      subroutine rdabsf(luni,a,l,iword)
c;*
c;* author: h.j. werner
c;*
c;* read of l double precision words into a. iword is zero adjusted adress
c;* on file. l and iword should be multiple of lseg, otherwise unefficient
c;      implicit double precision (a-h,o-z)
c;      character*14 nam
c;      common /cobuf/ lbuf,ibuf(1)
c;      common/disc/ ipos(98),iun(98),iostat(98),icatf(98),nam(98)
c;      common/clseg/ lseg
c;      dimension a(1),buf(lbuf)
c;      if(lseg.gt.lbuf) stop 'lbuf too small in rdabsf'
c;      ibl=iword/lseg+1
c;      m=iword-(ibl-1)*lseg
c;      n=0
c;      if(m.eq.0) goto 20
c;      n=min0(lseg-m,l)
c;      read(luni,rec=ibl) (buf(i),i=1,lseg)
c;      ibl=ibl+1
c;      do 15 i=1,n
c;  15  a(i)=buf(m+i)
c;  20  l1=l-n
c;      if(l1.eq.0) return
c;      lbl=(l1-1)/lseg+1
c;      do 25 ib=1,lbl
c;      read(luni,rec=ibl) (a(n+i),i=1,min0(l1,lseg))
c;      ibl=ibl+1
c;      l1=l1-lseg
c;  25  n=n+lseg
c;      return
c;c
c;      entry wrabsf(luni,a,l,iword)
c;c
c;      if(lseg.gt.lbuf) stop 'lbuf too small in wrabsf'
c;      ibl=iword/lseg+1
c;      m=iword-(ibl-1)*lseg
c;      n=0
c;      if(m.eq.0) goto 40
c;c.....first address not on sector boundary
c;      n=min0(lseg-m,l)
c;      if(ibl.le.ipos(luni))
c;     1   read(luni,rec=ibl) (buf(i),i=1,lseg)
c;      do 30 i=1,n
c;  30  buf(m+i)=a(i)
c;      write(luni,rec=ibl) (buf(i),i=1,lseg)
c;      ipos(luni)=max(ipos(luni),ibl)
c;      ibl=ibl+1
c;c.....write from next sector boundary
c;  40  l1=l-n
c;      if(l1.eq.0) return
c;      lbl=l1/lseg
c;      do 45 ib=1,lbl
c;      write(luni,rec=ibl) (a(n+i),i=1,lseg)
c;      ipos(luni)=max(ipos(luni),ibl)
c;      ibl=ibl+1
c;      l1=l1-lseg
c;45    n=n+lseg
c;      if(n.eq.l) return
c;c.....last address not at end of sector
c;      if(ibl.le.ipos(luni))
c;     1   read(luni,rec=ibl) (buf(i),i=1,lseg)
c;      do 50 i=1,l-n
c;50    buf(i)=a(n+i)
c;      write(luni,rec=ibl) (buf(i),i=1,lseg)
c;      ipos(luni)=max(ipos(luni),ibl)
c;      return
c;c
c;      entry fwait(luni)
c;      return
c;c
c;      end
cend
cstart none
c;      subroutine tmpnm(ifil,name)
c;c.....should generate a unique filename
c;      character*14 name
c;c     integer getpid
c;      id=0
c;c     id=getpid()
c;      write(name,10) ifil,id
c;10    format('T',i2.2,i5.5,'.TMP  ')
c;      return
c;      end
cend
* ---------------------------------------------------------------
cstart unix
*********************************************************************
*                                                                   *
*             utility programs for unix systems                     *
*                                                                   *
*********************************************************************
*                     routines included:                            *
*                                                                   *
*   1. tmpnam   generates unique temp file names                    *
*   2. timdat   returns time, date and machine type as char. strings*
*   3. exit     stops calculation                                   *
*   4. assgn    assigns direct access files                         *
*               these routines need c-routines in hiunix.c          *
*         entries: finit (closes all files)                         *
*                  closf (closes one file)                          *
*                  chklun (checks if unit is open)                  *
*                  fwait (waits for completion of i/o)              *
*                                                                   *
*********************************************************************
      subroutine tmpnm(ifil,name)
      character*14 name
      integer getpid
      id=getpid()
      write(name,10) ifil,id
10    format('t',i2.2,i5.5,'.tmp  ')
      return
      end
      subroutine assgn(luni,name,lenn,icat)
      implicit double precision (a-h,o-z)
      character*(*) name
      character*14 nam,blank
      character*15 namx
      common/disc/ ipos(98),iun(98),iostat(98),icatf(98),nam(98)
      common/clseg/ lseg
      data lendef/2000/,blank/'              '/
      if (iun(luni).ne.0) return
      if(name.eq.blank) call tmpnm(luni,name)
      if(lenn.eq.0) lenn=lendef
      ll=len(name)
      if(ll.gt.14) name(15:)=' '
      ll=min(ll,14)
      nam(luni)=name
      namx=name(1:ll)
      icat1=index(name,'.TMP')
      if(icat.eq.0.or.icat1.ne.0) then
      icatf(luni)=-1
c temp file
      call openc(luni,namx,lenn,0)
      else
      icatf(luni)=1
c catalog file
      call openc(luni,namx,lenn,1)
      end if
      goto 80
80    iun(luni)=luni
      return
c
      entry finit
c close all files
      do 90 i=1,98
      ipos(i)=1
      iun(i)=0
      nam(i)=' '
      iostat(i)=0
      icatf(i)=0
90    continue
      return
c
      entry closf(luni)
c close file luni
      la=luni
      le=luni
      do 44 l=la,le
      if(iun(l).eq.0) goto 44
      call closec(l)
      icatf(l)=0
      nam(l)=' '
      iostat(l)=0
      ipos(l)=1
      iun(l)=0
44    continue
      return
      entry chklun(luni,iflag,name)
      name=' '
      iflag=0
      if(iun(luni).eq.0) return
      iflag=icatf(luni)
      name=nam(luni)
      return
      end
cend
