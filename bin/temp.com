#!/bin/sh

# HIBRIDON makeconfig $Revision: 99.1 $ --- set up default CONFIG file
# the PATH variable should be set to include the hibridon bin directory

bin=`hibriddir`/bin
PATH=$PATH:$bin
export PATH
cd $bin/..

no_opt='hinput.f himain.f hiversion.f';

if test ! -r $bin/machine.exe ; then echo makemachine ; $bin/makemachine; fi;
machine=`machine.exe`
gmalloc=0 fortran=f77 copt=-O opt=-O debug=-g cnopt=" " noopt=" " static=" " \
i64=" " linkopt=" " profile="-p" libraries=" " seek=.001 speed=200000 cc=cc
if test "$machine" = "unix unix-sun" ; then
  fortran="f77 -Nx500" linkopt="-L/usr/5lib/" gmalloc=1
  opt="-fast"
  machine="$machine unix-noblas"
fi

if test "$machine" = "unix unix-ibm" ; then 
  compiler_level=0.0
  syslevel=`uname -vr |awk '{printf "%d.%d\n",$2,$1}'`
  if test $syslevel = '4.1';then
    compiler_level=`lslpp -h -qcOu xlfcmp 2>/dev/null| awk -F: '/COMMIT/ {print $3}' | head -1 | awk -F. '{printf "%d.%d\n", $1, $2}'`
  fi
  if test $syslevel = '3.2';then 
    compiler_level=`lslpp -h -qcOu xlfcmp.obj 2>/dev/null|awk -F:  '/COMPLETE/{print $4}'|head -1|awk -F. '{printf "%d.%d\n", $1, $2}'`
  fi
  opt=-O; 
  if test $compiler_level = '2.3' -o $compiler_level = '3.2' ; then opt=-O3 ; fi
  if test $compiler_level = '3.1' ; then opt=-O2 ; fi
  fortran="xlf -NA16384 -NQ20000 -qcharlen=1024 -qdpc=e -qxlf77=leadzero -qextname" 
  debug=-g ;
  copt="-O3"
  static="-qsave"
  gmalloc=0
  machine="unix unix-aix"
  extra='# AIX System Level: '$syslevel', xlf Compiler Level:  '$compiler_level
  if test -r /usr/lib/libesslp2.a;
    then machine="unix unix-ibm" libraries="$libraries -lesslp2" opt="$opt -qarch=pwr2";fi
  if test -r /usr/lib/libesslp1.a;
    then machine="unix unix-ibm" libraries="$libraries -lesslp1";fi
  if test -r /usr/lib/libessl.a;
    then machine="unix unix-ibm" libraries="$libraries -lessl";fi
  if test ! -r /usr/lib/libesslp2.a;then 
    extra=$extra'
# If you have a POWER2 machine, consider adding -qarch=pwr2 to the FOPT options'
  fi
fi
if test "$machine" = "unix unix-aix" ; then
   libraries="$libraries -lblas"
fi
if test "$machine" = "unix unix-sequent" ; then
  fortran="fortran -fpa -mp" opt=-O debug=-g ;
fi

if test "$machine" = "unix unix-convex" ; then
  fortran=fc noopt=-no opt=-O2 debug=-g seek=0.0004 speed=800000 gmalloc=1;
fi

if test "$machine" = "unix unix-dec" ; then
  fortran=f77 noopt=-O0 opt="-O -real_size 64" debug=-g gmalloc=1;
  if test -r /usr/lib/libdxml.a;
  then machine="$machine unix-blas3" libraries="$libraries -ldxml";fi
  if test -r /usr/lib/libblas.a;
  then machine="$machine unix-blas3" libraries="$libraries -lblas";fi
fi

if test "$machine" = "unix unix-hp unix-hp300" ; then
  fortran="f77 -U +ffpa" opt="-O -R8" debug=-g gmalloc=1
fi;


if test "$machine" = "cray cray-unicos cray-ymp" ; then
  fortran="cf77 -Wf'-dp '" opt=" " debug="-g"  static="-Wf'-a static'" i64="-Wf'-i 64'" ;
fi

if test "$machine" = "cray cray-unicos cray-c90" ; then
  fortran="cf77 -Wf'-dp '" opt=" " debug="-g"  static="-Wf'-a static'" i64="-Wf'-i 64'" ;
fi

if test "$machine" = "cray cray-unicos cray-2"   ; then
  fortran="cf77 -Wf'-dp -eP '" opt=" " debug="-Wf'-ez'"  static="-Wf'-a static'" i64="-Wf'-i 64'" ;
fi

if test "$machine" = "unix unix-nec" ; then
  fortran="f77sx -pvctl noassume" opt=-O debug="-g -NO -Ni -Nv" noopt="-NO -Nv" linkopt="-l -t";
fi

if test "$machine" = "unix unix-iris" ; then
  fortran="f77 -mips2 -static" opt="-O2 -r8";
  machine="$machine unix-blas3"
  libraries="-lblas "
fi

if test "$machine" = "unix unix-iris unix-r8000" ; then
  cc="cc -G0"
  fortran="f77 -G0 -OPT:roundoff=3:IEEE_arithmetic=3:fast_sqrt=OFF -TENV:X=3 -static"
  static=-static
  seek=.001 speed=2000000
  i64=-i8
  opt="-O3"
  copt="-O3"
  cnopt="-O0"
  machine="$machine unix-blas3"
  libraries="-lfastm -lblas "
fi
# following for old sgi machines only
# if test "$machine" = "unix unix-iris" ; then fortran="f77 -G 0 -static" opt=-O2; fi;

if test "$machine" = "unix unix-fujitsu" ; then
  fortran="frt -Abe" opt="-Wv,-pa" debug="-Wv, -ad" noopt="-Wv,-sc -Ob" linkopt="-J -t -lm";
  extra='FUJITSU_TYPE=VPX
# If your Fujitsu is VPP not VPX, change the above to VPP,
# uncomment the following, and comment out corresponding stanzas below
# CC=fccpx
# FC="frtpx -Ab"
# FOPT=" "
# FDEBUG="-Wv,-ad"
# CDEBUG="-Wv,-ad"
# FNOPT="-Wv,-sc -Ob"
# LINKOPT="-t -lm"'
fi
# reduce machine type
machine=`echo $machine|sed "s/unix //"`

# save old CONFIG file if necessary
if test -r CONFIG ; then mv CONFIG CONFIG.BAK; echo old CONFIG saved as CONFIG.BAK; fi

date=`date` version=`cat VERSION`
cat <<! >CONFIG
# HIBRIDON CONFIG generated at $date for release $version
#               for architecture $machine

$extra
# Machine type ..
MACH="$machine"
# Compilers ..
CC="$cc"
FC="$fortran"
# C defines
# nb	SEEK is disk seek speed in seconds
#	SPEED is disk transfer rate in 8 bytes/second
CDEF="-DSEEK=$seek -DSPEED=$speed"
# compiler optimisation
COPT="$copt"
FOPT="$opt"
# compiler explicit no optimisation
CNOPT="$cnopt"
FNOPT="$noopt"
# compiler debug flag
CDEBUG="$debug"
FDEBUG="$debug"
# static variables
FSTATIC="$static"
# 64-bit integers (deferred implementation)
FI64="$i64"
# profiling
FPROFILE="$profile"
CPROFILE="$profile"
# additional libraries and link options
LIBS="$libraries"
LINKOPT="$linkopt"
# files in src/ to be compiled without optimization
NO_OPT="$no_opt"
# use malloc replacement?
GMALLOC="$gmalloc"
!

echo "CONFIG file created (`hibriddir`/CONFIG)"
echo "  please inspect and edit if necessary before proceeding"
export FC
