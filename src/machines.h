/* machines.h $Revision: 2002.6 Patch(2002.6): opteron2 $ */
/* MOLPRO machine characteristics for C programs */
#ifndef __MACHINES_H__
#define __MACHINES_H__
/* disk parameters for printing i/o times
... override at compile time with -DSEEK=.....  etc.! */


#ifndef SEEK
#define SEEK		0.001		/* average seektime in seconds */
#endif
#ifndef SPEED
#define SPEED		307250		/* speed in words per second */
#endif

/* only HP-UX10.20 and later is supported */
#ifdef __hpux
#define NOGETWD
#define MVBITs	MVBITS
#define HAS_UTSNAME
#ifdef ppu
#define NAME_LU
#else
#define NAME_L
#endif
#ifdef I64
#ifdef __ia64
#define MACHINE_TYPE	"unix unix-i8 unix-hp unix-ia64"
#else
#define MACHINE_TYPE	"unix unix-i8 unix-hp unix-hp800"
#endif
#ifndef FORTINT
#define FORTINT long long
#endif
#ifndef FORTINTC
#define FORTINTC long
#endif
#else
#ifdef __ia64
#define MACHINE_TYPE	"unix unix-i4 unix-hp unix-ia64"
#else
#define MACHINE_TYPE	"unix unix-i4 unix-hp unix-hp800"
#endif
#define OPEN open64
#define LSEEK lseek64
#define TRUNCATE ftruncate64
#define OFFSET          off64_t
#endif
#define FORTCL_END
#endif

#ifdef _AIX
#define NAME_LU
#define OPEN open
#ifdef AIX42
#undef OPEN
#define OPEN open64
#define LSEEK lseek64
#define TRUNCATE ftruncate64
#define OFFSET          off64_t
#endif
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-ibm"
#ifndef FORTINT
#define FORTINT long long
#endif
#ifndef FORTINTC
#define FORTINTC long
#endif
#else
#define MACHINE_TYPE	"unix unix-i4 unix-ibm"
#ifndef AIX_SHM
#define AIX_SHM 16777216 /* threshold for using shared memory segments */
#endif
#ifndef AIX_SHM_SEG_MAX
#define AIX_SHM_SEG_MAX 10
#endif
#endif
#include <memory.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <sys/errno.h>
#include <sys/mode.h>
#define HAS_UTSNAME
#define FORTCL_END
#ifdef PAM6000
#define CLSEG clseg_
#endif
#ifndef RESERVE
#define RESERVE 10000
#endif
#endif

#ifdef ultrix
#define NAME_LU
#define REVERSE_BYTE
#define MACHINE_TYPE	"unix unix-i4 unix-dec"
#define HAS_UTSNAME
#define FORTCL_END
#endif

#ifdef __alpha
#ifndef linux
#include <sys/types.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#include <db.h>
#endif
#define NAME_LU
#define REVERSE_BYTE
#define HAS_UTSNAME
#define FORTCL_END
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-dec"
#ifndef FORTINT
#define FORTINT long 
#endif
#ifndef FORTINTC
#define FORTINTC int
#endif
#else
#define MACHINE_TYPE	"unix unix-i4 unix-dec"
#endif
#define MEMALLOC(n)	(char *)malloc(n * sizeof(double))
#ifndef FORTINT
#define FORTINT int
#endif
#endif

#ifdef __convexc__
#define convex
#endif
#ifdef convex
#define NOMALLOCH
/* char *malloc(); */
#define TRACEBACK traceback_()
#define	NAME_LU
#define CLSEG	_clseg_
#define MACHINE_TYPE	"unix unix-i4 unix-convex"
#define HAS_UTSNAME
#define FORTCL_END
#define	OPENFLAGS	O_LARGEFILE
#define OFFSET		off64_t
#define	LSEEK		lseek64
#define TRUNCATE 	ftruncate64
#endif

#ifdef sgi
#include <string.h>
#include <strings.h>
#define	NAME_LU
#define HAS_UTSNAME
#define FORTCL_END
#ifndef FORTINTC
#define FORTINTC long
#endif
#ifndef FORTINT
#ifdef I64 /* integer*8 */
#define FORTINT long long
#else
#define FORTINT int
#endif /*I64 */
#endif /* FORTINT */

#ifdef R10000
#ifdef I64 /* integer*8 */
#define MACHINE_TYPE    "unix unix-i8 unix-iris unix-r10000"
#else
#define MACHINE_TYPE    "unix unix-i4 unix-iris unix-r10000"
#endif /*I64 */
#endif /*R10000*/
#ifdef R8000
#ifdef I64 /* integer*8 */
#define MACHINE_TYPE    "unix unix-i8 unix-iris unix-r8000"
#else
#define MACHINE_TYPE    "unix unix-i4 unix-iris unix-r8000"
#endif /*I64 */
#endif /*R8000*/
#ifndef MACHINE_TYPE
#ifdef I64 /* integer*8 */
#define MACHINE_TYPE    "unix unix-i8 unix-iris"
#else
#define MACHINE_TYPE    "unix unix-i4 unix-iris"
#endif /*I64 */
#endif /*MACHINE_TYPE*/

#if (_MIPS_ISA > 2)
#define OFFSET off64_t
#define LSEEK lseek64
#define TRUNCATE ftruncate64
#endif /* mips3 or higher */

#endif /*sgi*/

#ifdef	linux

#ifdef LARGEFILES
#define _USE_LARGEFILE64
#define _LARGEFILE64_SOURCE
#define OPENFLAGS O_LARGEFILE
#define OPEN open64
#define LSEEK lseek64
#define TRUNCATE ftruncate64
#define OFFSET __off64_t
#else  /* LARGEFILES */
#define TRUNCATE ftruncate
#endif /* LARGEFILES */

#include <string.h>
#define	NAME_LU
#define FORTCL_END
#define REVERSE_BYTE
#ifdef __alpha
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-linux unix-linux-alpha"
#else  /* I64 */
#define MACHINE_TYPE	"unix unix-i4 unix-linux unix-linux-alpha"
#endif /* I64 */
#endif /* __alpha */
#ifdef __ia64
#define MACHINE_TYPE	"unix unix-i8 unix-linux unix-linux-ia64"
#define MEMALLOC(n)       calloc(n , sizeof(double))
#ifdef I64
#ifndef FORTINT 
#define FORTINT long long
#endif /* FORTINT */
#ifndef FORTINTC
#define FORTINTC long
#endif /* FORTINTC */
#endif /* I64 */
#endif /* __ia64 */
#ifdef __x86_64__
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-linux unix-linux-x86_64"
#ifndef FORTINT
#define FORTINT long long
#endif /* FORTINT */
#ifndef FORTINTC
#define FORTINTC long
#endif /* FORTINTC */
#else  /* I64 */
#define MACHINE_TYPE	"unix unix-i4 unix-linux unix-linux-x86_64"
#endif /* I64 */
#endif /* __x86_64__ */
#if ! defined (__ia64) && ! defined (__alpha) && ! defined (__x86_64__)
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-linux"
#ifndef FORTINT
#define FORTINT long long
#endif /* FORTINT */
#ifndef FORTINTC
#define FORTINTC long
#endif /* FORTINTC */
#else  /* I64 */
#define MACHINE_TYPE	"unix unix-i4 unix-linux"
#endif /* I64 */
#endif /* ! defined (__ia64) && ! defined (__alpha)  && ! defined (__x86_64__)*/
#define HAS_UTSNAME
#endif /* linux */

#ifdef __ppc__ 

#define TRUNCATE ftruncate
#include <string.h>
#define NAME_LU
#define FORTCL_END
#define REVERSE_BYTE
#define MACHINE_TYPE    "unix unix-i4 unix-darwin"
#define HAS_UTSNAME
#ifndef FORTINT
#define FORTINT int
#endif 

#endif

#ifdef	sun
/* the next 4 lines might not be applicable on all OS */
#define OPEN open64
#define LSEEK lseek64
#define TRUNCATE ftruncate64
#define OFFSET          off64_t
#ifdef	hal
/* This seems to be needed for Solaris on HAL */
#define NOGETWD
#define	MVBITs	mvbits_
#define	NAME_LU
#define FORTCL_END
/* this doesn't seem to be needed any more
typedef long clock_t;
*/
#define MACHINE_TYPE	"unix unix-i4 unix-sun unix-hal"
#define HAS_UTSNAME
#define KEEPTEMP
#else
/* This seems to be needed for Solaris */
#define NOGETWD
#define	MVBITs	mvbits_
#define	NAME_LU
#define FORTCL_END
/* this doesn't seem to be needed any more
typedef long clock_t;
*/
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-sun"
/* following apparently not needed with -xarch=v9a
#ifndef FORTINT
#define FORTINT long long
#endif
#ifndef FORTINTC
#define FORTINTC long
#endif
*/
#else
#define MACHINE_TYPE	"unix unix-i4 unix-sun"
#endif
#define HAS_UTSNAME
#endif
#endif

#ifdef	sequent
#define NOMALLOCH
#define NOTRAP
#define	MEMALLOC(n)	shmalloc(n * sizeof(double))
char *shmalloc();
#define	free shfree
#define	NAME_LU
#define FORTCL_END
#define REVERSE_BYTE
typedef long clock_t;
#define MACHINE_TYPE	"unix unix-i4 unix-sequent"
#endif

#ifdef SX
typedef char int8_t;
typedef short int16_t;
/* typedef int int32_t; */
#define KEEPTEMP
#define NAME_LU
#define FORTCL_END
#define HAS_UTSNAME
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-nec"
#ifndef FORTINT
#define FORTINT long
#endif
#ifndef FORTINTC
#define FORTINTC long
#endif
#else
#define MACHINE_TYPE	"unix unix-i4 unix-nec"
#ifndef FORTINT
#define FORTINT int
#endif
#ifndef FORTINTC
#define FORTINTC int
#endif
#endif
#endif

#ifdef _UNICOS
#include <bind/bitypes.h>
#define KEEPTEMP
#include <sys/machd.h>
#define HAS_UTSNAME
#ifdef CRAY2
#define MACHINE_TYPE	"cray cray-unicos cray-2"
#define TRACEBACK tracebk()
#define USE_FCDLEN
#else
#ifdef _CRAYMPP
#define MACHINE_TYPE	"cray cray-unicos cray-mpp"
#define FORTCL_NEXT
/* Cray T3D, but not T3E, multiplies character string lengths by 8 */
#ifdef _CRAYT3D
#define FORTINTC_DIVIDE	8
#endif
#define GETENv	GETENV
#else
#ifdef _CRAYIEEE
/* This is IEEE T90.  Amazing combo of c90 and t3e stuff --- and do not define USE_FCDLEN! */
#define MACHINE_TYPE	"cray cray-unicos cray-t90"
#define FORTCL_NEXT
#define FORTINTC_DIVIDE	8
#else
#define USE_FCDLEN
#ifdef  CRAYC90
#define MACHINE_TYPE	"cray cray-unicos cray-c90"
#else
#define MACHINE_TYPE	"cray cray-unicos cray-ymp"
#endif
#endif
#endif
#define TRACEBACK trbk()
#endif
#endif

#ifdef __uxp__
typedef char int8_t;
typedef short int16_t;
typedef int int32_t;
#define NOGETWD
#define	MVBITs	mvbits_
#define NAME_LU
#define FORTCL_END
#define HAS_UTSNAME
#define _LLTYPES
#define TRACEBACK errtra_()
#ifdef I64
#define MACHINE_TYPE    "unix unix-i8 unix-fujitsu"
#define FORTINT long
#define FORTINTC int
#else
#define MACHINE_TYPE	"unix unix-i4 unix-fujitsu"
#define FORTINT int
#define FORTINTC int
#endif
#define NOINDEX
#define NOTRAP
#endif

#ifndef MACHINE_TYPE
/*#error "Sorry, machines.h hasn't been set up for your machine"*/
#define MACHINE_TYPE	"unknown"
#endif

/* fortran integer type */
#ifndef FORTINT
#define FORTINT	long
#endif
typedef FORTINT fortint ;
/* fortran character string length type */
#ifndef FORTINTC
#define FORTINTC FORTINT
#endif
typedef FORTINTC fortintc ;

#ifndef OPENFLAGS
#define OPENFLAGS	0
#endif
#ifndef OFFSET
#define OFFSET		off_t
#endif
#ifndef OPEN
#define OPEN 		open
#endif
#ifndef LSEEK
#define LSEEK		lseek
#endif
#ifndef TRUNCATE
#define TRUNCATE	ftruncate
#endif

#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>


#ifndef MEMALLOC
#define MEMALLOC(n)	malloc(n * sizeof(double))
/*#define MEMALLOC(n)	calloc(n , sizeof(double))*/
#endif
#ifndef NOMALLOCH
#ifdef __ppc__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif
#endif

/* sysv/ucb nonsense */
#ifdef NOINDEX
#define index	strchr
#define rindex	strrchr
#endif

/* SJM 1/3/1999 got rif of getwd */ 
/*
  #ifdef NOGETWD
  char *getcwd();
*/
#define Getcwd(a,b)	getcwd(a,b);
/*
  #else
  char *getwd();
  #define Getcwd(a,b)	getwd(a);
  #endif
*/

#include <sys/time.h>
#include <sys/times.h>
#include <sys/param.h>


#ifdef _CRAYMPP
#define TIMER_TICK	CLK_TCK
#endif
#ifndef TIMER_TICK
#ifndef HZ
#define TIMER_TICK      60
#else
#define TIMER_TICK      HZ
#endif
#endif


#ifdef FORTCL_END
#define FORTCL
#endif
#ifdef FORTCL_NEXT
#define FORTCL
#endif

#ifndef MAXPATHLEN
#define MAXPATHLEN	1024
#endif

#define MAXPARIO	8
typedef struct {
	int fd; OFFSET addr, size; int nd;
	char *fn, *dn[MAXPARIO];
} FILE_DEFINITION;
/* PARIOTHRESH=1: parallel I/O stuff used for one thread too; otherwise set 2 */
#ifndef PARIOTHRESH
#define PARIOTHRESH	2
#endif


#ifdef USE_FCDLEN
#include <fortran.h>
#endif
#endif

#if defined(cray) || defined(linux) || defined(__sun) || defined(__uxp__) || defined(__ppc__)
#include <stdlib.h>
#endif
#include <errno.h>
