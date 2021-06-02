*comdeck machines.h
/* machines.h $Revision: 95.1 $ */
/* MOLPRO machine characteristics for C programs */

/* disk parameters for printing i/o times
... override at compile time with -DSEEK=.....  etc.! */
#ifndef SEEK
#define SEEK		0.001		/* average seektime in seconds */
#endif
#ifndef SPEED
#define SPEED		307250		/* speed in words per second */
#endif


#ifdef hpux
#define NOGETWD
#define MVBITs	MVBITS
#define HAS_UTSNAME
#ifdef hp9000s300
#define MACHINE_TYPE	"unix unix-i4 unix-hp unix-hp300"
#define FORTCL_END
#define SYSTEm	SYSTEM
#define GETENv	GETENV
#endif
#ifdef hp9000s800
#define NAME_L
#define MACHINE_TYPE	"unix unix-i4 unix-hp unix-hp800"
#define FORTCL_END
#endif
#endif


#ifdef _AIX
#define NAME_L
#ifdef I64
#define MACHINE_TYPE	"unix unix-i8 unix-ibm"
#else
#define MACHINE_TYPE	"unix unix-i4 unix-ibm"
#endif
#ifndef AIX_SHM
#define AIX_SHM 16777216 /* threshold for using shared memory segments */
#endif
#ifndef AIX_SHM_SEG_MAX
#define AIX_SHM_SEG_MAX 10
#endif
#include <memory.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <sys/errno.h>
#include <sys/mode.h>
#define HAS_UTSNAME
#define FORTCL_END
#ifdef PAM6000
#define CLSEG clseg
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
#define NAME_LU
#define REVERSE_BYTE
#define HAS_UTSNAME
#define FORTCL_END
#define MACHINE_TYPE	"unix unix-i4 unix-dec"
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
char *malloc();
#define TRACEBACK traceback_()
#define	NAME_LU
#define CLSEG	_clseg_
#define MACHINE_TYPE	"unix unix-i4 unix-convex"
#define HAS_UTSNAME
#define FORTCL_END
#define	OPENFLAGS	O_LARGEFILE
#define OFFSET		off64_t
#define	LSEEK		lseek64
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
#if (_MIPS_SZLONG == 64)
#ifdef I64
#ifndef FORTINT
#define FORTINT long
#endif
#define MACHINE_TYPE	"unix unix-i8 unix-iris unix-r8000"
#else
#ifndef FORTINT
#define FORTINT int
#endif
#define MACHINE_TYPE	"unix unix-i4 unix-iris unix-r8000"
#endif
#define IRIX64
#else
#ifndef FORTINT
#define FORTINT int
#endif
#define MACHINE_TYPE	"unix unix-i4 unix-iris"
#endif
#endif

#ifdef	sun
/* This seems to be needed for Solaris */
#define NOGETWD
#define	MVBITs	mvbits_
#define	NAME_LU
#define FORTCL_END
typedef long clock_t;
#define MACHINE_TYPE	"unix unix-i4 unix-sun"
#define HAS_UTSNAME
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
#define NAME_LU
typedef long clock_t;
#define FORTCL_END
#define HAS_UTSNAME
#define MACHINE_TYPE	"unix unix-i4 unix-nec"
#endif

#ifdef _UNICOS
#include <sys/machd.h>
#define HAS_UTSNAME
#ifdef CRAY2
#define MACHINE_TYPE	"cray cray-unicos cray-2"
#define TRACEBACK tracebk()
#else
#ifdef _CRAYMPP
#define MACHINE_TYPE	"cray cray-unicos cray-mpp"
#define FORTCL_NEXT
#define FORTINTC_DIVIDE	8
#define GETENv	GETENV
#else
#ifdef  CRAYC90
#define MACHINE_TYPE	"cray cray-unicos cray-c90"
#else
#define MACHINE_TYPE	"cray cray-unicos cray-ymp"
#endif
#endif
#define TRACEBACK trbk()
#endif
#endif

#ifdef __uxp__
#define NOGETWD
#define	MVBITs	mvbits_
#define NAME_LU
#define FORTCL_END
#define HAS_UTSNAME
#define TRACEBACK errtra_()
#define MACHINE_TYPE	"unix unix-i4 unix-fujitsu"
#define NOINDEX
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
#ifndef LSEEK
#define LSEEK		lseek
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
#include <malloc.h>
#endif

/* sysv/ucb nonsense */
#ifdef NOINDEX
#define index	strchr
#define rindex	strrchr
#endif

#ifdef NOGETWD
char *getcwd();
#define Getcwd(a,b)	getcwd(a,b);
#else
char *getwd();
#define Getcwd(a,b)	getwd(a);
#endif

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


