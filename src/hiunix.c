/* hiunix.c $Revision: 2000.6 $ */
/* hibridon c files adapted from molpro2002.3 */
#include "machines.h"

#if defined(cray) || defined(linux) || defined(__sun) || defined(__uxp__)
#include <stdlib.h>
#endif
#include <errno.h>

#ifdef	NAME_LU
#define IZLOC izloc_
#define MPPSETMEM mppsetmem_
#define MCTYPE	mctype_
#define	XQUIT	xquit_
#define	INC_CLOCK	inc_clock_
#define	TRUNC	trunc_
#define	OPENC	openc_
#define	DELETC	deletc_
#define	FWAIT	fwait_
#define	WRABSF	wrabsf_
#define	RDABSF	rdabsf_
#define	CLOSEC	closec_
#define TRUNCATEC truncatec_
#define	SECOND	second_
#define	TIMING	timing_
#define	LEADZ	leadz_
#define	POPCNT	popcnt_
#define GMAINV	gmainv_
#define FMAIN	fmain_
#define	FLUSH6	flush6_
#define	GETPID	getpid_
#define	UNLINK	unlink_
#define	EXIT	exit_
#define	ABORT	abort_
#define	KILLME	killme_
#define	TIMDAT	timdat_
#define	TRAP	trap_
#define WALLCL	wallcl_
#define RAND	rand_
#define SRAND	srand_
#endif
#ifdef	NAME_L
#define IZLOC izloc
#define MPPSETMEM mppsetmem
#define	XQUIT	xquit
#define	INC_CLOCK	inc_clock
#define MCTYPE	mctype
#define	TRUNC	trunc
#define	OPENC	openc
#define	DELETC	deletc
#define	FWAIT	fwait
#define	WRABSF	wrabsf
#define	RDABSF	rdabsf
#define	CLOSEC	closec
#define TRUNCATEC truncatec
#define	SECOND	second
#define	TIMING	timing
#define	LEADZ	leadz
#define	POPCNT	popcnt
#define GMAINV	gmainv
#define FMAIN	fmain
#define	FLUSH6	flush6
#define	GETPID	getpid
#define	UNLINK	unlink
#define	EXIT	exit
#define	ABORT	abort
#define	KILLME	killme
#define	TIMDAT	timdat
#define	TRAP	trap
#define WALLCL	wallcl
#define RAND	rand
#define SRAND	srand
#endif


extern char *getenv();

#ifndef RESERVE
#define RESERVE		2048000
#endif

#define MAXUNIT		49
#define WORT		8

#define UNKNOWN		1
#define SCRATCH		0
#define NEW		2
#define OLD		3

#define L0 (long) 0


FILE_DEFINITION files[MAXUNIT]; int first_of_all=1;
double nwords=0, nio=0;

void MCTYPE(ltype,type) char *type; fortint *ltype;
{
strcpy(type,MACHINE_TYPE);
*ltype = strlen(type);
}
void crashhhh() {
#ifdef TRACEBACK
        TRACEBACK;
#endif
        kill(getpid(),SIGIOT);
}


void crasherr(format,string,d1,d2,d3,d4,d5,d6)
char *format, *string; fortint d1,d2,d3,d4,d5,d6; {
        fprintf(stderr,format,string,d1,d2,d3,d4,d5,d6);
        crashhhh();
}

#ifndef cray_old
#ifdef CLSEG
extern struct { fortint lseg, ispace[4], ltrack, jspace[2], icache, iasync;} CLSEG;
/* common/clseg/lseg,ispace(4),ltrack,jspace(2),icache,iasync
 ltrack is disk blocking in 8-byte units
 if (iasync.ne.0) enable asynchronous i/o
 if (icache.ne.0) disable i/o cache           */
#endif


void f_unit_test(u,ochk) fortint u; int ochk; {
	int i;
	if (first_of_all) { /* initialization of files structures */
		first_of_all=0;
		for (i=0;i<MAXUNIT;i++) files[i].nd=-1; /* flag files as not open */
	}
	if (u>MAXUNIT || u<0) crasherr("Unit number out of range: %s%ld\n","",(long)u,L0,L0);
	if (ochk && files[u].nd < 0) crasherr("Attempted I/O before open on unit %s%ld\n","",(long)u,L0,L0);
}

void OPENC(unit,fname
#ifdef FORTCL_NEXT
,lfname
#endif
,size,status
#ifdef FORTCL_END
,lfname
#endif
) fortint *unit, *size, *status; char *fname;
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc lfname;
#endif

{
	char *cp,*cpp,*oldwd,*filewd,*filedir, fnamec[MAXPATHLEN];
	int i, l, lfn, lnb, flags;
#ifdef __hp9000s800 
        OFFSET LSEEK();
#endif
	
#ifdef PAM6000
	/* special filename structure for PAM/6000 .. character*100 array */
#define PAMPATHLEN 100
	char pampaths[MAXPARIO][PAMPATHLEN];
#endif

	f_unit_test((*unit),0);

/* validate file name and make one up if empty */
#ifdef FORTINTC_DIVIDE
  lfname=lfname/FORTINTC_DIVIDE;
#endif
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
	lfn=lfname;
#else
	lfn=strlen(fname);
#endif
#ifdef USE_FCDLEN
	lfn=_fcdlen(fname);
#endif
	for (lnb=lfn;fname[lnb-1] == ' ' && lnb >0; lnb--) ;
	if (lnb == 1 && *fname==' ') {
		sprintf(fname,"Tmp%d",getpid());
		for (cp=fname+lfn-1;*cp != ' '; cp--) lnb=cp-fname+1;
	}
	strncpy((files[*unit].fn = (char *) malloc(lnb+1)),fname,lnb);
	files[*unit].fn[lnb]=(char)0;

	Getcwd((oldwd=(char *) malloc(MAXPATHLEN+1)),MAXPATHLEN); /* Save current directory */
	/* see if there are special instructions for allocating this stream */
#define FILEDIRBASE "MOLPRO_SCR_"
	filedir=(char *) malloc(strlen(FILEDIRBASE)+3);
	sprintf(filedir,"%s%d",FILEDIRBASE,(int)*unit);
	if ((filewd=getenv(filedir)) == NULL) filewd=oldwd;
#ifdef PAMDEBUG
	fprintf(stderr,"filedir=@%s@, filewd=@%s@, oldwd=@%s@\n",filedir,filewd,oldwd);
#endif
	/* parse directory1:directory2:... */
	files[*unit].nd=0; for (cp=filewd; *cp; cp+=l+1) {
		if ((cpp=(char *)index(cp,':')) == NULL) l=strlen(cp); else l=cpp-cp;
		strncpy((files[*unit].dn[files[*unit].nd]=(char *)malloc(l+1)),cp,l);
		files[*unit].dn[files[*unit].nd][l]=(char) 0;
		if (l == strlen(cp)) l--;
		for(i=0;i<files[*unit].nd;i++)
			if (!strcmp(files[*unit].dn[files[*unit].nd],files[*unit].dn[i])) crasherr("Illegal repeated directory in parallel I/O, %s\n",files[*unit].dn[i],L0,L0,L0);
		files[*unit].nd++;
	}

#ifdef PAMDEBUG
	for (i=0; i<files[*unit].nd; i++) fprintf(stderr,"directory @%s@\n",files[*unit].dn[i]);
#endif
	if (files[*unit].nd < PARIOTHRESH) {
		/* regular sequential I/O */
		chdir(files[*unit].dn[0]);
		flags=O_CREAT|OPENFLAGS;
		if (*status == NEW) flags=flags|O_TRUNC;
		if (*status == OLD) flags=OPENFLAGS;
		if ((files[*unit].fd=OPEN(files[*unit].fn,O_RDWR|flags,0666))==-1)
		if ((files[*unit].fd=OPEN(files[*unit].fn,O_RDONLY|flags,0666))==-1)
			crasherr("openc: Error in opening file %s\n",files[*unit].fn,L0,L0,L0);
		*size=(fortint)(LSEEK(files[*unit].fd,(OFFSET) 0,2)+(OFFSET)511)/(OFFSET)512;
	} else {
		/* parallel I/O */
#ifdef PAM6000
		for (l=0;l<files[*unit].nd;l++) {
			if ((cp=rindex(files[*unit].fn,'/')) == NULL) cp=files[*unit].fn; else cp++;
			sprintf(pampaths[l],"%s/%s",files[*unit].dn[l],cp);
			for (i=strlen(pampaths[l]);i<PAMPATHLEN;i++) pampaths[l][i]=' ';
#ifdef PAMDEBUG
			fprintf(stderr,"%d:%s pampath @%s@\n",*unit,files[*unit].fn,pampaths[l]);
#endif
		}
		i=CLSEG.lseg*WORT;
		pamo(unit,&files[*unit].nd,&i,pampaths,&l,PAMPATHLEN);
		if (l!=0)crasherr("openc: Error in opening PAM file %s, error number=%ld\n",files[*unit].fn,(long)l,L0,L0);
		*size=0;
#else
		crasherr("Unimplemented parallel I/O\n","",L0,L0,L0);
#endif
	}
#ifdef NOKEEPTEMP
#undef KEEPTEMP
#endif
#ifndef KEEPTEMP
	/* ensure temporary files vanish on closing */
	if (*status == SCRATCH)	for (l=0;l<files[*unit].nd;l++) {
		chdir(files[*unit].dn[l]);
		unlink(files[*unit].fn);
	}
#endif

	files[*unit].size= *size;
	files[*unit].addr= -1;
#ifdef convex
	if (CLSEG.iasync) fcntl(files[*unit].fd,F_SETFL,FASIO);
	if (CLSEG.icache) fcntl(files[*unit].fd,F_SETFL,FNCACHE);
#endif
	chdir(oldwd); free(oldwd);
}

void FWAIT(unit) fortint *unit;
{
	f_unit_test(*unit,1);
#ifdef convex
	if (CLSEG.iasync) if (asiostat(files[*unit].fd) < 0) crasherr("fwait: I/O error on unit %s%ld\n","",(long)*unit,L0,L0);
#endif
}

void f_seek(unit,p) fortint *unit, *p;{
	fortint ierr, block;
	f_unit_test(*unit,1);
	if (*p==files[*unit].addr) return;
	if (*p < 0)
	  crasherr("f_seek: illegal negative file offset on file %s, unit %ld, offset=%ld\n",files[*unit].fn,(long)*unit,(long)*p,L0);
	if (sizeof(OFFSET) < 8 && *p >= 268435455)
	  crasherr("f_seek: file offset on file %s, unit %ld overflows maximum; offset=%ld\n",files[*unit].fn,(long)*unit,(long)*p,L0);
	if (sizeof(block) < 8 && *p > 2147483647)
	  crasherr("f_seek: file offset on file %s, unit %ld overflows maximum; offset=%ld\n",files[*unit].fn,(long)*unit,(long)*p,L0);
	if (files[*unit].nd < PARIOTHRESH) {
		/* regular I/O */
		if (LSEEK(files[*unit].fd,(OFFSET)*p * (OFFSET)WORT,0)==-1)
                        crasherr("f_seek: Error in lseek on file %s, unit %ld, at offset %ld\n",files[*unit].fn,(long)*unit,(long)*p,L0);
	} else {
		/* parallel I/O */
#ifdef PAM6000
		block=*p / CLSEG.lseg;
		if (block*CLSEG.lseg != *p)
			crasherr("f_seek: Illegal non-block-aligned call\n","",L0,L0,L0);
		block++;  /* 1st block is 1 not 0 */
		pamp(unit,&block,&ierr);
		if (ierr) crasherr("f_seek: Error in seeking PAM file %s at offset=%ld, block=%ld, error number=%ld\n",files[*unit].fn,(long)*p,(long) block,(long) ierr);
#else
		crasherr("Unimplemented parallel I/O\n","",L0,L0,L0);
#endif
	}
	files[*unit].addr=*p;
}

void WRABSF(unit,a,l,p) fortint *unit, *l, *p; char *a;
{
	fortint block, ierr;
	f_seek(unit,p);
	if (files[*unit].nd < PARIOTHRESH) {
		/* regular I/O */
		ierr = write (files[*unit].fd,a,*l * WORT) - *l * WORT;
	} else {
		/* parallel I/O */
#ifdef PAM6000
		block = *l / CLSEG.lseg;
		if (block*CLSEG.lseg != *l)
			crasherr("wrabsf: Illegal non-block-aligned call; length=%ld, lseg=%ld\n",(long)*l, (long)CLSEG.lseg,L0);
		pamw(unit,a,&block,&ierr);
#else
		crasherr("Unimplemented parallel I/O\n","",L0,L0,L0);
#endif
	}
	if (ierr) crasherr("wrabsf: Error in writing to file %s (unit %ld), %ld words at word offset %ld\n",files[*unit].fn,(long)*unit,(long)*l,(long)*p);
	nio++; nwords+= *l;
	files[*unit].addr += *l;
}

void RDABSF(unit,a,l,p) fortint *unit, *l, *p; char *a;
{
	fortint block, ierr;
	f_seek(unit,p);
	if (files[*unit].nd < PARIOTHRESH) {
		/* regular I/O */
		ierr = read (files[*unit].fd,a,*l * WORT) - *l * WORT;
	} else {
		/* parallel I/O */
#ifdef PAM6000
		block = *l / CLSEG.lseg;
		if (block*CLSEG.lseg != *l)
			crasherr("rdabsf: Illegal non-block-aligned call\n","",L0,L0,L0);
		pamr(unit,a,&block,&ierr);
#else
		crasherr("Unimplemented parallel I/O\n","",L0,L0,L0);
#endif
	}
        if (*p != 0) {
	if (ierr) crasherr("rdabsf: Error in reading from file %s (unit %ld), %ld words at word offset %ld\n",files[*unit].fn,(long)*unit,(long)*l,(long)*p);
        }
	nio++; nwords+= *l;
	files[*unit].addr += *l;
}


/* mgs stuff */
void TRUNCATEC(unit,length) fortint *unit, *length;
{     
        f_unit_test(*unit,0);
        if (files[*unit].nd < 1) return;
        if (files[*unit].nd < PARIOTHRESH) {
                /* regular I/O */
          if (TRUNCATE(files[*unit].fd,(OFFSET)*length * (OFFSET)WORT))  {
            perror("error in truncatec: \n");
            crasherr("truncatec: Error in truncating file %s\n",files[*unit].fn,L0,L0,L0); 
          }
        }
}
/* end mgs stuff */


void CLOSEC(unit) fortint *unit;
{
	int i; fortint ierr, iscratch;
	f_unit_test(*unit,0);
	if (files[*unit].nd < 1) return;
	if (files[*unit].nd < PARIOTHRESH) {
		/* regular I/O */
		close(files[*unit].fd);
	} else {
		/* parallel I/O */
#ifdef PAM6000
		iscratch=0;
		pamc(unit,&iscratch,&ierr);
		if (ierr != 0) crasherr("closec: Error in closing PAM file %s, error number=%ld\n",files[*unit].fn,(long) ierr,L0,L0);
#else
		crasherr("Unimplemented parallel I/O\n","",L0,L0,L0);
#endif
	}
	free(files[*unit].fn);
	for (i=0;i<files[*unit].nd;i++) free(files[*unit].dn[i]);
	files[*unit].fd=files[*unit].addr=files[*unit].size=0;
}
#endif
/*###################  timing etc ###############*/

/*
      this routine will return the user-time in second
*/
static clock_t it; struct tms itt;

#ifndef cray
double SECOND()

{ it=times(&itt);
return (double) itt.tms_utime / (double) TIMER_TICK;
}

void TIMING(cpu,sys,io) double *cpu,*sys,*io;

{ it=times(&itt);
*cpu=(double) itt.tms_utime / (double) TIMER_TICK; *sys=(double) itt.tms_stime / (double) TIMER_TICK;
*io=(double) nio * SEEK + (double) nwords / SPEED;
}

#endif
/* subroutine timdat(tim,dat,mach)
   character*(*) tim,dat,mach */
#ifdef HAS_UTSNAME
#include <sys/utsname.h>
#endif
#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN	255
#endif
struct timeval tpp;
char *tval, *ctime();
void TIMDAT(tim
#ifdef FORTCL_NEXT
,ltim
#endif
,dat
#ifdef FORTCL_NEXT
,ldat
#endif
,mach
#ifdef FORTCL_NEXT
,lmach
#endif
#ifdef FORTCL_END
,ltim,ldat,lmach
#endif
) char *tim, *dat, *mach;
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc ltim,ldat,lmach;
#endif
{	
	char *os, *osver, *host, *hardw, buffer[256];
#if defined(IRIX64) || defined(linux)
	time_t timdum;
#endif
#ifdef _AIX
#define	nibmtab	35
	char ibmtabm[nibmtab][3]={
		"31"    ,"35"    ,"30"    ,"10"    ,"18"    ,"14"    ,"1C"    ,
		"5C"    ,"20"    ,"2E"    ,"11"    ,"37"    ,"38"    ,"41"    ,
		"63"    ,"64"    ,"66"    ,"34"    ,"43"    ,"67"    ,"75"    ,
		"76"    ,"77"    ,"78"    ,"80"    ,"71"    ,"70"    ,"47"    ,
		"49"    ,"46"    ,"58"    ,"82"    ,"57"    ,"48"    ,"43"};
	char ibmtabv[nibmtab][7]={
		"320"   ,"320H"  ,"520"   ,"530"   ,"530H"  ,"540"   ,"550"   ,
		"560"   ,"930"   ,"950"   ,"540"   ,"340"   ,"350"   ,"220"   ,
		"970"   ,"980"   ,"580"   ,"520H"  ,"M20"   ,"570"   ,"370"   ,
		"360"   ,"350"   ,"315"   ,"990"   ,"580H"  ,"590"   ,"230"   ,
		"250-80","250-66","380"   ,"R24"   ,"390"   ,"C10"   ,"M20"};
	int i;
#endif
#ifdef HAS_UTSNAME
	struct utsname n;
#else
#endif
#ifdef FORTINTC_DIVIDE
  ltim=ltim/FORTINTC_DIVIDE;ldat=ldat/FORTINTC_DIVIDE;lmach=lmach/FORTINTC_DIVIDE;
#endif
#ifdef SX
        gettimeofday(&tpp);
#else
        gettimeofday(&tpp,0);
#endif
#if defined(IRIX64) || defined(linux)
	timdum = tpp.tv_sec;
	tval=ctime(&timdum);
#else
	tval=ctime(&tpp.tv_sec);
#endif
	sprintf(tim,"%8.8s",tval+11);
	sprintf(dat,"%3.3s-%3.3s-%2.2s",tval+7,tval+4,tval+22);
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
	for (os=tim+strlen(tim);os<tim+ltim;os++) *os=' ';
	for (os=dat+strlen(dat);os<dat+ldat;os++) *os=' ';
#endif
#ifdef HAS_UTSNAME
	uname(&n);
	os=n.sysname; osver=n.release; host=n.nodename; hardw=n.machine;
#ifdef _CRAYMPP
	hardw="Cray MPP";
#endif
#ifdef _AIX
	osver=(char *) malloc(strlen(n.release)+strlen(n.version)+1);
	strcpy(osver,n.version);strcat(osver,".");strcat(osver,n.release);
	for(i=0;i<nibmtab;i++) if (!strncmp(ibmtabm[i],n.machine+8,2)) hardw=ibmtabv[i];
#endif
#ifdef convex
	osver=n.version;hardw=n.release;
#endif
#else
	os="Unix"; osver=index(MACHINE_TYPE,'-')+1;hardw=(char *) 0;
	gethostname((host=(char *) malloc(MAXHOSTNAMELEN+1)),MAXHOSTNAMELEN);
#endif
	strcpy(buffer,os); strcat(buffer,"-"); strcat(buffer,osver);
	strcat(buffer,"/");strcat(buffer,host);
	if (hardw) {strcat(buffer,"(");strcat(buffer,hardw);strcat(buffer,")");}
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
	strncpy(mach,buffer,lmach);
	for (os=mach+strlen(mach);os<mach+lmach;os++) *os=' ';
#else
	strcpy(mach,buffer);
#endif
}


double WALLCL()
{
/* On Solbourne (Sun?) must link with sys V library (-L/usr/5lib/ ) */
it=times(&itt);
return (double) it / (double) TIMER_TICK;
}

#ifdef SYSTEm
SYSTEm(s) char *s; { return system(s); }
#endif
#ifdef GETENv
void GETENv(s
#ifdef FORTCL_NEXT
,l
#endif
,rs
#ifdef FORTCL_NEXT
,rl
#endif
#ifdef FORTCL_END
,l,rl
#endif
) char *rs, *s;
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc rl, l;
#else
Error: program not completed!
#endif
{ char *bla, *p1, *p2; int i;

#ifdef FORTINTC_DIVIDE
  rl=rl/FORTINTC_DIVIDE;l=l/FORTINTC_DIVIDE;
#endif
  bla=(char *) calloc(l+1,1);
  for (p1=bla, p2=s, i=0; (i<l && *p2 != ' '); i++) *p1++ = *p2++; *p1=(char)0;
  for (p1=rs; p1<rs+rl; p1++) *p1=' ';
  if ((p2=getenv(bla)) != NULL)
    for (p1=rs; *p2 && p1<rs+rl; p1++, p2++) *p1 = *p2;
  free(bla);
}
#endif

#ifdef PUT_ENV
void PUT_ENV(s
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
,l
#endif
) char *s;
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc l;
#else
Error: program not completed!
#endif
{ char *bla, *p1, *p2; int i;

#ifdef FORTINTC_DIVIDE
  l=l/FORTINTC_DIVIDE;
#endif
  bla=(char *) calloc((int)l+1,1);
  for (p1=bla, p2=s, i=0; (i<(int)l && *p2 != ' '); i++) *p1++ = *p2++; *p1=(char)0;
  putenv(bla);
#if !defined(_AIX) && !defined(__sun) && !defined(linux) && !defined(__alpha)
    free(bla);
#endif
}
#endif

/*###################  bit handling ###############*/

#ifndef cray
static unsigned char tab1[]=
      { 8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3,
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      };

static unsigned char tab2[]=	
      { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3,
	2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4,
	2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5,
	4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
	3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
	4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
	3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5,
	4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
	3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
	5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
      };

fortint LEADZ(pp) register unsigned char *pp;
{ register fortint i,j; register unsigned char *p=pp;
#ifndef REVERSE_BYTE
i=tab1[*p++]; if (i<8) return i;
i+=tab1[*p++]; if (i<16) return i;
i+=tab1[*p++]; if (i<24) return i;
#ifdef I64
i+=tab1[*p++]; if (i<32) return i;
i+=tab1[*p++]; if (i<40) return i;
i+=tab1[*p++]; if (i<48) return i;
i+=tab1[*p++]; if (i<56) return i;
#endif
#else
#ifdef I64
p+=7;
i=tab1[*p--]; if (i<8) return i;
i+=tab1[*p--]; if (i<16) return i;
i+=tab1[*p--]; if (i<24) return i;
i+=tab1[*p--]; if (i<32) return i;
i+=tab1[*p--]; if (i<40) return i;
i+=tab1[*p--]; if (i<48) return i;
i+=tab1[*p--]; if (i<56) return i;
#else
p+=3;
i=tab1[*p--]; if (i<8) return i;
i+=tab1[*p--]; if (i<16) return i;
i+=tab1[*p--]; if (i<24) return i;
#endif
#endif
i+=tab1[*p];
return i;
}

fortint POPCNT(p) register unsigned char *p;

{ register fortint i;

#ifdef I64
i=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p++];
i+=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p];
#else
i=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p++]; i+=tab2[*p];
#endif
return i;
}
#ifdef MVBITs
#ifndef REVERSE_BYTE
static unsigned int bits[]={
          1,          2,          4,          8,
         16,         32,         64,        128,
        256,        512,       1024,       2048,
       4096,       8192,      16384,      32768,
      65536,     131072,     262144,     524288,
    1048576,    2097152,    4194304,    8388608,
   16777216,   33554432,   67108864,  134217728,
  268435456,  536870912, 1073741824, 2147483648};
#else
static unsigned int bits[]={
2147483648, 1073741824,  536870912,  268435456,
  134217728,   67108864,   33554432,   16777216,
    8388608,    4194304,    2097152,    1048576,
     524288,     262144,     131072,      65536,
      32768,      16384,       8192,       4096,
       2048,       1024,        512,        256,
        128,         64,         32,         16,
          8,          4,          2,          1};
#endif

void MVBITs(m,i,len,n,j) int *m,*i,*len,*n,*j;
{
	int k;
	for (k=0;k<*len; k++){
        	*n=(bits[*i+k] & (*m)) ? *n | bits[*j+k] : *n & (~bits[*j+k]) ;
	}
}

#endif

#endif

/*###################  memory allocation  ###############*/

#ifdef MA_ALLOC
#define EXT_INT
#include "macdecls.h"
#include "matypes.h"
#define MA_GLOBALSIZE 1000000      /* Default Global memory allocation on heap 
                                      needs to determined dynamically */
static int ma_init_ok = 0;
static Integer ma_hndl = -1;
static Pointer ma_ptr = (Pointer) NULL;

#define GLOBAL_MEM_DEFAULT   "1000000"
#define STACK_MEM_DEFAULT    "1000000"

void MPPSETMEM()
{
	char	str[1024], *ptr;
	int	gmem_req, smem_req;

	strcpy(str,GLOBAL_MEM_DEFAULT);
	if ((ptr=getenv("MOLPRO_GLOBAL_MEM")))
		strcpy(str,ptr);
	gmem_req = atoi(str);
	strcpy(str,STACK_MEM_DEFAULT);
	if ((ptr=getenv("MOLPRO_STACK_MEM")))
		strcpy(str,ptr);
	smem_req = atoi(str);
#ifdef MA_ALLOC_DEBUG
	printf("\n\nInit MA\n Stack: %d\n Global: %d\n",smem_req, gmem_req);
#endif
	if (MA_init(MT_DBL,(Integer) smem_req,(Integer) gmem_req)==MA_FALSE)
		crasherr("Cannot initialize MA\n","",L0,L0,L0);
}
#endif

/*  #ifndef cray  */

#ifdef AIX_SHM
static int n_shm_seg=0;
static int shmid[AIX_SHM_SEG_MAX],nseg, nbyte;
static char *start_address, *address, *seg_address[AIX_SHM_SEG_MAX];
#endif

#ifndef MA_ALLOC
static char *pointer;
static char *reserve=0;
void GMAINV(q,ibase,n) double *q; fortint *ibase, *n;
{
	double *qq; int reserving=0, i, lenseg; fortint inttmp;
#ifdef AIX_SHM
	struct shmid_ds buffer;
	if (*n > AIX_SHM) {
		/* used shared memory routines for large request */
		nbyte = (*n+3)*WORT;
		n_shm_seg = (nbyte-1) / SHMLBA + 1;
		if (n_shm_seg > AIX_SHM_SEG_MAX) crasherr("Request for %s%ld shared memory segments is too many\n","",(long) n_shm_seg,L0,L0);
		start_address=0;
		for (i=0;i<n_shm_seg;i++) {
			lenseg=nbyte-i*SHMLBA; if (lenseg > SHMLBA) lenseg = SHMLBA;
			shmid[i] = shmget(IPC_PRIVATE,lenseg,(IPC_CREAT|S_IRUSR|S_IWUSR));
			if (shmid[i] == -1) {
				if (errno==EINVAL) crasherr("Memory allocation: invalid identifier in shmget\n","",L0,L0,L0);
				if (errno==EACCES) crasherr("Memory allocation: permission denied in shmget\n","",L0,L0,L0);
				if (errno==ENOMEM) crasherr("Memory allocation: insufficient paging space in shmget\n","",L0,L0,L0);
				crasherr("Memory allocation: unknown error in shmget\n","",L0,L0,L0);
			}
			seg_address[i]=(char *) shmat(shmid[i],start_address,SHMLBA);
			if (shmctl(shmid[i],IPC_RMID,&buffer)) {
				fprintf(stderr,"gmainv: warning: error in removing shared memory ID\n");
				perror("gmainv");
			}
			if (seg_address[i] < 0) {
				fprintf(stderr,"Memory allocation: shmat failure %d\n",seg_address[i]);
				perror(" ");
				crasherr("Memory allocation problem\n","",L0,L0,L0);
			}
			if (!i) pointer=seg_address[i];
			start_address = seg_address[i]+lenseg;
		}
  	qq=(double *) pointer;
  	*ibase=qq-q+1;
	return;
	}
#endif
#ifdef MA_ALLOC
/*
	Use MA to allocate analogous to malloc() call)
*/
/* 
    This is replaced with MPPSETMEM()

	if (!ma_init_ok) {
		{
			int	nstack;

                 	nstack = *n + MA_sizeof_overhead(MT_DBL);
			if (MA_init(MT_DBL,nstack,MA_GLOBALSIZE)==MA_FALSE)
				crasherr("Cannot initialize MA\n","",L0,L0,L0);
			++ma_init_ok;
		}
	}

*/
#ifdef debugma
	fprintf(stderr,"gmainv %d\n",*n);
	MA_trace(1);
#endif
	if (MA_push_stack(MT_DBL,(Integer) (*n),"MOLPRO segment",&ma_hndl)==MA_FALSE)
		crasherr("Memory allocation failed with MA\n","",L0,L0,L0);
#ifdef debugma
	fprintf(stderr,"Following request for %d MA returns handle=%d\n",*n,ma_hndl);
#endif
	if (MA_get_pointer(ma_hndl,&ma_ptr)==MA_FALSE)
		crasherr("MA cannot get pointer for block\n","",L0,L0,L0);
#ifdef debugma
	fprintf(stderr,"Following request for %d MA returns handle=%d, pointer=%d\n",*n,ma_hndl,ma_ptr);
#endif
	*ibase = (((double *)ma_ptr) - q) + 1;
#else
/* 
   Default: use malloc()
 */
#ifdef RESERVE
	/* carve out a low memory block for other users of malloc
   	this will ensure that gmainv block is always on top, so long as RESERVE
   	is large enough */
  	if (! reserve) {
		reserving=RESERVE; if (reserving > *n-1) reserving=*n-1;
		reserve = MEMALLOC(reserving);
		MEMALLOC(1); /* barrier */
	}
#endif
	inttmp = (fortint) (*n)+3;
  	pointer = MEMALLOC((*n+3));
  	if (pointer == NULL) {
	  printf("Failure in attempting memory allocation of %ld words (%ld Mbyte)\n",(long) (*n), (long)(*n /(128*1024)));
	  printf("This error has been generated by the operating system,\n");
	  printf("and may be the result of insufficient system memory or paging space\n");
	  printf("In order to avoid the problem in the MOLPRO context,\n");
	  printf("consider also reducing the requested memory through\n");
	  printf("the MEMORY input command, or the -m command line option\n");
	  crasherr("gmainv1 failure to allocate %s%ld\n","",(long) inttmp,L0,L0);
	}
  	qq=(double *) pointer;
  	*ibase=qq-q+1;
  	if (reserving) free(reserve);
#endif
}

void FMAIN(q,n) double *q; int *n;
{
	int i,ierr;
	*n=0;
#ifdef AIX_SHM
	ierr=0;
	for (i=0;i<AIX_SHM_SEG_MAX;i++) {
		if (shmid[i]) {
			if (shmdt((char *)seg_address[i]) < 0) ierr++;
			shmid[i] = 0;
		}
		if (ierr) crasherr("Problems found in cleaning up shared memory;\nuse the ipcs and iprm commands to clean up\n","",L0,L0,L0);
	}
	if (n_shm_seg) { n_shm_seg = 0; return;}
#endif
#ifdef MA_ALLOC
#ifdef debugma
	fprintf(stderr,"fmain: MA handle=%d\n",ma_hndl);
#endif
	if (ma_hndl<0)
		crasherr("No MA block to free?\n","",L0,L0,L0);
	MA_pop_stack(ma_hndl);
	ma_hndl = -1;
	ma_ptr = (Pointer) NULL;
#else
	free(pointer);
#endif
}
#endif
/* #endif   */

/*###################  error handling  ###############*/

/* flush all fortran units */
void FLUSH6() {fflush(stdout);fflush(stderr);}

void handler(sig)
{char string[10];
FLUSH6();
/* If the signal was SIGUSR1 then just flush buffers and continue */
if (sig == SIGUSR1) {signal(SIGUSR1,handler); return;}
#ifdef __alpha
/*if (sig == SIGFPE) {
          fprintf(stdout,">>> Floatingpoint exception encountered. <<<\n");
          fprintf(stderr,"WARNING: floatingpoint exception. Continuing. \n");
	  signal(SIGFPE,handler); return;} */
#endif
sprintf(string,"Signal %d",sig);perror(string);
#ifdef TRACEBACK
TRACEBACK;
#endif
/*signal(SIGTERM,SIG_DFL); kill(getpid(),SIGTERM);*/
/*signal(sig,SIG_DFL); kill (getpid(),sig);*/
exit(sig);
}
void TRAP()
{
/* if environment variable NOTRAP is set, do not set up traps */
/* ul environment variable   TRAP is set, do not set up traps */
#ifdef NOTRAP
      if (! getenv("TRAP")) return;
#else
      if (getenv("NOTRAP")) return;
#endif
      signal(SIGHUP,handler);
      signal(SIGINT,handler);
      signal(SIGQUIT,handler);
      signal(SIGILL,handler);
      signal(SIGTRAP,handler);
      signal(SIGIOT,handler);
#ifdef SIGEMT
      signal(SIGEMT,handler);
#endif
      signal(SIGFPE,handler);
      signal(SIGBUS,handler);
      signal(SIGSEGV,handler);
#ifdef SIGSYS
      signal(SIGSYS,handler);
#endif
      signal(SIGPIPE,handler);
#ifndef NOALRM
      signal(SIGALRM,handler);
#endif
      signal(SIGTERM,handler);
#ifdef SIGXCPU
      signal(SIGXCPU,handler);
#endif
#ifdef SIGXFSZ
      signal(SIGXFSZ,handler);
#endif
#ifdef SIGVTALRM
      signal(SIGVTALRM,handler);
#endif
#ifdef SIGPROF
      signal(SIGPROF,handler);
#endif
#ifdef SIGLOST
      signal(SIGLOST,handler);
#endif
#ifdef SIGUSR1
      signal(SIGUSR1,handler);
#endif
#ifdef SIGUSR2
      signal(SIGUSR2,handler);
#endif
#ifdef SIGTTIN
      signal(SIGTTIN,handler);
#endif
#ifdef SIGTTOU
      signal(SIGTTOU,handler);
#endif
#ifdef SIGABRT
      signal(SIGABRT,handler);
#endif
}

void KILLME() {
#ifndef CRASH_SIGNAL
#define CRASH_SIGNAL SIGINT
#endif
  signal(CRASH_SIGNAL,SIG_DFL);
#ifdef TRACEBACK
  TRACEBACK;
#endif
  kill(getpid(),CRASH_SIGNAL);
  exit(1);
}

#ifndef	NAME_L
#ifndef cray
int GETPID() { return getpid();}
int UNLINK(file) char *file; { return unlink(file);}

void ABORT(a
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
,la
#endif
) char *a;
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc la;
#endif
{
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
#ifdef FORTINTC_DIVIDE
  la=la/FORTINTC_DIVIDE;
#endif
  a[la]=(char) NULL; 
#endif
  fprintf(stdout,"%s\n",a);
  fflush(stdout);
  fprintf(stderr,"%s\n",a);
  fflush(stderr);
  KILLME();
}

#endif
#endif


#ifndef cray_old
void DELETC(fname
#ifdef FORTCL_NEXT
,lfname
#endif
#ifdef FORTCL_END
,lfname
#endif
) 
char *fname; 
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
fortintc lfname;
#endif


{ char *cp; int l, lfn;

#ifdef FORTINTC_DIVIDE
  lfname=lfname/FORTINTC_DIVIDE;
#endif
#if defined(FORTCL_END) || defined(FORTCL_NEXT)
	lfn=lfname;
#else
	lfn=strlen(fname);
#endif
#ifdef USE_FCDLEN
	lfn=_fcdlen(fname);
#endif
  fname[lfn]=(char) NULL; cp = (char *) (fname+lfn-1);
  if (*cp==' ') { while (*(cp--)==' '); *(cp+2) = (char) NULL; }
  if (!access(fname,F_OK)) unlink(fname);
}
#endif

/******************************************************************************/
/*                                                                            */
/*                                 I Z L O C                                  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* This routine converts the pointer to a memory address into an integer      */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Author:  Per-Olof Widmark                                                  */
/*          IBM Sweden                                                        */
/* Written: September 1993                                                    */
/*                                                                            */
/******************************************************************************/


fortint IZLOC(i)
 int i;
{ 
  return(i);
}

#if defined(cray) || defined(linux) || defined(__sun) || defined(__uxp__)
double RAND () {
	return ( (double) rand()) / ( (double) RAND_MAX);
}

void SRAND (seed) fortint *seed; {
	srand(*seed);
}
#endif


void XQUIT(iRc)   

fortint *iRc;   

{

 int i;

 i=(int)*iRc;
 (void)printf("STOP %d \n",i);
 exit(i);

}


fortint INC_CLOCK()

{

  fortint cpu;

   cpu=sysconf((fortint) _SC_CLK_TCK);
/*
  cpu=CLOCKS_PER_SEC;
*/

  return(cpu);
}

#if defined(pgf90)
/* fix for portland/blas compatibility */

void MAIN__() { printf("MAIN__ should not be called\n"); }
#endif

