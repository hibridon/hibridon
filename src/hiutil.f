**********************************************************************
*                                                                    *
*                          utility library                           *
*                                                                    *
**********************************************************************
*                         routines included:                         *
*  1. vaxhlp/upcase,getwor     subroutine to mimic vaxhelp
*  2. parse      parser for analyzing command strings                *
*  3. getval     assigns numeric command strings to variables        *
*  4. upper      converts string to upper case string                *
*  5. lower      converts string tp lower case string                *
*  6. dater      subroutine to return calendar date and clock time   *
*  7. mtime      routine to return cpu and wall clock time in seconds*
*  8. second     function to return current cpu time
*  9. gettim     converts seconds to time string: hhh:mm:ss          *
* 9a. timdat     returns time and date
* 10. gennam     routine to generate filenames                       *
* 11. factlg     calculates factorials used in xf3j, xf6j, xf9j      *
* 12. xf3j       function, returns value of 9j-symbol                *
* 13. xf6j       function, returns value of 6j-symbol                *
* 14. xf9j       function, returns value of 9j-symbol                *
* 15. f3j0,f6j   3j and 6j routines for two molecule program
* 16. intairy    evaluate s integrals of airy functions              *
* 17.  cheby     evaluate expansion in chebyshev polynomials         *
**********************************************************************
* NB cstart ultrix for fortran rather than c utilities
**********************************************************************
      subroutine vaxhlp(line1)
* latest revision 19-aug-1991
      implicit character*10(z)
      character*(*) line1
      common /coipar/ jtot1, jtot2, jtotd, jlpar, nerg,numax,numin,nud,
     :                lscreen, iprint
cstart cray
c;* return for these machines
c;        write (6,2)
c;2       format (' Help not installed on this machine')
c;        return
c;        end
cend
cstart unix mac
c simulates vax vms help command
      parameter (inp=5,iout=6,ihlpmain=51,ihlpalt=52)
      parameter (maxlev=9,maxsub=50)
      character*(*) helpdir,helptail
cend
      include "common/parhlp"
cstart unix max
      character*80 helpfile,line,lin,lin2
      character*1 pat
cend
cstart unix mac
      logical exist
cend
cstart unix mac
      dimension zreq(maxlev),lreq(maxlev),zsub(maxsub)
c...  by default, looks in file help. if that is not available, then
c     tries file name of first key in the tree.
      helpfile=helpdir//'hibrid'//helptail
cend
cstart unix-convex mac
c;      open (ihlpmain,file=helpfile,err=299,readonly,status='old')
cend
cstart unix-sun unix-hp unix-dec unix-ibm unix-aix unix-iris
c      write(6,*) 'open ',helpfile
      open (ihlpmain,file=helpfile,err=299,status='unknown')
cend
cstart unix mac
      zalt=' '
      level=0
      line = line1
      if (line.ne.' ') goto 1111
1     if (level.lt.0) goto 99
      levini = level
      zreq(level+1)='Topic?'
      lreq(level+1)=6
      lin2=' '
      write (iout,'(1x)')
      ll = 0
      do 1112 i=1,level+1
      lin2=lin2(1:ll+1)//zreq(i)(1:lreq(i))
      ll=lenstr(lin2)
      if (ll.ge.72) write (iout,'(a)') lin2(1:ll)
      if (ll.ge.72) lin2=' '
      if (ll.ge.72) ll = 0
1112  continue
      if (ll.ge.1) write (iout,'(1x,a)') lin2(1:ll)
      lin2=' '
c      write (iout,111) (zreq(i)(1:lreq(i)),i=1,level+1)
c111   format(/1x,a,8(' ',a,:))
      read (inp,'(a)',end=99) line
1111  call upcase(line)
c     if (line(1:3).eq.'---'.or.line(1:3).eq.'end') goto 99
c...  parse the key
      if (line.eq.' ') then
c..     null input; backspace one level
        level = level-1
        goto 1
      else if (line.eq.'?') then
c...    repeat last request
      else
      ipos=1
      do 10 nreqq=1,maxlev-level
      level = level+1
      call getwor(line,ipos,zreq(level))
      lreq(level)=lenstr(zreq(level))
      if (line(ipos:).eq.' ') goto 12
10    continue
12    continue
      end if
c
c...  open and position file. tries help first
      isear=1
      if (levini.eq.0) then
c...    top level, so need to open/switch file possibly
        ihlp = ihlpmain
        rewind ihlp
        line=' '
211     read (ihlp,'(a)',end=291) line(1:12)
        if (line(1:1).ne.'1') goto 211
        call upcase(line)
        lin=line
        if (level.eq.0) goto 40
        if (line(3:lreq(1)+2).ne.zreq(1)) goto 211
c...    key found in help
        ipos=3
        call getwor(line,ipos,zreq(1))
        lreq(1)=lenstr(zreq(1))
        isear=2
        goto 215
c...    key not found in help. attempt to open alternate file.
291     ihlp=ihlpalt
        ilev=1
        if (zalt.ne.zreq(1)) then
          if (zalt.ne.' ') close(ihlp)
          zalt=zreq(1)
cend
cstart unix mac
          helpfile='ls '//helpdir//zalt(1:lenstr(zalt))//'*'//helptail
     >     //' >helphelphelphelp'
          call system(helpfile)
          inquire (file='helphelphelphelp',exist=exist)
          if (.not.exist) goto 29
          open (ihlpalt,file='helphelphelphelp',status='old')
          ifound=0
          read (ihlpalt,'(a)',end=295) helpfile
          ifound=1
295       close (ihlpalt,status='delete')
          if (ifound.eq.0) goto 29
cend
cstart unix-convex
c;          open (ihlpalt,file='helphelphelphelp',status='old')
c;          open (ihlp,file=helpfile,err=29,readonly,status='old')
cend
cstart unix-sun unix-hp unix-dec unix-ibm unix-aix unix-iris
          open (ihlp,file=helpfile,err=29,status='unknown')
cend
cstart unix mac
        end if
      end if
      rewind ihlp
215   do 20 ilev=isear,level
      pat=char(ilev+ichar('0'))
21    read (ihlp,'(a)',end=29) line(1:12)
      if (line(1:1).ne.pat) goto 21
      call upcase(line)
      if (line(3:lreq(ilev)+2).ne.zreq(ilev)) goto 21
      ipos=3
      call getwor(line,ipos,zreq(ilev))
      lreq(ilev)=lenstr(zreq(ilev))
20    continue
      goto 28
29    ll = 0
      lin2=' '
      do 1122 i=1,ilev
      lin2=lin2(1:ll+1)//zreq(i)(1:lreq(i))
      ll=lenstr(lin2)
1122  continue
      write (iout,293) lin2(1:ll)
      lin2=' '
293   format(' key not found:',a)
      level=ilev-1
      goto 1
28    continue
c
c...  print the text from the file
c     write (iout,'(1x,a)') line(3:lenstr(line))
c     write (iout,'(10(1x,a))') (zreq(i)(1:lreq(i)),i=1,level)
      ll = 0
      lin2=' '
      do 3087 i=1,level
      lin2 = lin2(1:ll+1)//zreq(i)
3087  ll = lenstr(lin2)
      write (iout,'(a)') lin2(1:ll)
      write (iout,'(1x)')
      iscreen=1
      do 30 irec=1,9999
      read (ihlp,'(a)',end=45) lin
      if (lin(1:1).ge.'0'.and.lin(1:1).le.'9') goto 40
      if (iscreen.ge.lscreen) then
        write (iout,'(/,'' Press RETURN to continue'')')
        read (inp,'(a)') line
        level = level-1
        if (line.ne.' ') goto 1111
        level = level+1
        iscreen=0
      end if
      write (iout,'(1x,a)') lin(1:lenstr(lin))
30    iscreen = iscreen+1
c...  print sub level keys
40    zmatch=char(ichar('0')+level+1)
      level = level-1
      if (lin(1:1).ne.zmatch) goto 1
      level = level+1
      nsub=1
      zsub(1)=lin(3:)
      do 41 irec=1,9999
      read (ihlp,'(a1,1x,a10)',end=49) zmat,zsubk
      if (zmat.lt.'0'.or.zmat.gt.'9') goto 41
      if (zmat.lt.zmatch) goto 49
      if (zmat.eq.zmatch) then
        nsub=nsub+1
        if (nsub.gt.maxsub) goto 49
        zsub(nsub) = zsubk
      end if
41    continue
49    write (iout,48) (zsub(i),i=1,nsub)
48    format(/' Additional information:'/30(/7a11))
      goto 1
c...  no sub topics
45    level = level-1
      goto 1
c
99    continue
      close (ihlpmain)
      if (zalt.ne.' ') close(ihlpalt)
      return
299   write (iout,'('' Primary help file '',a,'' is missing'')')
     >  helpfile(1:lenstr(helpfile))
      end
cend
* ----------------------------------------------------------------
      subroutine upcase(l)
      character*(*) l
      ishift = ichar('A')-ichar('a')
      ia = ichar('a')
      iz = ichar('z')
      do 1 i=1,len(l)
      j=ichar(l(i:i))
      if (j.ge.ia.and.j.le.iz) l(i:i)=char(j+ishift)
1     continue
      return
      end
      subroutine getwor (line,ipos,crd)
      character*80 line
      character*(*) crd
      i2=81
      if(line(ipos:80).eq.' ') then
        crd=' '
        ipos=81
        return
      end if
      do 5 i1=ipos,80
      if (line(i1:i1).ne.' ') goto 6
5     continue
6     ii=index(line(i1+1:80),' ')
      if(ii.eq.0) goto 10
      i2=ii+i1
10    ipos = i2
      crd = line(i1:i2-1)
      return
      end
* ----------------------------------------------------------------
      subroutine parse(line,i,code,lcod)
*   current revision date: 9-apr-90
*   returns string between delimiters "," or ";" in code(1:lcod)
*   on input, iabs(i) points to first character to be searched in line
*   on output, i points to first character after next delimiter
*   if remainder of line is blank, i=0 is returned
*   if i.lt.0 on entry, first delimiter may also be '='
*   if i.eq.0 on entry, blank is returned and i unchanged
      character*(*) line,code
      code=' '
      lcod=0
      if(i.eq.0) return
      i1=iabs(i)
      k1=index(line(i1:),',')
      k2=index(line(i1:),';')
      k3=index(line(i1:),'=')
      k=k1
      if(k1.eq.0.or.(k2.ne.0.and.k2.lt.k1)) k=k2
      if(i.lt.0.and.k3.ne.0.and.(k3.lt.k.or.k.eq.0)) k=k3
      i=i1+k
      if(k.eq.0) then
        i=0
c       k=index(line(i1:),';')
        if(k.eq.0) k=len(line)-i1+2
      end if
      lcod=k-1
      if(lcod.le.0) then
        code=' '
        lcod=0
        return
      end if
      i2=i1+lcod-1
      do 10 k1=i1,i2
10    if(line(k1:k1).ne.' ') goto 20
      code=' '
      lcod=0
      return
20    do 30 k2=i2,k1,-1
30    if(line(k2:k2).ne.' ') goto 40
40    lcod=k2-k1+1
      code=line(k1:k2)
      return
      end
* ----------------------------------------------------------------
      subroutine getval(code,clist,nlist,i,val)
*   current revision date: 23-sept-87
*   searches strings in clist to match code
*   returns associated value in val and position in clist in i
*   values may be specified as integers or reals
*   if(code=t(rue)  is specified, val=1 is returned
*   if(code=f(alse) is specified, val=0 is returned
*   delimiter between code and value is equal sign or blank
      implicit double precision (a-h,o-z)
      character*(*) code,clist(1)
      character*80 line
      l=-1
      if(nlist.eq.0) goto 30
      l=index(code,'=')
      l=l-1
      if(l.le.0) then
        write(6,5)
 5      format(' equal sign missing in specification')
        i=-1
        return
      end if
      do 10 i=1,nlist
10    if(code(1:l).eq.clist(i)(1:l)) goto 30
      i=0
      return
30    l1=l+2
      l2=len(code)
      val=0
      if(l2.lt.l1) return
      do 40 j=l1,l2
40    if(code(j:j).ne.' ') goto 50
      val=0
      return
50    if(code(j:j).eq.'T') then
        val=1.0
        return
      end if
      if(code(j:j).eq.'F') then
        val=0.0
        return
      end if
      k=index(code(j:l2),'.')
      line=code
      if(k.eq.0) then
        do 60 l=l2,j,-1
60      if(line(l:l).ne.' ') goto 70
70      line(l+1:)='.'
      end if
      read(line(j:),80) val
80    format(f40.5)
      return
      end
* ----------------------------------------------------------------
      subroutine upper(line)
      character*(*) line
      l = len(line)
      iac = ichar('A')
      ial = ichar('a')
      izl = ichar('z')
      idf = iac-ial
      do 10 i = 1,l
      icr = ichar(line(i:i))
      if(icr.lt.ial.or.icr.gt.izl) goto 10
      icr = icr+idf
      line(i:i) = char(icr)
10    continue
      return
      end
* ----------------------------------------------------------------
      subroutine lower(line)
*  subroutine to convert the character string 'line'
*  to lower case
*  only alphabetic characters are changed
*  current revision date: 23-sept-87
      character*(*) line
      integer l, i, idf, icr
      l=len(line)
      iac=ichar('A')
      ial=ichar('a')
      izc=ichar('Z')
*  idf is the offset between upper case and lower case characters
      idf=ial-iac
      do 10 i=1,l
        icr=ichar(line(i:i))
        if(icr.ge.iac.and.icr.le.izc) then
          icr=icr+idf
          line(i:i)=char(icr)
        end if
10    continue
      return
      end
* -----------------------------------------------------------------
      subroutine dater (cdate)
*  subroutine to return calendar date and clock time
*  current revision date: 1-jan-1996 by mha
*  -----------------------------------------------------------------
*  variables in call list
*
*  on return:  cdate is a character string containing the date
*              in format dd-mmm-yy followed by the time in
*              format hh:mm:ss
*  -----------------------------------------------------------------
      character*20 cdate
*  this is specifically adapted to the cray xmp-48 at the san diego
*  supercomputer center running under ctss and fortlib
cstart unix mac
      character*10 ctime
      character*10 cday
      character*32 mach
      call timdat(ctime,cday,mach)
      cdate = cday//'  '//ctime//'   '
cend
      return
      end
c-------------------------------------------------------------------
      subroutine mtime(t1,t2)
*
*  subroutine to return cpu and wall clock time in seconds
*  current revision date: 3-jun-1991 by mha
*
*  -----------------------------------------------------------------
*  variables in call list
*
*    t1:     on return:  contains cpu clock time in seconds
*    t2:     on return:  contains wall clock time in seconds
*
*  -----------------------------------------------------------------
        implicit double precision(a-h,o-z)
cstart mac
c;      t1=second()
c;      t2=t1
cend
cstart unix cray
      double precision cpu,sys,tio
      call timing(cpu,sys,tio)
      t1=cpu
      t2=sys+tio+cpu
cend
cstart ultrix-dec
c;      common /codec/ ttim(2)
c;      real et, etime, secnds
c;      dimension et(2)
c;      t1=etime(et(1),et(2))
c;      tt2=secnds(0.0)
c;      delt=tt2-ttim(2)
c;      if (delt .lt. 0.0) delt=delt+86400.
c;      ttim(1)=ttim(1)+delt
c;      ttim(2)=tt2
c;      t2=ttim(1)
cend
      return
      end
*------------------------------------------------------------------------
cstart mac
c;      function second()
c;      double precision second
c;      second = long(362)/60.0d0
c;      end
cend
cstart ultrix-dec
c;      double precision function second()
c;      implicit double precision (a-h,o-z)
c;      real et(2)
c;      fe=etime(et(1),et(2))
c;      second=fe
c;      return
c;      end
cend
*------------------------------------------------------------------------
      subroutine gettim(sec,time)
*
*  author: b. follmeg
*  current revision date: 23-sept-87
*
*  subroutine to convert timing in seconds into string which contains
*  time in format hour:min:sec:millisec (000:00:00.00)
*  on input:  sec   -> seconds
*  on output: time  -> string containing time (must have been dimensioned
*                      as character*12 at least)
      integer ihour, imin, isec, imilli
      double precision sec
      character*(*) time
      isecn=100*sec
      imilli=mod(isecn,100)
      isecn=isecn/100
      isec=mod(isecn,60)
      isecn=isecn/60
      imin=mod(isecn,60)
      ihour=isecn/60
      write(time,10) ihour,imin,isec,imilli
10    format(i3.2,':',i2.2,':',i2.2,'.',i2.2)
      return
      end
* --------------------------------------------------------------------
cstart none
c;ccstart mac
c;      subroutine timdat(tim,dat,mach)
c;c.....should return time, date and machine type
c;      character *(*) tim,dat,mach
c;      call time(tim)
c;      call date(dat)
c;      mach='MAC IIci'
c;      return
c;      end
c;ccend
c;ccstart ultrix-dec
c;      subroutine timdat(tim,dat,mach)
c;      character*8 tim
c;      character*9 dat
c;      character*19 mach
c;      mach='DECStation  ULTRIX'
c;      call time(tim)
c;      call date(dat)
c;      return
c;      end
c;ccend
c;ccstart unix-dec unix-sun
c;      subroutine timdat(tim,dat,mach)
c;      character *(*) tim,dat,mach
c;      mach='DECStation  ULTRIX'
c;      call timec(tim)
c;      call datec(dat)
c;      l=len(tim)
c;      if(l.gt.8) tim(9:l)=' '
c;      l=len(dat)
c;      if(l.gt.9) dat(10:l)=' '
c;      return
c;      end
c;ccend
cend
cstart unix-irix
c;      subroutine timdat(tim,dat,mach)
c;      character *(*) tim,dat,mach
c;      mach='SGI IRIX'
c;      call timec(tim)
c;      call datec(dat)
c;      l=len(tim)
c;      if(l.gt.8) tim(9:l)=' '
c;      l=len(dat)
c;      if(l.gt.9) dat(10:l)=' '
c;      return
c;      end
cend
cstart none
c;ccstart unix-hp unix-convex unix-ibm unix-aix
c;      subroutine timdat(tim,dat,mach)
c;      character *(*) tim,dat,mach
c;      character*9 a(7)
c;      integer uname
c;      call timec(tim)
c;      call datec(dat)
c;      l=len(tim)
c;      if(l.gt.8) tim(9:l)=' '
c;      l=len(dat)
c;      if(l.gt.9) dat(10:l)=' '
c;      i=uname(a)
c;      mach=' '
c;      mach(1:5)=a(1)(1:5)
c;      mach(7:9)=a(3)(1:3)
c;      mach(11:18)=a(5)(1:8)
c;      return
c;      end
c;ccend
cend
      subroutine gennam(xname,jobnam,ifil,xtens,ln)
*
* subroutine to generate filnames
* on input: jobnam ->  filename
*           ifil ->    if > 0 number is concatened to filename (e.g. filnam1)
*                      if = 0 no effect
*           xtens ->   extension
*
* on output: xname ->  generated filnam (e.g. filnam1.ext)
*            ln    ->  length of xname
*
* current revision date: 4-jun-1991
*
      character*(*) xname,jobnam,xtens
      character*1 dot
      include "common/comdot"
      xname=jobnam
      i=index(xname,dot)
      if(i.eq.0) i=index(xname,' ')
      if(ifil.eq.0) then
        ln=i-1
      else if(ifil.lt.10) then
        write(xname(i:),10) ifil
10      format(i1)
        ln=i
      else
        write(xname(i:),20) ifil
20      format(i2)
        ln=i+1
      end if
      xname(ln+1:)=dot
      xname(ln+2:)=xtens
      ln=ln+len(xtens)+1
      return
      end
* ----------------------------------------------------------------
      subroutine factlg(n)
c
c...subroutine to compute factorials needed in xf3j,xf6j and xf9j
c   functions.
c   this subroutine must be called before calling any of the
c   j-symbols routines
c   the common block cofact must be dimensioned in the main program !
c   on input: n  maximum number  of factorials to be calculated
c
      double precision si
      common /cofact/ si(1)
      si(1) = 0.0d0
      do 5  k = 2,n
        si(k) = si(k-1) + dlog (dble(k-1))
5     continue
      return
      end
      function xf3jm0(a,b,c)
*     program to compute the 3j symbol (a b c / 0 0 0)
*     author:  b. follmeg
*     current revision date: 14-dec-87
*     revised by mha 4-may-1997
*
* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      data zero,two /0.d0, 2.d0/
      x = zero
* check for triangular conditions
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 100
      iabc = nint(a + b + c)
* check for even sum
      ig = iabc / 2
      if (ig*2 .ne. iabc) goto 100
* compute delta
      j1 = iabc - nint(two*a) + 1
      j2 = iabc - nint(two*b) + 1
      j3 = iabc - nint(two*c) + 1
      j4 = iabc + 2
      delta = 0.5d0 * (si(j1) + si(j2) + si(j3) - si(j4))
      j1 = ig + 1
      j2 = nint(ig - a) + 1
      j3 = nint(ig - b) + 1
      j4 = nint(ig - c) + 1
      x = delta + si(j1) - si(j2) - si(j3) - si(j4)
      x = dexp(x)
      ip = nint(ig + c + a - b)
      if ((ip/2)*2 .ne. ip) x = -x
100   xf3jm0 = x
      return
      end
* --------------------------------------------------------------------------
      function xf3j(a,b,c,am,bm,cm)
*     program to compute the 3j symbol (a b c / am bm cm)
*     authors:  t. orlikowski and b. follmeg
*     current revision date: 4-may-1997

* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      data tol /1.d-10/
      data zero,one,two /0.d0, 1.d0,2.d0/
      x = zero
* check for triangular conditions
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 3
      abc = a + b + c
      iabc = nint(abc)
      if (abs(am + bm + cm) .gt. tol) goto 3
      if (abs(am) .gt. a) goto 3
      if (abs(bm) .gt. b) goto 3
      if (abs(cm) .gt. c) goto 3
      if (am.eq.0 .and. bm.eq.0 .and.((iabc/2)*2).ne.iabc) goto 3
      iacbm = nint(a - c + bm)
      ibcam = nint(b - c - am)
      iabmc = nint(a + b - c)
      iamam = nint(a - am)
      ibpbm = nint(b + bm)
      iapam = nint(a + am) + 1
      ibmbm = nint(b - bm) + 1
      icpcm = nint(c + cm) + 1
      icmcm = nint(c - cm) + 1
      minchi = max0(0,ibcam,iacbm)
      maxchi = min0(iabmc,iamam,ibpbm) - minchi
      iabmc = iabmc - minchi
      iaam = iamam - minchi
      ibbm = ibpbm - minchi
      ibcam = minchi - ibcam
      iacbm = minchi - iacbm
      iamam = iamam + 1
      ibpbm = ibpbm + 1
* compute delta
      j1 = iabc - nint(two*a) + 1
      j2 = iabc - nint(two*b) + 1
      j3 = iabc - nint(two*c) + 1
      j4 = iabc + 2
      delta = si(j1) + si(j2) + si(j3) - si(j4)
      ll=0
      x = one
      if (maxchi.le.0) goto 2
      a1 = dble(ibcam + maxchi + 1)
      a2 = dble(iacbm + maxchi + 1)
      a3 = dble(minchi + maxchi + 1)
      b1 = dble(iabmc - maxchi)
      b2 = dble(iaam - maxchi)
      b3 = dble(ibbm - maxchi)
      do 1 ichi = 1, maxchi
      xchi = dble(ichi)
      xb = (a1 - xchi) * (a2 - xchi) * (a3 - xchi)
1     x = one - (b1 + xchi) * (b2 + xchi) * (b3 + xchi) * x / xb
      if(x) 4,3,2
4     x = -x
      ll = 1
2     x = dlog(x) - si(iabmc+1) - si(iaam+1) - si(ibbm+1) - si(ibcam+1)
     :            - si(iacbm+1) - si(minchi+1)
      x = two*x + si(iapam) + si(iamam) + si(ibpbm) + si(ibmbm) +
     :           si(icpcm) + si(icmcm) + delta
      x = dexp( x / two)
      l = ll + minchi + nint(b - a + cm)
      if( 2 * (l/2) .ne. l) x = -x
3     xf3j = x
      return
      end
      function xf6j(a,b,e,d,c,f)
*
*                                    | a  b  e |
*   program to compute the 6j symbol {         }
*                                    | d  c  f |
*   author: b. follmeg
*   current revision date: 4-may-1997
*
* -------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      data tol,zero,one,two /1.d-10,0.d0,1.d0,2.d0/
      x=zero
* check triangular conditions for triad ( a b e)
      if ((e .gt. (a + b)) .or. (e .lt. abs(a - b))) goto 40
      sum = a + b + e
      if (mod(sum,one) .gt. tol) goto 40
      iabe = nint(sum)
* check triangular conditions for triad ( d c e)
      if ((e .gt. (c + d)) .or. (e .lt. abs(c - d))) goto 40
      sum = d + c + e
      if (mod(sum,one) .gt. tol) goto 40
      idce = nint(sum)
* check triangular conditions for triad ( a c f)
      if ((f .gt. (a + c)) .or. (f .lt. abs(a - c))) goto 40
      sum = a + c + f
      if (mod(sum,one) .gt. tol) goto 40
      iacf = nint(sum)
* check triangular conditions for triad ( d b f)
      if ((f .gt. (d + b)) .or. (f .lt. abs(d - b))) goto 40
      sum = d + b + f
      if (mod(sum,one) .gt. tol) goto 40
      idbf = nint(sum)
      iabdc = nint(a + b + d + c)
      iaedf = nint(a + e + d + f)
      ibecf = nint(b + e + c + f)
      minchi = max(iabe,idce,iacf,idbf,-1)
      maxchi = min(iabdc,iaedf,ibecf) - minchi
* indices for deltas
      delta = zero
      i2a = nint(two*a) - 1
      i2b = nint(two*b) - 1
      i2c = nint(two*c) - 1
      i2d = nint(two*d) - 1
      i2e = nint(two*e) - 1
      i2f = nint(two*f) - 1
* delta(abe)
      j1 = iabe - i2a
      j2 = iabe - i2b
      j3 = iabe - i2e
      j4 = iabe + 2
      delta = delta + si(j1) + si(j2) + si(j3) - si(j4)
* delta(dce)
      j1 = idce - i2d
      j2 = idce - i2c
      j3 = idce - i2e
      j4 = idce + 2
      delta = delta + si(j1) + si(j2) + si(j3) - si(j4)
* delta(acf)
      j1 = iacf - i2a
      j2 = iacf - i2c
      j3 = iacf - i2f
      j4 = iacf + 2
      delta = delta + si(j1) + si(j2) + si(j3) - si(j4)
* delta(dbf)
      j1 = idbf - i2d
      j2 = idbf - i2b
      j3 = idbf - i2f
      j4 = idbf + 2
      delta = delta + si(j1) + si(j2) + si(j3) - si(j4)
      delta = 0.5d0 * delta
      iabdc = iabdc - minchi
      iaedf = iaedf - minchi
      ibecf = ibecf - minchi
      iabe = minchi - iabe
      idce = minchi - idce
      iacf = minchi - iacf
      idbf = minchi - idbf
      abdc = dble(iabdc - maxchi)
      aedf = dble(iaedf - maxchi)
      becf = dble(ibecf - maxchi)
      abe = dble(maxchi + iabe + 1)
      dce = dble(maxchi + idce + 1)
      acf = dble(maxchi + iacf + 1)
      dbf = dble(maxchi + idbf + 1)
* loop over chi
      x = 1.d0
      ipower = 0
      if (maxchi .le. 0) goto 30
      ii = minchi + maxchi + 2
      do 10 ichi = 1, maxchi
         xchi = dble(ichi)
         xa = (abdc + xchi) * (aedf + xchi) * (becf + xchi)
         xb = (abe - xchi) * (dce - xchi) * (acf - xchi) * (dbf - xchi)
         x = 1.d0 - xa * (ii - ichi) * x / xb
10    continue
      if (x) 20,40,30
20    x = -x
      ipower = 1
30    x = dlog(x) + si(minchi+2)-si(iabdc+1)-si(iaedf+1)-si(ibecf+1)
     :            - si(iabe+1)-si(idce+1)-si(iacf+1)-si(idbf+1) + delta
      ipower = ipower + minchi
      x = (-1.d0)**ipower * dexp(x)
40    xf6j = x
      return
      end
* -------------------------------------------------------------------
      double precision function xf9j(a,b,c,d,e,f,g,h,i)
*
*                                        | a b c |
*... function to compute the 9j - symbol { d e f }
*                                        | g h i |
*
*    this is done by contracting three 6j - symbols.
*    (see edmonds, p. 101, eq. (6.4.3))
*     current revision date: 3-apr-91 by mha
*
      implicit double precision (a-h,o-z)
      double precision i
      logical flaghf
      common /cofact/ si(1)
      dimension xj1(6),xj2(6)
      data tol,one /1.d-5,1.d0/
      x=0.d0
* correct for possible round-off error in angular momenta
      ia=2*(a+tol)
      a=ia/2.d0
      ib=2*(b+tol)
      b=ib/2.d0
      ic=2*(c+tol)
      c=ic/2.d0
      id=2*(d+tol)
      d=id/2.d0
      ie=2*(e+tol)
      e=ie/2.d0
      if=2*(f+tol)
      f=if/2.d0
      ig=2*(g+tol)
      g=ig/2.d0
      ih=2*(h+tol)
      h=ih/2d0
      ii=2*(i+tol)
      i=ii/2d0
* check triangular conditions for triad ( a b c)
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 150
      sum = a + b + c
      if (mod(sum,one) .gt. tol) goto 150
* check triangular conditions for triad ( d e f)
      if ((f .gt. (d + e)) .or. (f .lt. abs(d - e))) goto 150
      sum = d + e + f
      if (mod(sum,one) .gt. tol) goto 150
* check triangular conditions for triad ( g h i)
      if ((i .gt. (g + h)) .or. (i .lt. abs(g-h))) goto 150
      sum = g + h + i
      if (mod(sum,one) .gt. tol) goto 150
* check triangular conditions for triad ( a d g)
      if ((g .gt. (a + d)) .or. (g .lt. abs(a - d))) goto 150
      sum = d + a + g
      if (mod(sum,one) .gt. tol) goto 150
* check triangular conditions for triad ( b e h)
      if ((h .gt. (b + e)) .or. (h .lt. abs(b - e))) goto 150
      sum = b + e + h
      if (mod(sum,one) .gt. tol) goto 150
* check triangular conditions for triad ( c f i)
      if ((i .gt. (c + f)) .or. (i .lt. abs(c - f))) goto 150
      sum = c + f + i
      if (mod(sum,one) .gt. tol) goto 150
* restrict sum
* j1 and j2 are pairs that are coupled with kappa in triangular relations
      xj1(1)=a
      xj1(2)=h
      xj1(3)=b
      xj1(4)=d
      xj1(5)=a
      xj1(6)=f
      xj2(1)=i
      xj2(2)=d
      xj2(3)=f
      xj2(4)=h
      xj2(5)=i
      xj2(6)=b
      xmx=1.e+20
      xmn=-1.d0
      do 100 is=1,6
        sum=xj1(is)+xj2(is)
        dif=abs(xj1(is)-xj2(is))
        if (sum .lt. xmx) xmx=sum
        if (dif .gt. xmn) xmn=dif
100   continue
      if (xmn .gt. xmx) then
        print *,' kmn = ',i3,' .gt. kmx = ',i3,' in xf9j; abort'
        call exit
      endif
* main sum
      k2=2*xmn+tol
      k1=xmn+tol
      if (k2 .ne. 2*k1) then
        flaghf=.true.
        kmn=xmn-0.5d0+tol
        kmx=xmx-0.5d0+tol
      else
        flaghf=.false.
        kmn=xmn+tol
        kmx=xmx+tol
      endif

      do 120 kap=kmn,kmx
        if (flaghf) then
          xkap=kap+0.5d0
          iphase=(-1)**(2*kap+1)
          xnorm=2*kap+2
        else
          xkap=kap
          iphase=(-1)**(2*kap)
          xnorm=2*kap+1
        endif
        t1=xf6j(a,d,g,h,i,xkap)
        t2=xf6j(b,e,h,d,xkap,f)
        t3=xf6j(c,f,i,xkap,a,b)
        x=x+iphase*xnorm*t1*t2*t3
120   continue
150   xf9j=x
      return
      end
* -------------------------------------------------------------------
      function f3j0(j1,j2,j3)
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      j=j1+j2+j3
      jp=j/2
      jpp=jp+1
      jp1=jpp-j1
      jp2=jpp-j2
      jp3=jpp-j3
      x=si(jpp)-si(jp1)-si(jp2)-si(jp3)+del(j1,j2,j3)
      f3j0=(-1)**jp*dexp(x)
      return
      end
* -------------------------------------------------------------------
      function f6j(j1,j2,j3,l1,l2,l3)
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      f6j=0.d0
      if(j3.gt.(j1+j2)) return
      if(j3.lt.iabs(j2-j1)) return
      if(l3.gt.(j1+l2)) return
      if(l3.lt.iabs(l2-j1)) return
      if(l3.gt.(l1+j2)) return
      if(l3.lt.iabs(j2-l1)) return
      if(j3.gt.(l1+l2)) return
      if(j3.lt.iabs(l2-l1)) return
      ia1=j1+j2+l1+l2
      ia2=j2+j3+l2+l3
      ia3=j3+j1+l3+l1
      ib1=j1+j2+j3
      ib2=j1+l2+l3
      ib3=l1+j2+l3
      ib4=l1+l2+j3
      k1=max0(ib1,ib2,ib3,ib4,-1)
      k2=min0(ia1,ia2,ia3)-k1
      ia1=ia1-k1
      ia2=ia2-k1
      ia3=ia3-k1
      ib1=k1-ib1
      ib2=k1-ib2
      ib3=k1-ib3
      ib4=k1-ib4
      ll=0
      x=1.0d0
      if(k2.le.0) go to 2
      do 1 i=1,k2
      xl=float(k2-i+1)
      xa=(ia1-xl+1.0d0)*(ia2-xl+1.0d0)*(ia3-xl+1.0d0)
      xb=(ib1+xl)*(ib2+xl)*(ib3+xl)*(ib4+xl)
1     x=1.0d0-xa*(k1+xl+1.0d0)*x/xb
      if(x) 4,3,2
4     x=-x
      ll=1
2     x=dlog(x)+si(k1+2)-si(ia1+1)-si(ia2+1)-si(ia3+1)-si(ib1+1)
     #-si(ib2+1)-si(ib3+1)-si(ib4+1)
      x=x+del(j1,j2,j3)+del(j1,l2,l3)+del(l1,j2,l3)+del(l1,l2,j3)
      ll=ll+k1
      x=(-1)**ll*dexp(x)
3     f6j=x
      return
      end
* -------------------------------------------------------------------
      function del(j1,j2,j3)
      implicit double precision (a-h,o-z)
      common /cofact/ si(1)
      j=j1+j2+j3+1
      k1=j-2*j1
      k2=j-2*j2
      k3=j-2*j3
      k=k1
      l1=k2
      l2=k3
      if(k.ge.k2) go to 7
      k=k2
      l1=k1
      l2=k3
7     if(k.ge.k3) go to 4
      k=k3
      l1=k1
      l2=k2
4     d=0.
      do 3 l=k,j
3     d=d+dlog(dfloat(l))
      del=(si(l1)+si(l2)-d)/2.0d0
      return
      end
* -------------------------------------------------------------------
      integer function lenstr(string)
      character*(*) string
      l=len(string)
      do 10 i=l,1,-1
10    if(string(i:i).ne.' ') goto 20
      lenstr=0
      return
20    lenstr=i
      return
      end
cstart unix-convex cray mac
c;      subroutine exit
c;      write (6, 100)
c;100   format (' *** EXIT FROM HIBRIDON ***')
c;      stop
c;      end
cend
* ----------------------------------------------------------------------------
      subroutine intairy(x, xairy, xbiry)
* ----------------------------------------------------------------------------
* subroutine to evaluate integrals of airy functions using chebyshev
* polynomial expansions from G. Nemeth, Magyar Tud. Akad. Mat. Fiz. Oszt.
* Kozl. 20, 13 (1971).

*  returns xairy and xbiry, where
* ----------------------------------------------------------------------------
*  for x .le. 0
* ----------------------------------------------------------------------------
*      intai = Int[Ai(-x),0,x] - 2/3
*            = -xairy*cos(zeta-pi/4) + xbiry*sin(zeta+pi/4)
*      xbiry = Int[Bi(-x),0,x]
*            = xairy*sin(zeta-pi/4) + xbiry*cos(zeta+pi/4)
*      where zeta = (2/3)*(-x)^1.5
* ----------------------------------------------------------------------------
*  for x > 0
* ----------------------------------------------------------------------------
*      xairy = exp(zeta)*{1/3 - Int[Ai(x),0,x]}
*      xbiry = exp(-zeta)*Int[Bi(x),0,x]
*      where xi = (2/3)*x^1.5
* ----------------------------------------------------------------------------
*  written by:  millard alexander
*  current revision date:  11-feb-1991 (modified for sun by
*  alexandra borysow)
* ----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a3(22),r1(17),b3(29),s1(36)
      dimension c3(25),d3(26),p3(21),q3(25)
      dimension cn(20)
      data x2 / 5.24148278841779d0/
      data zero, one, two /0.d0, 1.d0, 2.d0/
      data oneth / 0.33333333333333333d0/
      data twoth / .66666666666666666666d0/
      data pi /3.1415926535897932d0/
c this is expansion of luke for Int[K1/3,{x,Inf}]
      data cn / 9.52775168162172d-1,  -4.445349990836d-2,
     :   2.533588840117d-3,  -2.1194358768d-4,
     :   2.2493477764d-5,  -2.82531583d-6,
     :   4.03411595d-7,   -6.3794575d-8,
     :   1.097293d-8,   -2.026008d-9,
     :   3.97563d-10,   -8.227d-11,
     :   1.7842d-11,   -4.035d-12,
     :   9.47d-13,   -2.3d-13,    5.8d-14,
     :  -1.5d-14,    4.d-15,   -1.d-15/
c here follow nemeth's coefficients
      data a3 /1.6730995554933d-01,
     :  -1.38833749391644d-01,    4.269397609363d-02,
     :  -7.77632245783d-03,   -3.92114983777d-04,
     :   8.41351286511d-04,   -3.04480235087d-04,
     :   5.3284908589d-05,    9.30472131d-07,
     :  -3.321295855d-06,    9.42808894d-07,
     :  -1.08124674d-07,   -1.2158869d-08,
     :   7.397362d-09,   -1.327282d-09,
     :   4.9601d-11,    3.1827d-11,
     :  -8.216d-12,    7.68d-13,
     :   6.9d-14,   -3.2d-14,    4.d-15/
      data  r1 /9.68611686973613d-01,
     :  -3.012965177208d-02,    1.182866477802d-03,
     :  -6.9897206095d-05,    5.346267929d-06,
     :  -4.91980038d-07,    5.2176903d-08,
     :  -6.200693d-09,    8.09646d-10,
     :  -1.14492d-10,    1.7342d-11,
     :  -2.79d-12,    4.73d-13,
     :  -8.4d-14,    1.6d-14,
     :  -3.d-15,    1.d-15/
      data b3 /1.9746423175415449d+01,
     :   3.3421598170057534d+01,    2.3112132392259135d+01,
     :   1.3441327140061466d+01,  6.749710207607081d+00,
     :   3.007841304580219d+00,  1.205402544150081d+00,
     :   4.40612841135375d-01,  1.48452873274722d-01,
     :   4.6447156901552d-02,   1.359618591953d-02,
     :   3.743995685863d-03,    9.74415692922d-04,
     :   2.4074747858d-04,  5.665796275d-05,
     :   1.2741213076d-05,  2.745709407d-06,
     :   5.68338025d-07,  1.13250222d-07,
     :   2.176799d-08,  4.042866d-09,
     :   7.26726d-10,   1.2662d-10,
     :   2.1412d-11,  3.519d-12,
     :   5.63d-13,  8.8d-14,  1.3d-14,  2.d-15/
      data s1 /1.043373261895484d0,
     :   4.6562348283283d-02,  3.621330668957d-03,
     :   4.83205510468d-04,  4.1127965419d-05,
     :  -1.8136189822d-05,  -9.703694198d-06,
     :  -7.34800816d-07,   9.61099744d-07,
     :   2.3505703d-07, -1.08774363d-07,
     :  -4.0117513d-08,   1.6681504d-08,
     :   6.11259d-09, -3.30861d-09,
     :  -7.72088d-10,  7.27928d-10,
     :   3.2063d-11,  -1.53606d-10,
     :   2.6875d-11,  2.6486d-11, -1.3029d-11,
     :  -2.217d-12,  3.714d-12, -7.1d-13,
     :  -6.49d-13,  4.16d-13,  1.d-15, -1.07d-13,
     :   4.5d-14,  9.d-15, -1.6d-14,  5.d-15,  2.d-15,
     :  -2.d-15,   1.d-15/
      data c3 /
     :   2.90568567087045d-1,   -1.68969223877061d-1,
     :  -4.271408503181d-2,    8.4978293684633d-2,
     :   9.065203040869d-3,   -2.0322262881791d-2,
     :  -5.068718344054d-3,    2.113265160901d-3,
     :   1.036734509965d-3,   -1.7469500078d-5,
     :   -9.5657221195d-5,   -1.65848004d-5,
     :   3.52834694d-6,    1.642011754d-6,
     :   8.8939642d-8,   -6.7283625d-8,
     :  -1.5048855d-8,   4.4065d-10,
     :   6.52841d-10,    8.2431d-11,
     :  -9.827d-12,   -3.998d-12,
     :  -2.7d-13,    7.4d-14,    1.7d-14/
      data d3 /
     :   1.87484069565758d-1,   -2.94866612241399d-1,
     :   1.66193041700529d-1,    3.519023817395d-3,
     :  -4.3311223754907d-2,   -8.318876773784d-3,
     :   7.527677093084d-3,    2.526578667711d-3,
     :  -4.0984074063d-4,   -3.47588883407d-4,
     :  -2.9121350908d-5,    2.1320634725d-5,
     :   5.927200433d-6,   -2.58263507d-7,
     :  -3.70409669d-7,   -4.9404676d-8,
     :   8.914115d-9,    3.488812d-9,
     :   1.93554d-10,   -9.6744d-11,
     :  -2.0696d-11,    7.2d-14,
     :   6.2d-13,    8.3d-14,
     :  -5.d-15,   -3.d-15/
      data p3 /
     :   9.93614142516261d-1,   -6.213543340449d-3,
     :   1.61900652955d-4,   -9.418135576d-6,
     :   8.66969841d-7,   -1.07894731d-7,
     :   1.6651496d-8,   -3.023183d-9,
     :   6.23691d-10,   -1.42722d-10,
     :   3.5598d-11,   -9.551d-12,
     :   2.729d-12,   -8.24d-13,
     :   2.61d-13,   -8.6d-14,
     :   3.d-14,   -1.1d-14,
     :   4.d-15,   -1.d-15,    1.d-15/
      data q3 /
     :   1.108268927151576d0,   -2.9247192676291d-2,
     :   1.261229956734d-3,   -9.8492949041d-5,
     :   1.1090415008d-5,   -1.604311887d-6,
     :   2.78860541d-7,   -5.5829252d-8,
     :   1.2511344d-8,   -3.07548d-9,
     :   8.16959d-10,   -2.31862d-10,
     :   6.9687d-11,   -2.2024d-11,
     :   7.277d-12,   -2.502d-12,
     :   8.92d-13,   -3.28d-13,
     :   1.24d-13,   -4.8d-14,
     :   1.9d-14,   -8.d-15,
     :   3.d-15,   -1.d-15,    1.d-15/
* here for nemeth calculation
*      if (x .gt. 3.83154716d0) then
*        chi=twoth*x*sqrt(x)
*        xx=5/chi
*        call cheby(fairy,19,xx,cn,2)
*        xintk=sqrt(pi/(two*chi))*exp(-chi)*fairy
*        xairy=oneth-(one/(pi*sqrt(3.d0)))*xintk
*        write (6, 10) x, chi,  xairy, oneth-xairy
*      endif
      if (x .eq. zero) then
        xairy=-twoth
        xbiry=zero
      else if (x .gt. zero) then
        chi=twoth*x*sqrt(x)
        if (x .le. x2) then
          xx=x/x2
          exx=exp(chi)
          call cheby(fairy, 21, xx, a3, 2)
          call cheby(fbiry, 28, xx, b3, 2)
          xairy=x*fairy
          xbiry=x*fbiry
          xairy=oneth-xairy
          xairy=xairy*exx
          xbiry=xbiry/exx
        else
          xx=8/chi
          call cheby(fairy,16,xx,r1,2)
          xairy=fairy/sqrt(6.d0*pi*chi)
          call cheby(fbiry,35,xx,s1,2)
          xbiry=fbiry/sqrt(pi*1.5d0*chi)
        endif
      else if (x .lt. zero) then
        x32=-x*sqrt(abs(x))
        chi=twoth*x32
        if (x .gt. -x2) then
          xx=-x/x2
          call cheby(fairy, 24, xx, c3, 2)
          call cheby(fbiry, 25, xx, d3, 2)
          xairy=-x*fairy
          xbiry=-x*fbiry
          xairy=-twoth+xairy
          arg=chi-0.25d0*pi
          cs=cos(arg)
          sn=sin(arg)
          qq3=-xairy*cs+xbiry*sn
          pp3=xairy*sn+xbiry*cs
        else if (x .le. -x2) then
          xx=8/chi
          call cheby(pp3,20,xx,p3,1)
          call cheby(qq3,24,xx,q3,1)
          qq3=qq3/(two*chi)
          xmod=one/sqrt(pi*x32)
          qq3=qq3*xmod
          pp3=pp3*xmod
        endif
        xairy=qq3
        xbiry=pp3
      endif
      end
* ------------------------------------
      subroutine cheby(function, n_max, x, ac, iexp)
*  subroutine to evaluate a function which is expanded in a finite series of
*  Chebyshev functions (if iexp = 0)
*             n_max
*     f(x) =  Sigma  [ ac   T (x) ]
*              n=0       n   n
*  as shown in Y. L. Luke, "Mathematical functions and their approximations,"
*  (Academic, 1975) p. 478, this can be rewritten as
*     f(x) = b   - x  b
*             0        1
*  where the b coefficients satisfy the downward recursion relation:
*     b  = 2x b    - b     + ac ,  0 .le. n .le n_max
*      n       n+1    n+2      n

*    with b        = b        = 0
*          n_max+1    n_max+2
*  or (if iexp = 1) even Chebyshev functions

*             n_max
*     f(x) =  Sigma  [ ac   T  (x) ]
*              n=0       n   2n
*  as shown in Y. L. Luke, "Mathematical functions and their approximations,"
*  (Academic, 1975) p. 478, this can be rewritten as
*                     2
*     f(x) = b   - (2x  - 1) b
*             0               1
*  where the b coefficients satisfy the downward recursion relation:
*               2
*     b  = 2 (2x  - 1)b    - b     + ac ,  0 .le. n .le n_max
*      n               n+1    n+2      n

*    with b        = b        = 0
*          n_max+1    n_max+2
*  or (if iexp = 2)
*             n_max          *
*     f(x) =  Sigma  [ ac   T (x) ]
*              n=0       n   n
*  as shown in Y. L. Luke, "Mathematical functions and their approximations,"
*  (Academic, 1975) p. 478, this can be rewritten as

*     f(x) = b   - (2x  - 1) b
*             0               1
*  where the b coefficients satisfy the downward recursion relation:

*     b  = 2 (2x  - 1)b    - b     + ac ,  0 .le. n .le n_max
*      n               n+1    n+2      n

*    with b        = b        = 0
*          n_max+1    n_max+2
*
* ----------------------------------------------------------------------------
*  written by:  millard alexander
*  current revision date:  19-nov-1991
* ----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension ac(0:n_max)
      data zero,half, one, two /0.d0,0.5d0,1.d0,2.d0/
* carry out downward recursion for b-coefficients
      if (iexp .eq. 0) then
        txsq = two * x
      elseif (iexp .eq. 1) then
        txsq = two*(two * x * x - one)
      else if (iexp .eq. 2) then
        txsq = two*(two * x - one)
      endif
      b1 = zero
      b0 = zero
*      print *, ' txsq: ', txsq
      do 10  i = n_max, 0, -1
        b2 = b1
        b1 = b0
        b0 = txsq * b1 - b2 + ac(i)
*  at this point b0, b1, and b2 contain, respectively, b(i), b(i+1), and
*  b(i+2)
*      print *, 'i, b0, b1, b2, ac(i): ', i, b0, b1, b2, ac(i)
10    continue
      function = b0 - half* txsq * b1
      return
      end
