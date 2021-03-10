C THIS PROGRAM CAN BE CONVERTED TO THE APPROPRIATE MACHINE
c TYPE BY ITSELF. FIRST COMPILE AND LINK, AND THEN RUN AND
c USE THIS CODE AS INPUT.
      program ftconv
      parameter (mxtype=10)
      character*30 mx(mxtype)
      integer lx(mxtype)
      character*30 inclh,inclt,comfil
      logical resolv,lowerc,stripb
      character*132 l
      character*4 ext
      common/cext/ ext
      print*,'*** fortran conversion program ***'
      ext='.new'
      stripb=.false.
      mx(1)=' '
      mx(2)=' '
      nx=1
      nxx=1
cstart unknown
      nx = 0
1     write(6,10) 'Machine type:  '
10    format(1x,a,$)
      nx = nx+1
      if (nx.gt.mxtype) stop 'too many types'
      read(5,20) mx(nx)
20    format(a)
      call lower(mx(nx))
      if (mx(nx).ne.' ') goto 1
      nx = nx-1
cend
cstart molpro
c;      mx(2) ='molpro'
c;      nx=2
c;      nxx=1
cend
cstart disco-u
c;      mx(1) ='disco'
c;      nx=2
c;      nxx=2
cend
cstart disco-f
c;      mx(2) ='disco'
c;      nx=2
c;      nxx=1
cend
cstart dec-risc
c;      mx(nxx)='dec-risc'
c;      ext='.dec'
cend
cstart unix-hp
c;      mx(nxx)='unix-hp'
c;      ext='.hpu'
cend
cstart unix-mac
c;      mx(nxx)='unix-mac'
c;      ext='.mac'
cend
cstart unix-convex
c;      mx(nxx)='unix-convex'
c;      ext='.cvx'
cend
cstart unix-sequent
c;      mx(nxx)='unix-sequent
c;      ext='.seq'
cend
cstart ibm-risc
c;        mx(nxx)='ibm-risc'
c;        ext='.ibm'
cend
cstart ibm-vm
c;      mx(nxx)='ibm-vm'
c;      ext='.ibm'
cend
cstart vax
c;      mx(nxx)='vax'
c;      ext='.vax'
cend
cstart cray-cos
c;      mx(nxx)='cray-cos'
c;      ext='.cos'
cend
cstart cray-unicos
c;      mx(nxx)='cray-unicos'
c;      ext='.ucs'
cend
cstart fps
c;      mx(nxx)='fps'
c;      ext='.fps'
cend
      if(mx(nxx).eq.' ') then
        write(6,25)
25      format(' NO MACHINE GIVEN')
        stop
      end if
      do 40 i=1,nx
      if(mx(i).eq.'cray') then
         write(6,30)
30       format(' machine=cray is insufficient machine type'/
     1          ' specify cray-cos or cray-unicos')
cstart unknown
c;         nx = 0
c;         goto 1
cend
         stop
      end if
      if(mx(i).eq.'unix') then
         write(6,31)
31       format(' machine=unix is insufficient machine type'/
     1          ' specify unix-hp, unix-mac, unix-convex,',
     2          ' unix-sequent, ibm-risc, or dec-risc')
cstart unknown
c;         nx = 0
c;         goto 1
cend
         stop
      else if(mx(i).eq.'ibm') then
         write(6,32)
32       format(' machine=ibm is insufficient machine type'/
     1          ' specify ibm-risc, or ibm-vm')
cstart unknown
c;         nx = 0
c;         goto 1
cend
         stop
      end if
40    continue
      nxx=nx
      do 22 i=1,nxx
      lx(i) = lenstr(mx(i))
      do 22 k=3,4
      if (mx(i)(k+1:k+1).eq.'-') then
        nx = nx+1
        mx(nx) = mx(i)(1:k)
        lx(nx) = k
      end if
22    continue
      do 24 k=1,nx
24    write (6,26) mx(k)
26    format(1x,'machine=',a)
      l=mx(1)
      call lower(l)
c
      if(index(l,'unix-mac').ne.0) then
c
      inclh='      include '':common:'
      inclt=''''
      resolv=.false.
      lowerc=.true.
      comfil='common.mac'
c
      else if(index(l,'unknown').ne.0) then
c
      inclh='      include "common/'
      inclt='"'
      resolv=.false.
      lowerc=.true.
      comfil='common.all'
c
      else if(index(l,'unix').ne.0 .or. index(l,'dec') .ne. 0
     :        .or. index(l,'ibm-risc') .ne. 0) then
c
      inclh='      include "common/'
      inclt='"'
      resolv=.false.
      lowerc=.true.
      comfil='common.unix'
c
      else if(index(l,'cray').ne.0) then
c
      inclh='      include ''('
      inclt=')'''
      resolv=.true.
      comfil='common.cray'
      lowerc=.false.
c
      else if(index(l,'vax').ne.0) then
c
      inclh='      include ''('
      inclt=')'''
      resolv=.false.
      comfil='common.vax'
      lowerc=.true.
c
      else if(index(l,'ibm').ne.0 .and. index(l,'ibm-risc') .eq.0)
     :     then
c
      inclh='      include ''('
      inclt=')'''
      resolv=.true.
      comfil='common.ibm'
      lowerc=.false.
c
      else if(index(l,'disco').ne.0) then
c
      inclh='*CALL '
      inclt=' '
      resolv=.false.
      comfil='common.cdk'
      lowerc=.false.
c
      else
c
      print*,'unknown machine type: ',l(1:30)
      stop
      end if
      i1=lenstr(inclh)
      i2=lenstr(inclt)
      i3=lenstr(comfil)
      call ftcnv(mx,lx,nx,inclh(1:i1),inclt(1:i2),
     1   comfil(1:i3),resolv,lowerc,stripb)
      stop
      end
      subroutine ftcnv(mx,lx,nx,inclh,inclt,comfil,resol,lowerc
     >  ,stripb)
      character*(*) inclh,inclt,comfil,mx(*)
      integer lx(*)
      logical resolv,lowerc,resol,stripb
      character*132 l,str
      character*4 ext
      logical exist
      logical force
      common/cr/ ieof,nr,nw,icond
      common/cl/ force
      common/cext/ ext
c
      nr = 0
      nw = 0
      resolv=resol
      inp=5
20    format(a)
cstart .not. batch
1     if(inp.eq.5) write(6,10) 'file: or *option:  '
10    format(1x,a,$)
cend
2     read(inp,20,end=101) l
      if(l.eq.' ') stop
      if (l.eq.'quit' .or. l.eq.'exit' .or. l.eq.'q' .or.
     :    l.eq.'q' ) stop
      if (l(1:1).eq.'*') then
c...  this escape allows introduction of further options
	if (l(1:6).eq.'*inclh') inclh = l(8:)
	if (l(1:6).eq.'*inclt') inclt = l(8:)
	if (l(1:7).eq.'*lowerc') lowerc = .true.
	if (l(1:7).eq.'*upperc') lowerc = .false.
	if (l(1:6).eq.'*strip') stripb = .true.
	if (l(1:8).eq.'*nostrip') stripb = .false.
	if (l(1:7).eq.'*resolv') resolv = .true.
	if (l(1:9).eq.'*noresolv') resolv = .false.
	goto 1
      end if
      i=index(l,'.list')
      if(i.ne.0) then
        if(inp.eq.15) stop 'Only one list file!'
        inp=15
        open(15,file=l,status='old')
        goto 2
      end if
      open (unit=1,file=l,form='formatted',access='sequential',
     >  status='old')
      rewind 1
      inquire(1,name=l)
cstart batch
c;      write(6,*) 'input file  ',l(1:lenstr(l))
cend
      i=index(l,'.for')
      if(i.ne.0) then
        l(i:)='.f'
      else
        i=index(l,'.')
        l(i:)=ext
      end if
      inquire(file=l,exist=exist)
      if(exist) then
cstart batch
c;        stop 'file exists'
celse
        write(6,30) l(1:lenstr(l))
30      format(/1x,'file ',a,' exists'/
     1  ' overwrite? (y/n)[y] ',$)
        read(5,20,end=101) str
        if(str(1:1).eq.'n') stop
cend
      end if
      write(6,*) 'output file:  ',l(1:lenstr(l))
      open (2,file=l,access='sequential',status='unknown',
     1   form='formatted')
      rewind 2
c
      if(resolv) then
        inquire(file=comfil,exist=exist)
        if(.not.exist) then
          call upper(comfil)
          inquire(file=comfil,exist=exist)
        end if
        if(.not.exist) then
          print*,'common file ',comfil,' missing'
          stop
        end if
        open(3,file=comfil,access='sequential',status='old',
     1       form='formatted')
        rewind 3
        icond=0
        force=.true.
        call comdek(l,lowerc,mx,lx,nx)
      end if
      icond=0
      force=.true.
      lhh = lenstr(inclh)
      if (inclh(1:6).eq.'*CALL ') lhh=6
      ltt = lenstr(inclt)
c
c...  main code .. read input, process, write output
c...  convert to force case except when quoted string is open
100   call linein(1,l,lowerc,mx,lx,nx)
      if (ieof.lt.0) goto 999
      str=l(1:13)
      call lower(str(1:13))
      if (str(1:5).eq.'*deck'.or.str(1:8).eq.'*fordeck')
     >    write (6,*) l(1:lenstr(l))
      if (str(1:5).eq.'*call'.or.str(7:13).eq.'include'.or.
     2    str(1:8).eq.'#include') then
        call includx(l,resolv,inclh(1:lhh),inclt(1:ltt),lowerc)
      else
c...  end of processing, write line
      call linout(l)
      end if
      goto 100
c...  end of file exit here
999   continue
      close (1)
      close (2)
      print*,'ftconv processing completed;  lines read: ',nr,
     1    '    written: ',nw
cstart .not. batch
      goto 1
cend
  101 continue
c...  end of files
      end
      subroutine includx(l,resolv,inclh,inclt,lowerc)
      character*(*) l,inclh,inclt
      character*132 temp
      logical resolv,lowerc
c...  locate name of include block
      ifound=0
      do 1 ipos=1,80
      if (l(ipos:ipos).ne.' '.and.ifound.eq.0) ifound=1
      if (l(ipos:ipos).eq.' '.and.ifound.eq.1) ifound=2
      if (ifound.eq.2.and.l(ipos:ipos).ne.' '.and.l(ipos:ipos).ne.'('
     > .and.l(ipos:ipos).ne.''''.and.l(ipos:ipos).ne.'"') goto 2
1     continue
      stop 'includ error'
2     do 3 jpos=ipos+1,80
      if (l(jpos:jpos).eq.')'. or. l(jpos:jpos).eq.'''' .or.
     > l(jpos:jpos).eq.' '.or.l(jpos:jpos).eq.'"') goto 4
3     continue
      stop 'include error'
4     continue
c     print*,'after 4; ipos,jpos,string ',ipos,jpos,l(ipos:jpos-1)
c...  strip off any path name that might be present
66    jjpos=max(index(l(ipos:jpos-1),'/'),index(l(ipos:jpos),':'))
c     print*,'after 66; ipos,jpos,jjpos ',ipos,jpos,jjpos
      if (jjpos.eq.0) goto 67
      ipos=ipos+jjpos
      goto 66
67    continue
c     print*,'after 67; ipos,jpos,jjpos ',ipos,jpos,jjpos
c     print*,'l(ipos:jpos-1) ',l(ipos:jpos-1)
      if (resolv) then
        call comwr(l(ipos:jpos-1))
      else
      temp = l(ipos:jpos-1)
      call fcase(temp(1:jpos-ipos),lowerc)
        l = inclh//temp(1:jpos-ipos)//inclt
      call linout(l)
      end if
      return
      end
      subroutine comdek(l,lowerc,mx,lx,nx)
cstart unix-hp unix-mac vax ibm eta unix-sequent unix-convex ibm-risc
      parameter (maxd=200,maxl=10)
cend
cstart cray unknown
c;      parameter (maxd=400,maxl=10000)
cend
      character*(*) l,mx(*)
      integer lx(*)
      character*132 buff(maxl),temp
      character*16 name(maxd)
      common/cr/ ieof,nr,nw,icond
      logical lowerc
      integer istart(maxd),iend(maxd)
      save buff,istart,iend,nd,nl,name
      data nd,nl/0,0/
c...  first parse name from line
      nd=0
      nl=0
      nrrr=nr
      call linein(3,l,lowerc,mx,lx,nx)
      if (ieof.lt.0) goto 12
3     do 1 ipos=9,80
      if (l(ipos:ipos).ne.' ') goto 2
1     continue
2     temp = l(ipos:)
      l = temp
c...  store name
      nd = nd+1
      if (nd.gt.maxd) stop 'maxd'
      name(nd) = l
      call upper(name(nd))
      istart(nd) = nl+1
c...  now read lines into memory
5     nl = nl+1
      if (nl.gt.maxl) stop 'maxl'
      call linein(3,buff(nl),lowerc,mx,lx,nx)
      if(ieof.lt.0) goto 12
      if(buff(nl)(1:8).eq.'*comdeck'.or.
     >   buff(nl)(1:8).eq.'*COMDECK') goto 11
10    format((a))
      goto 5
11    l=buff(nl)
      nl=nl-1
      iend(nd)=nl
      goto 3
12    nl = nl-1
      iend(nd) = nl
      nr=nrrr
      return
      entry comwr (l)
c...  write common l to output stream
      call upper (l)
      do 110 id=1,nd
      if (l.eq.name(id)) goto 132
110   continue
      print*,l
      stop 'unknown common'
132   do 130 i=istart(id),iend(id)
130   call linout(buff(i))
      return
      end
      function lenstr(l)
      character*(*) l
      do 10 i=len(l),1,-1
10    if(l(i:i).ne.' ') goto 20
      lenstr=1
      return
20    lenstr=min(132,i)
      return
      end
      subroutine linein(ifil,l,lowerx,mx,lx,nx)
      character*(*) mx(*)
      integer lx(*)
      logical lowerx,lowerc,logic
      character*(*) l
      character*134 lll,str
      logical force
      common/cr/ ieof,nr,nw,icond
      common/cl/ force
c...  force controls whether case forcing is currently on
c     it is off inside a quoted string
c...  icond>0 cstart current, code required
c          =0 cstart not current
c          <0 cstart current, code not required
      ieof=0
      lowerc=lowerx
1     read (ifil,'(a)',end=999) l
      nr = nr+1
      if (l.eq.' ') then
        if(ifil.eq.1) write(2,'()')
        goto 1
      end if
c...  is this a cstart/cend?
      str=l(1:6)
      call lower(str(1:6))
      i=1
      if (str(1:4).eq.'cend') then
        if (icond.eq.0) stop 'illegal cend'
        icond = 0
        i=4
      else if (str(1:5).eq.'celse') then
        if (icond.eq.0) stop 'illegal celse'
        icond = -icond
        i=5
      else if (str(1:6).eq.'cstart') then
        if (icond.ne.0) stop 'illegal cstart'
        icond = -1
        str=l(1:80)
        call lower(str)
        if(logic(str(7:),mx,nx)) icond=1
        i=7
      else if (icond.gt.0.and.str(1:2).eq.'c;') then
        lll = l(3:)
        l = lll
        i=1
      else if (icond.lt.0.and.str(1:2).ne.'c;') then
        lll='c;'//l(1:lenstr(l))
        l=lll
        i=3
      end if
c..  break line into segments delimited by single quote, and forcecase
      if(l(i:i).eq.'c'.or.l(i:i).eq.'C'.or.l(i:i).eq.'*') then
        call fcase(l(1:i),lowerc)
        return
      end if
      i=1
88    j=index(l(i:),'''')
      if(j.eq.0) then
c...    no quote found
        if(force) then
          if(i.gt.80) return
          call fcase(l(i:80),lowerc)
        end if
        return
      else
c...    quote found
        n=1
        k=i+j
89      if(k.gt.80) goto 90
c...    check for more quotes
        if(l(k:k).eq.'''') then
          n=n+1
          k=k+1
          goto 89
        end if
90      if (force) then
        call fcase(l(i:i-1+j),lowerc)
        end if
      end if
      if(mod(n,2).ne.0) force=.not.force
      i=k
      if(i.gt.80) return
      goto 88
999   ieof=-1
      return
      entry linout(l)
      nw = nw+1
      write (2,'(a)') l(1:lenstr(l))
      return
      end
      subroutine fcase(l,lowerc)
      character*(*) l
      logical lowerc
      if(lowerc) then
        call lower(l)
      else
        call upper(l)
      end if
      return
      end
      subroutine upper(l)
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
      subroutine lower(l)
      character*(*) l
      ishift = ichar('a')-ichar('A')
      ia = ichar('A')
      iz = ichar('Z')
      do 1 i=1,len(l)
      j=ichar(l(i:i))
      if (j.ge.ia.and.j.le.iz) l(i:i)=char(j+ishift)
1     continue
      return
      end
      function logic(line,mx,n)
      logical lneg,land,logic,logo
      character*(*) line,mx(n)
      character*5 not,or,and
      data not,or,and/'.not.','.or.','.and.'/
      l=lenstr(line)
      logic=.false.
      land=.false.
      lneg=.false.
      k=1
1     do 10 i=k,l
10    if(line(i:i).ne.' ') goto 20
      goto 50
20    do 21 ie=i+1,l
21    if(line(ie:ie).eq.' '.or.line(ie:ie).eq.'.') goto 22
      ie=l+1
22    ie=ie-1
       if(line(i:i+4).eq.not) then
        lneg=.true.
        k=i+5
        goto 1
      end if
      if(line(i:i+4).eq.and) then
        land=.true.
        k=i+5
        goto 1
      end if
      if(line(i:i+3).eq.or(1:4)) then
        land=.false.
        k=i+4
        goto 1
      end if
      do 30 m=1,n
      ml=lenstr(mx(m))
30    if(line(i:ie).eq.mx(m)(1:ml)) goto 40
      logo=lneg
34    lneg=.false.
      if(land) then
         logic=logic.and.logo
         land=.false.
      else
         logic=logic.or.logo
      end if
      do 35 k=i,l
35    if(line(k:k).eq.' '.or.line(k:k).eq.'.') goto 1
      goto 50
40    logo=.not.lneg
      goto 34
50    return
      end
