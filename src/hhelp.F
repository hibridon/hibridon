* driver for hibridon help facility
      character*80 line
      common /coipar/ junk(8), lscreen
      lscreen=25
      line ='HELP'
      print *, ' '
      print *, 'HIBRIDON ON-LINE HELP FACILITY'
      print *, ' '
      print *, 'to exit, hit <RETURN> at the "Topic?" prompt'

      call vaxhlp(line)
      end
#if defined(unix-xlf)
@process fixed(132)
#endif
      subroutine vaxhlp(line1)
#if defined(unix-ifort)
cdec$ fixedformlinesize:132
#endif
* latest revision 24-feb-2004
      implicit character*10(z)
      character*(*) line1
      common /coipar/ junk(8),lscreen
#if defined(HIB_CRAY)
c;* return for these machines
c;        write (6,2)
c;2       format (' Help not installed on this machine')
c;        return
c;        end
#endif
#if defined(HIB_UNIX) || defined(HIB_UNIX_DEC)  || defined(HIB_UNIX_IBM)
c simulates vax vms help command
      parameter (inp=5,iout=6,ihlpmain=51,ihlpalt=52)
      parameter (maxlev=9,maxsub=50)
      character*(*) helpdir,helptail
#endif
      include "common/parhlp"
#if defined(HIB_UNIX) || defined(HIB_UNIX_DEC)  || defined(HIB_UNIX_IBM)
      character*80 helpfile,line,lin,lin2
      character*1 pat
#endif
#if defined(HIB_UNIX)
      logical exist
#endif
#if defined(HIB_UNIX) || defined(HIB_UNIX_DEC)  || defined(HIB_UNIX_IBM)
      dimension zreq(maxlev),lreq(maxlev),zsub(maxsub)
c...  by default, looks in file help. if that is not available, then
c     tries file name of first key in the tree.
      helpfile=helpdir//'hibrid'//helptail
      print *, helpfile
#endif
#if defined(HIB_UNIX_CONVEX)
c;      open (ihlpmain,file=helpfile,err=299,readonly,status='old')
#endif
#if defined(HIB_UNIX_SUN) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_AIX) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
c      write(6,*) 'open ',helpfile
      open (ihlpmain,file=helpfile,err=299,status='unknown')
#endif
#if defined(HIB_UNIX) || defined(HIB_UNIX_DEC)  || defined(HIB_UNIX_IBM)
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
        line='HELP'
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
#endif
#if defined(HIB_UNIX)
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
#endif
#if defined(HIB_UNIX_CONVEX)
c;          open (ihlpalt,file='helphelphelphelp',status='old')
c;          open (ihlp,file=helpfile,err=29,readonly,status='old')
#endif
#if defined(HIB_UNIX_SUN) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_AIX) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
          open (ihlp,file=helpfile,err=29,status='unknown')
#endif
#if defined(HIB_UNIX) || defined(HIB_UNIX_DEC)  || defined(HIB_UNIX_IBM)
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
#endif
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
