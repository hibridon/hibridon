      program hkey
      character*1 response
      character*80 reward1, sscript1, sscript2, sscript3, sscript4,
     :   sscript5, punish, reward2, reward3, noexist, sscript2a,
     :   reward0
      logical existf
      sscript1='echo "new hibridon user:"> reg.file;exit'
      sscript2=
*     + 'csh set who="$user@`hostname`"; echo $who >> reg.file;exit'
cstart unix-aix unix-sun unix-hp 
     + 'who=`logname`"@"`hostname`; echo $who>>reg.file ;exit'
cend
cstart unix-dec
c;     + 'who=$USER"@"`hostname`; echo $who>>reg.file ;exit'
cend
      sscript2a='grep domain /etc/resolv.conf >>reg.file;exit'
      sscript3=
     +  'echo `uname -a` >>reg.file;exit'
      sscript5=
     + 'mail -s hibridon mha@mha-ibm4.umd.edu <reg.file;rm reg.file'
      punish='cat key |sed 1,1d >key1;mv key1 key;exit'
      reward0='cat key |sed "snxxxn199n" > key1;exit'
      reward1=
     : 'cat key1|tr ''[5-9][0-4]'' ''[0-4][5-9]'' > hkey.uu;exit'
      reward2='uudecode hkey.uu;rm hkey.uu;rm key1;exit'
      reward3='mv key.f ../src/hiamp.f;exit'
      inquire(file='key',exist=existf)
      noexist=
     + 'echo "`hibriddir`/bin/key not present; abort install"'
      if (.not.existf) then
        call system(noexist)
        stop
      endif
10    write (6,13) 
13    format($,'do you accept the license agreement [y/n]? ')
      read (5, 14, end=99) response
14    format(a1)
      if (response .eq. 'y') go to 130
      If (response .eq. 'n') go to 100
      print *, 'response must be either y or n'
      go to 10
100   call system(punish)
      sscript4='echo "license not accepted" >>reg.file;exit'
      goto 150
130   if (response .eq. 'y') 
     +   sscript4='echo "license accepted" >> reg.file;exit'
150   call system(sscript1)
      call system(sscript2)
      call system(sscript2a)
      call system(sscript3)
      call system(sscript4)
      call system(sscript5)
      if (response .eq. 'n') stop
      call system(reward0)
      call system(reward1)
      call system(reward2)
      call system(reward3)
99    end      
