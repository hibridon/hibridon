      program hkey
      character*1 response
      character*80 reward1, sscript1, sscript2, sscript3, sscript4,
     :   sscript5, punish, reward2, reward3
      logical existf
      sscript1='echo "new hibridon user:"> reg.file;exit'
      sscript2=
*     + 'csh set who="$user@`hostname`"; echo $who >> reg.file;exit'
     + 'who=`logname`"@"`hostname`; echo $who>>reg.file ;exit'
      sscript3=
     +  'echo `uname -a` >>reg.file;exit'
      sscript5=
     + 'mail mha@mha-ibm4.umd.edu <reg.file;unalias rm;rm reg.file'
      punish='unalias rm;rm reg.file;exit'
      reward1='unalias mv;mv key.file hkey.file.Z;exit'
      reward2='uncompress hkey.file.Z;exit'
      reward3='unalias mv;mv hkey.file ../src/hiamp.f;exit'
10    write (6,13) 
13    format($,'do you accept license agreement [y/n]? ')
      read (5, *, end=99) response
      if (response .eq. 'y') go to 130
      if (response .eq. 'n') go to 100
      print *, 'response must be either y or n'
      go to 10
100   call system(punish)
      sscript4='echo "license not accepted" >>reg.file;exit'
      goto 150
130   if (response .eq. 'y') 
     +   sscript4='echo "license accepted" >> reg.file;exit'
      call system(sscript1)
      call system(sscript2)
      call system(sscript3)
      call system(sscript4)
      call system(reward1)
      call system(reward2)
      call system(reward3)
150   call system(sscript5)
99    end      
