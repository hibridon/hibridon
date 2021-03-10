*fordeck comsplit.f $Revision: 92.1 $
        character*80 line,file,ext  
        write(6,5) 'Input file: '
        read(5,10) file
5       format(1x,a)
        write(6,5) 'Extension for output files: '
        read(5,10) ext  
        le=1
        if(ext.ne.' ')le=lenstr(ext)
        open(1,file=file,status='OLD')
1       read(1,10,end=20,err=20) line
10      format(a)
        iee=lenstr(line)
        if(line(1:8).eq.'*comdeck'.or.line(1:8).eq.'*COMDECK'.or.
     1     line(1:8).eq.'*fordeck'.or.line(1:8).eq.'*FORDECK') then
        close(2)
        file=line(10:)
        ie=lenstr(file)
        file=file(1:ie)//ext(1:le)
        ie=ie+le
        open(2,file=file(1:ie),status='unknown')
        write(6,15) file(1:ie)
15      format(' file ',(a),' created')
        end if
        write(2,25) line(1:iee)
25      format(10(a))
        goto 1
20      stop
        end
        function lenstr(l)
        character*(*) l
        do 10 i=len(l),1,-1
        if(l(i:i).ne.' ') goto 20
10      continue
        lenstr=0
        return
20      lenstr=i
        return
        end
