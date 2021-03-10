      implicit double precision (a-h,o-z)
      common /coc/ ncount
      common /col/ nlen
      parameter(kmax=400,kaux=100*kmax)
      dimension det(2)
      dimension aa(kmax,kmax),sc2(kmax)
      dimension a(kmax,kmax),b(1,1),sc(kaux)
      common /cosc11/ sc
      common /cokaux/ naux
c      open (unit=5,file='infile')
      open (unit=7,file='outfile')
      nx=kmax
      naux=max(kaux,1800)
      nop=300000000
      izero=0
      tol=1.e-10
      write (6, 5)
5      format ('iuse, ndim:  ',$)
      read (5, *) iuse, ndim
      if (ndim .gt. 0) then
        nmin=ndim
        nmax=ndim
      else
        nmin=20
        nmax=400
      endif
      do 200 ndim = nmin,nmax,20
      do 10 i=1,ndim
      do 10 j=1,i
      a(i,j)=i*.1+j*j*.2
      a(j,i) = a(i,j)
      aa(j,i)=a(i,j)
10    aa(i,j)=a(i,j)
      cpu1=second()
c.....perform about nop floating point operations between each cpu measurment
      niter=max(1,nop*3/(8*ndim**3))
      flop3=1.d-6*dble(8*niter*(ndim**3)/3)

      do 499 i=1,niter
      do 49 ii=1,ndim
      call dcopy(ndim,aa(1,ii),1,a(1,ii),1)
49    continue
      ncount=0
      call smxinv (a, nx, ndim, sc, sc2, ierr)
*      call dgeicd(a,nx,ndim,izero,rcond,det,sc,naux)
499   continue
      cpu1=second()-cpu1
      write (6,101) ndim,cpu1,niter,flop3/cpu1
      write (7,101) ndim,cpu1,niter,flop3/cpu1
101   format (i4,1pe16.3,i13, 1pe16.3)
      if (ndim .lt. 20) then
        print *, ' a matrix'
        do 110 	i=1,ndim
110       write (6, 115) (aa(i,j), j=1,ndim)
115       format(7f11.5)
        print *, ' ainv matrix'
        do 120 	i=1,ndim
120       write (6, 115) (a(i,j), j=1,ndim)
      endif
200   continue
 
      end
