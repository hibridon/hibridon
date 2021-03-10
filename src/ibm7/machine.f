*fordeck machine.f $Revision: 94.1 $
      program machine
      common /clseg/ iclseg
      character*80 type
      call mctype(ltype,type)
      write (6,1) type(1:ltype)
1     format(a)
      end
