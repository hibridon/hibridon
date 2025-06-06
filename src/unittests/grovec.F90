!#define INT(x) int(x, 8) for some reason, the code was once written with integer8. how could it possibly work, since gfortran refuses to compile it, unless using -fdefault-integer-8.

#define INT(x) x

program test_grovec
use mod_grovec, only: dgrovec_type, igrovec_type
implicit none
type(dgrovec_type) :: dg1
type(igrovec_type) :: ig1
dg1 = dgrovec_type(block_size=1024*1024, num_blocks=1024)
call dg1%set_element(INT(1), 3.d0)
call dg1%append(4.d0)
call dg1%set_element(INT(1024), 1024.d0)
call dg1%set_element(INT(1024*1024+1), 1025.d0)
write (6,*) 'dg1(1)=', dg1%get_element(INT(1))
write (6,*) 'dg1(2)=', dg1%get_element(INT(2))
write (6,*) 'dg1(1024)=', dg1%get_element(int(1024))
write (6,*) 'dg1(1024*1024+1)=', dg1%get_element(INT(1024)*1024+1)
ig1 = igrovec_type(block_size=1024, num_blocks=16)
call ig1%append(42)
write (6,*) 'ig1(1)=', ig1%get_element(INT(1))
call ig1%empty()
call ig1%append(666)
write (6,*) 'ig1(1)=', ig1%get_element(int(1))
end program test_grovec
