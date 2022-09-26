! https://stackoverflow.com/questions/37500015/a-portable-way-to-suppress-an-unused-dummy-argument-warning-in-fortran
#define UNUSED(x) if (.false.) print*,shape(x)
