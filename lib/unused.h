/* use the following macro to explicitely declare a dummy argument or variable as unused, so that the warnings -Wunused-dummy-argument and -Wunused-variable don't trigger */
#define UNUSED_DUMMY(x) if (.false.) print*,shape(x)
#define UNUSED_VARIABLE(x) if (.false.) print*,shape(x)