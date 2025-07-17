inp=arno_test.inp
out=arno_test.out
wrpart=f
wrsmat=f
jmax=6
job=arno_test
show
label=CC Integral Cross Section tests
run
printc
printc,,0,1
printc,,1,2
label=CC Differential and Steric Effect Tests
wrsmat=t
run
difcrs,,0,-1,2,1,0,20,1,1,5,,,1,1.176,.785
difcrs,,0,-1,3,1,0,20,1,1,5
exit
