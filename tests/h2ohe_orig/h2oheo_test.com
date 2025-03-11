inp=h2oheo_test.inp
job=h2oheo
jtot2=20
run
partc,,1,-1
printc
exit
