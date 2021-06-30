inp=heco_test.inp
job=heco
jmax=5
wrpart=t
run
partc,,2,0
printc
exit
