inp=h2ohe_test.inp
job=h2ohe
jtot2=20
run
partc,,1,-1
printc
exit
