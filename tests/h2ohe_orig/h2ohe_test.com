inp=h2ohe_test.inp
job=h2ohe
jtot2=20
run
printc
partc,,1,-1
exit
