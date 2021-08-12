inp=Hco_testq.inp
job=Hcotrnq
jtot2=4
ener=300.
run
printc
trnprt
exit
