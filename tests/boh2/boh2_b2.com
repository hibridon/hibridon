inp=boh2_b2.inp
job=boh2_b2a
show
run
job=boh2_b2b
hsimp=-3
run
exit
