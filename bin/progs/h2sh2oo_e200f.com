inp=h2so_oh2f.inp
jtot2=57
jmax=6
ener=200
job=oof
run
show
intcrs
quit

