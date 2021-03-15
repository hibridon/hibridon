inp=Oxe_1d3p_dcs2.inp
job=E3400u
show
run
printc
difcrs,,0,1,0,3,0,180,.02,1,300
quit
