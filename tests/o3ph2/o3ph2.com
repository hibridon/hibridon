inp=O3ph2.inp
show
run
printc
difcrs,,0,0,0,1,0,180,1,1,50
;job=Job2
;difcrs,,0,0,0,2,0,180,1,1,50