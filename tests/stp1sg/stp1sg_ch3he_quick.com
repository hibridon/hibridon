inp=Stp_ch3he_quick.inp
job=ch3heq
run
difcrs,ch3heq,20,0,30,-3,0,180,5,1,40
intcrs
