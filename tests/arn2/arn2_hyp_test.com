inp=Arn2_hyp.inp
out=arn2_hyp.out
job=arn2_hyp
prxsec=t
show
label=Hyperfine structure with 2 spins
run
hypxsc,arn2_hyp,1,202,0,5
exit

