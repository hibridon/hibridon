inp=H2ohe_paraq.inp
emax=500
spac=0.25
rendld=100
airyfl=f
job=h2o1q
ener=337.159
run
inp=H2ohe_paraq.inp
emax=500
spac=0.25
rendld=100
airyfl=f
job=h2o2q
ener=370.134
run
prsbr,h2o1q,1,h2o2q,1,1,1,-1,2,0,1
