echo "comparing Arno results"
diff Arno_test1.ics ../tests
diff Arno_test1.xsc ../tests
echo "comparing B-oH2 bound state results"
diff Boh2_bound.evl ../tests
echo "comparing flux and wavefunction tests"
diff Cc.flx ../tests
diff Cc.psi ../tests
echo "comparing differential cross section tests"
diff Ccdxsec1.dcs ../tests
diff Ccdxsec1.xsc ../tests
echo "comparing close-coupling integral cross section tests"
diff Cctest1.ics ../tests
diff Cctest1.pcs ../tests
diff Cctest1.psc ../tests
diff Cctest1.xsc ../tests
echo "comparing restart tests"
diff Ccrstest1.ics ../tests/Cctest1.ics
diff Ccrstest1.xsc ../tests/Cctest1.xsc
diff Ccrstest1.pcs ../tests/Cctest1.pcs
diff Ccrstest1.psc ../tests/Cctest1.psc
echo "comparing CH3I photodissociation tests"
diff Ch3itest.flx ../tests
diff Ch3itest.psi ../tests
echo "comparing coupled-states Ar-N2 results"
diff Cstest1.ics ../tests
diff Cstest1.pcs ../tests
diff Cstest1.psc ../tests
diff Cstest1.xsc ../tests
echo "comparing Ar-N2 'big' results"
diff Ccbtest1.ics ../tests/Cctest1.ics
diff Ccbtest1.xsc ../tests/Cctest1.xsc
diff Ccbtest1.pcs ../tests/Cctest1.pcs
diff Ccbtest1.psc ../tests/Cctest1.psc
echo "comparing Ar-N2 'big' restart tests"
diff Ccbrstest1.ics ../tests/Cctest1.ics
diff Ccbrstest1.xsc ../tests/Cctest1.xsc
diff Ccbrstest1.pcs ../tests/Cctest1.pcs
diff Ccbrstest1.psc ../tests/Cctest1.psc
diff Csbtest1.ics ../tests/Cstest1.ics
diff Csbtest1.xsc ../tests/Cstest1.xsc
diff Csbtest1.pcs ../tests/Cstest1.pcs
diff Csbtest1.psc ../tests/Cstest1.psc
echo "comparing tests of Vfit potential"
diff Vfit_test1.ics ../tests
diff Vfit_test1.tcs ../tests
diff Vfit_test1.xsc ../tests
echo "comparing Ar-OH integral cross section tests"
diff Aroh_new1.xsc ../tests
diff Aroh_new1.psc ../tests
diff Aroh_new1.ics ../tests
exit
