echo "comparing Arno results"
mkdir testnew_hold
cp Arno_tes1.ics testnew_hold
cp Arno_tes1.xcs testnew_hold
cp Arno_tes1.dcs testnew_hold
cp Boh2_bound.evl testnew_hold
diff Cc.flx testnew_hold
diff Cc.psi testnew_hold
diff Ccdxsec1.dcs testnew_hold
diff Ccdxsec1.xsc testnew_hold
diff Cctest1.ics testnew_hold
diff Cctest1.pcs testnew_hold
diff Cctest1.psc testnew_hold
diff Cctest1.xsc testnew_hold
diff Ch3itest.flx testnew_hold
diff Ch3itest.psi testnew_hold
diff Cstest1.ics testnew_hold
diff Cstest1.pcs testnew_hold
diff Cstest1.psc testnew_hold
diff Cstest1.xsc testnew_hold
diff Vfit_test1.ics testnew_hold
diff Vfit_test1.xsc testnew_hold
diff Vfit_test1.tcs testnew_hold
diff Vfit_test1.mcs testnew_hold
diff Aroh_new1.xsc testnew_hold
diff Aroh_new1.psc testnew_hold
diff Aroh_new1.ics testnew_hold
diff H2ohe1.ics testnew_hold
diff H2ohe1.xsc testnew_hold
diff Heco1.ics testnew_hold
diff Heco1.xsc testnew_hold
diff Heco1.pcs testnew_hold
diff Heco1.psc testnew_hold
diff Hecn1.xsc testnew_hold
diff Hecn1.hfx testnew_hold
diff Hecn1.xms testnew_hold

exit
