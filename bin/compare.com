echo "comparing Arno results"
mkdir testnew_hold
cp Arno_tes1.ics testnew_hold
cp Arno_tes1.xcs testnew_hold
cp Arno_tes1.dcs testnew_hold
cp Boh2_bound.evl testnew_hold
cp Cc.flx testnew_hold
cp Cc.psi testnew_hold
cp Ccdxsec1.dcs testnew_hold
cp Ccdxsec1.xsc testnew_hold
cp Cctest1.ics testnew_hold
cp Cctest1.pcs testnew_hold
cp Cctest1.psc testnew_hold
cp Cctest1.xsc testnew_hold
cp Ch3itest.flx testnew_hold
cp Ch3itest.psi testnew_hold
cp Cstest1.ics testnew_hold
cp Cstest1.pcs testnew_hold
cp Cstest1.psc testnew_hold
cp Cstest1.xsc testnew_hold
cp Vfit_test1.ics testnew_hold
cp Vfit_test1.xsc testnew_hold
cp Vfit_test1.tcs testnew_hold
cp Vfit_test1.mcs testnew_hold
cp Aroh_new1.xsc testnew_hold
cp Aroh_new1.psc testnew_hold
cp Aroh_new1.ics testnew_hold
cp H2ohe1.ics testnew_hold
cp H2ohe1.xsc testnew_hold
cp Heco1.ics testnew_hold
cp Heco1.xsc testnew_hold
cp Heco1.pcs testnew_hold
cp Heco1.psc testnew_hold
cp Hecn1.xsc testnew_hold
cp Hecn1.hfx testnew_hold
cp Hecn1.xms testnew_hold

exit
