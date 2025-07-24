# Potential Subroutines

In the usual implementation of [close coupled](./Close-coupled-equations.md) methods the interaction potential $V(R)$ is expanded as in a set of products of functions of the internal degrees of freedom, multiplied by coefficients which depend on the separation coordinate $R$, namely

$$ V(R,q) = \sum V_\lambda(R) \times f_\lambda(q) $$

where $q$ designates the internal coordinates.  Often these are angles, in which case
the functions $f_\lambda$ are spherical harmonics, rotation matrix elements,
or products of these. The particular choice of functions depends on which kind of system you
are describing, and is controlled, within the Hibridon package, by the 
[`BASIS`](basis.md) subroutine.  At present, 18 different types of collision systems are
included.

For each particular choice of collision system, the user also needs to supply the crucial `POT`
subroutine, which delineates the range of values of $\lambda$ in the above equation
and which calculates the expansion coefficients, $V_\lambda(R)$.

## Potentials supplied with Hibridon < 5.0

At present, a number of different `POT` subroutines are supplied with the Hibridon package.
These are listed and described in the following table. Each listed subroutine is found in the Hibridon directory tree as

 - `src/pot_<potname>.f` (fortran subroutines)<p>
 - `bin/progs/potdata/<dataname>` (data files)<p>

where `<potname>` and `<dataname>` are given in a [table](#table-of-potential-subroutines) which
appears below.

To write your own `POT` subroutine, you can use any of the supplied subroutine as an example<p>
If your `POT` requires data from an additional file, make sure you set the flag 
[READPT](./param/readpt.md) to be `.true.`  Data files should be placed in the directory `bin/progs/potdata`.

To [link](./hib4only/linking.md) an executable program, you need to specify the name
of a particular `POT` subroutine.  Specifically, if you wish to link with the potential
subroutine `src/pot/pot_toto.f` you would need to use the [link](./hib4only/linking.md) command

```sh
makehib toto <n>
```

where `<n>` is the maximum number of internal states (channels) used in the expansion 
of the wavefunction. This command will produce the executable program `bin/progs/hib_toto_<n>`

It is also possible to create an executable program which will only test the `POT`
subroutine. To do so, execute the command [`makepot`](./hib4only/makepot.md)

```sh
makepot <potname>
```

This will create the file `bin/progs/runpot_<potname>`

---

## Table of Potential Subroutines


<table ID="1" BORDER=2>

<tr><td><b>potname<p><td><b>dataname<p><td><b>system<p></b>


<tr><td>alhara_mod<td>none<td>AlH(<i>A</i><sup>1</sup>&Pi;)+Ar [modified *ab initio*]

<tr><td><td><td>Yang, M. M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 102, 2413 (1994).

<tr><td>alharx_mod<td>none<td>AlH(<i>X</i><sup>1</sup>&Sigma;<sup>+</sup>)+Ar [modified *ab initio*]

<tr><td><td><td>Yang, M. M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 102, 2413 (1994).

<tr><td>alh2   <td>pot_alh2.dat<td>Al(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Williams-Alexander PES's]

<tr><td><td><td>J. Williams and M. H. Alexander, J. Chem. Phys. 112, 5722(2000); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994); M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995); X. Tan, P. J. Dagdigian, J. Williams, and M. H. Alexander, J. Chem. Phys. 114, 8938 (2001).

<tr><td>alh2_rdep   <td>pot_alh2.dat<td>Al(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Williams-Alexander PES's CBS averaged over v=0]

<tr><td><td><td>J. Williams and M. H. Alexander, J. Chem. Phys. 112, 5722(2000); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994); M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995); X. Tan, P. J. Dagdigian, J. Williams, and M. H. Alexander, J. Chem. Phys. 114, 8938 (2001).

<tr><td>arbha   <td>none<td>BH(<i>A</i><sup>1</sup>&Pi;)+Ar [*ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 101, 2887 (1994).

<tr><td>arbha_mod   <td>none<td>BH(<i>A</i><sup>1</sup>&Pi;)+Ar [modified *ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 101, 2887 (1994).

<tr><td>arbhx   <td>none<td>BH(<i>X</i><sup>1</sup>&Sigma;<sup>+</sup>)+Ar [*ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 101, 2887 (1994).

<tr><td>arbhx_mod   <td>none<td>BH(<i>X</i><sup>1</sup>&Sigma;<sup>+</sup>)+Ar [modified *ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 101, 2887 (1994).

<tr><td>archb   <td>none<td>CH(<i>B</i><sup>2</sup>&Sigma;<sup>+</sup>)+Ar [*ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

<tr><td>archb_mod   <td>none<td>CH(<i>B</i><sup>2</sup>&Sigma;<sup>+</sup>)+Ar [modified *ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

<tr><td>archx   <td>none<td>CH(<i>X</i><sup>2</sup>&Pi;)+Ar [*ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

<tr><td>archx_mod   <td>none<td>CH(<i>X</i><sup>2</sup>&Pi;)+Ar [modified *ab initio*]

<tr><td><td><td>M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

<tr><td>arcnb_avqz_mrci_scaled  <td>none<td>CN(B)+Ar [avqz MRCI+Q scaled]

<tr><td><td><td>M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>arcnb_avtz_mrci_scaled  <td>none<td>CN(B)+Ar [avtz MRCI+Q scaled]

<tr><td><td><td>M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>arcnx_avqz_mrci_scaled  <td>none<td>CN(X)+Ar [avqz MRCI+Q scaled]

<tr><td><td><td>M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>arcnxb <td>none<td>CN(X<sup>2</sup>&Sigma;)+Ar [scaled *ab initio* MRCI+Q]

<tr><td><td><td>M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>arcnx_avtz_mrci_scaled  <td>none<td>CN(X)+Ar [avtz MRCI+Q scaled]

<tr><td><td><td>M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>arc2h2   <td>none<td>C<sub>2</sub>H<sub>2</sub>(<i>X</i><sup>1</sup>&Sigma;)+Ar [*ab initio* CCSDT(1.1)]

<tr><td><td><td>M. Yang, M. H. Alexander, H.-J. Werner, and R. Bemish, J. Chem. Phys. 105, 10462 (1996).

<tr><td>arh2o <td>none<td>H<sub>2</sub>O + Ar  [*ab initio* CCSD(T) PES]

<tr><td><td><td>J. Makarewicz, J. Chem. Phys. 129, 184310 (2008), fitted by P. Dagdigian.

<tr><td>arnha   <td>none<td>NH(<i>a</i><sup>1</sup>&Delta;)+Ar [original MRCI *ab initio* PES]

<tr><td><td><td>G. Jansen and B. A. Hess, Chem. Phys. Lett. 192, 21 (1992).

<tr><td>arnha_mod   <td>none<td>NH(<i>a</i><sup>1</sup>&Delta;)+Ar [original MRCI *ab initio* PES]

<tr><td><td><td>M. Yang, M. H. Alexander, C.-C. Chuang, R. W. Randall, and M. I. Lester, J. Chem. Phys. 103, 905 (1995).

<tr><td>arnhc   <td>none<td>NH(c<sup>1</sup>&Pi;)+Ar [original *ab initio* MRCI+Q PES's]

<tr><td><td><td>M. Yang, M. H. Alexander, C.-C. Chuang, R. W. Randall, and M. I. Lester, J. Chem. Phys. 103, 905 (1995).

<tr><td>arnhc_mod   <td>none<td>NH(c<sup>1</sup>&Pi;)+Ar [modified *ab initio* MRCI+Q PES's]

<tr><td><td><td>M. Yang, M. H. Alexander, C.-C. Chuang, R. W. Randall, and M. I. Lester, J. Chem. Phys. 103, 905 (1995).

<tr><td>arno    <td>none<td>NO(<i>X</i><sup>2</sup>&Pi;)+Ar [*ab initio*]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 99, 7725 (1993).

<tr><td>arno_ccsdt    <td>none<td>NO(<i>X</i><sup>2</sup>&Pi;)+Ar [*ab initio* CCSD(T)]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 111, 7426 (1999).

<tr><td>arno_2sg_klos_scaled<td>none<td>NO(A<sup>2</sup>&Sigma;)+Ar  [RCCSD(T)/AVTZ+33221 PES scaled with 1.23 factor to recover experimental D<sub>0</sub>]

<tr><td><td><td>J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, J. Chem. Phys. 129, 244303 (2008).

<tr><td>arn2<td>none<td>N<sub>2</sub>(<i>X</i><sup>1</sup>&Sigma;<sup>+</sup>)+Ar [model]

<tr><td><td><td>D. Pattengill, R. A. LaBudde, and R. B. Bernstein, J. Chem. Phys. 55, 5517 (1971).

<tr><td>aroha_ho_rabitz<td>none<td>OH(A<sup>2</sup>&Sigma;<sup>+</sup>)+Ar  

<tr><td><td><td>T.-S. Ho, H. Rabitz, S. E. Choi and M. I. Lester, J. Chem. Phys. 102, 2282 (1995), and J. Chem. Phys. 104, 1187 (1996).

<tr><td>arohx_rccsdt_corrected  <td>none<td>OH(<i>X</i><sup>2</sup>&Pi;)+Ar [RCCSD(T) *ab initio*]

<tr><td><td><td>G. Paterson, S. Marinakis, M. L. Costen, K. G. McKendrick, J. Klos, and R. Tobola, J. Chem. Phys. 129, 074304 (2008); erratum: G. Paterson et al., J. Chem. Phys. 131, 159901 (2009).

<tr><td>arohx_uccsdt_cbs_v0aver  <td>none<td>OH(<i>X</i><sup>2</sup>&Pi;)+Ar [*ab initio* complete basis set extrapolation, UCCSD(T) averaged over v=0 of OH]

<tr><td><td><td>L. Scharfenberg, J. Klos, P. J. Dagdigian, M. H. Alexander, G. Meijer, S. Y. T. van der Meerakker, Phys. Chem. Chem. Phys. 12, 10660 (2010).

<tr><td>arohx_ump4   <td>none<td>OH(<i>X</i><sup>2</sup>&Pi;)+Ar [UMP4 *ab initio*]

<tr><td><td><td>J. Klos, G. Chalasinski, M. T. Berry, R. A. Kendall, R. Burcl, and M. M. Szczesniak, J. Chem. Phys. 112, 4952 (2000).

<tr><td>aroh_2sg_re_klos and aroh_2sg_klos<td>none<td>OH(A<sup>2</sup>&Sigma;<sup>+</sup>)+Ar  [Klos's *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>J. Klos, M. H. Alexander, M. Brouard, C. J. Eyles, and F. J. Aoiz,  J. Chem. Phys.  129, 054301 (2008).

<tr><td>aroh_2sg_schnupf<td>none<td>OH(A<sup>2</sup>&Sigma;<sup>+</sup>)+Ar  [PES of Schnupf, Fawzy & Heaven]

<tr><td><td><td>W. Fawzy and M. C. Heaven, J. Chem. Phys. 89, 7030 (1988); U. Schnupf, J. M. Bowman, and M. C. Heaven, Chem. Phys. Lett. 189,487 (1992).

<tr><td>aroh_2sg_v0aver_klos<td>none<td>OH(A<sup>2</sup>&Sigma;<sup>+</sup>)+Ar  [Klos's *ab initio* RCCSD(T) PES averaged over <i>v</i>=0 of OH(A)]

<tr><td><td><td>J. Klos, M. H. Alexander, M. Brouard, C. J. Eyles, and F. J. Aoiz,  J. Chem. Phys.  129, 054301 (2008).

<tr><td>ar2 <td>none<td>Ar + Ar [model]

<tr><td><td><td>
M. H. Alexander, notes on quantum mechanical collision theory (see
www2.chem.umd.edu/groups/alexander/teaching)

<tr><td>bar   <td>none<td>B(2<i>p</i><sup>2</sup><i>P</i>)+Ar [modified <i> ab initio</i>] 

<tr><td><td><td>E. Hwang, Y.-L. Huang, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 98, 8484 (1993).

<tr><td>barb   <td>none<td>B(3<i>s</i><sup>2</sup><i>S</i>)+Ar [modified *ab initio*]

<tr><td><td><td>E. Hwang, Y.-L. Huang, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 98, 8484 (1993).

<tr><td>bh2   <td>none<td>B(<sup>2</sup>P)+H<sub>2</sub> [*ab initio* PES, Dubernet-Hutson expansion]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); J. Williams and M. H. Alexander, J. Chem. Phys. 112, 5722 (2000).

<tr><td>boh2   <td>none<td>B(<sup>2</sup>P)+<i>ortho</i>-H<sub>2</sub> [*ab initio* PES, Dubernet-Hutson expansion]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994); M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).

<tr><td>b3sh2   <td>none<td>B(2s<sup>2</sup>3s)+H<sub>2</sub> [*ab initio* PES, Dubernet-Hutson expansion]

<tr><td><td><td>M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995); M.-L. Dubernet et al., J. Chem. Phys. 94, 7602 (1991).

<tr><td>cahe    <td>none<td>Ca(4<i>s</i>5<i>p</i>) + He [*ab initio*]

<tr><td><td><td>B. Pouilly, J.-M. Robbe, and M. H. Alexander, J. Chem. Phys. 91, 1658 (1989); T. Duhoo and B. Pouilly, J. Chem. Phys. 101, 7554 (1994).

<tr><td>cd3he_v0avg <td>cd3he_ccsdt_v0avg.dat<td>CH<sub>3</sub>(v=0) + He  [*ab initio* CCSD(T) PES]



<tr><td><td><td>
M.H.Alexander and P.J.Dagdigian, J. Chem. Phys. 135, 064306 (2011);
Q.Ma, P.J.Dagdigian, and M.H.Alexander, J. Chem. Phys. 138, 104317 (2013).

<tr><td>chhe<td>chhe_abin.dat<td>CH(<i>X</i><sup>2</sup>&Pi;)+He [*ab initio*]

<tr><td><td><td>A. F. Wagner, T. H. Dunning, and R. A. Kok, J. Chem. Phys. 100, 1326 (1994); M. H. Alexander, W. Kearney, and A. F. Wagner, ibid. 100, 1338 (1994).

<tr><td>chhe_cyb<td>none<td>CH(<i>X</i><sup>2</sup>&Pi;)+He [*ab initio*]

<tr><td><td><td>S. M. Cybulski, G. Chalasinski, M. M. Szczesniak, J. Chem. Phys. 105, 9525 (1996).

<tr><td>chhe_driss_rkhs <td>none<td>CH(X<sup>2</sup>&Pi;)+He  [new D. Ben Abdallah's PES (corrected at short range unphysical behavior by J. Klos)]

<tr><td><td><td>D. Ben Abdallah et al, Phys. Chem. News 12, 1 (2003).

<tr><td>ch2he_a190_c20<td>ch2he_pot2.dat<td>CH<sub>2</sub>(a)+He  [*ab initio* PES]

<tr><td><td><td>L. Ma, M. H. Alexander, and P. J. Dagdigian, J. Chem. Phys. 134, 154307 (2011).

<tr><td>ch2he_x52_c20_v0 <td>ch2he_x52_v0.dat<td>CH<sub>2</sub>(a)+He  [*ab initio* PES, averaged over v=0 bend level] 

<tr><td><td><td>L. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 136, 224306 (2012).

<tr><td>ch2he_x52_c20_v1 <td>ch2he_x52_v1.dat<td>CH<sub>2</sub>(a)+He  [*ab initio* PES, averaged over v=1 bend level]

<tr><td><td><td>L. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 136, 224306 (2012).

<tr><td>ch2he_x52_c20_v2 <td>ch2he_x52_v2.dat<td>CH<sub>2</sub>(a)+He  [*ab initio* PES, averaged over v=2 bend level]

<tr><td><td><td>L. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 136, 224306 (2012).

<tr><td>ch2he_x52_c20_v3 <td>ch2he_x52_v3.dat<td>CH<sub>2</sub>(X)+He  [*ab initio* PES, averaged over v=3 bend level]

<tr><td><td><td>L. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 136, 224306 (2012).

<tr><td>ch3he_ccsdt <td>ch3he_pot.dat<td>CH<sub>3</sub>(X)+He  [*ab initio* CCSD(T) PES]

<tr><td><td><td>P. J. Dagdigian and M. H. Alexander, J. Chem. Phys. 135, 064306 (2011).

<tr><td>ch3he_vib <td>pot_ch3he_vib_ylmsym, pot_ch3he_vib_ylmasym, pot_ch3he_vib_data<td>CH<sub>3</sub>(X)+He  [*ab initio* CCSD(T) PES, with umbrella mode vib. levels]

<tr><td><td><td>Q. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 138, 104317 (2013).

<tr><td>ch3i   <td>none<td>CH<sub>3</sub>I

<tr><td><td><td>M. Shapiro, J. Phys. Chem. 90, 3644 (1986); H. Guo and G. C. Schatz, J. Chem. Phys. 93,
 393 (1990); H. Guo, K. Q. Lao, G. C. Schatz, and A. D. Hammerich, *ibid.* 94, 6562 (1991); 
M. H. Alexander, C. Rist, and D. E. Manolopoulos, *ibid.* 97, 4836 (1992).

<tr><td>clh2  <td>none<td>Cl(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Capecchi-Werner PES's]

<tr><td><td><td>G. Capecchi and H.-J. Werner, Phys. Chem. Chem. Phys. 6, 4975 (2004).

<tr><td>clh2_aquilanti <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion]

<tr><td><td><td>V. Aquilanti, D. Cappelletti, S. Cavalli, F. Pirani, and A. Volpi, J. Phys. Chem. A 105, 2401 (2001).

<tr><td>clh2_rdep  <td>none<td>Cl(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Werner expansion of r-dependent PES's]

<tr><td><td><td>G. Capecchi and H.-J. Werner, Phys. Chem. Chem. Phys. 6, 4975 (2004).

<tr><td>clh2_rdep_negion  <td>none<td>Cl<sup>-</sup>(<sup>1</sup>S)+H<sub>2</sub> [Dubernet-Werner expansion of r-dependent PES]

<tr><td><td><td>G. Capecchi and H.-J. Werner, Phys. Chem. Chem. Phys. 6, 4975 (2004).

<tr><td>clh2minus_ccsdt  <td>none<td>Cl<sup>-</sup>(<sup>1</sup>S)+H<sub>2</sub> [avqz CCSD(T) *ab initio* PES]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 118, 9637 (2003).

<tr><td>clmin_h2_re_klospes<td>none<td>Cl<sup>(-)</sup>+H<sub>2</sub> [*ab initio* RCCSD(T) PES's]

<tr><td><td><td>A. A. Buchachenko, T. A. Grinev, J. Klos, E. J. Bieske, M. M. Szczesniak and G. Chalasinski, J. Chem. Phys. 119, 12931, (2003); M. H. Alexander, J. Klos and D. Manolopoulos, J. Chem. Phys. 128, 084312 (2008).

<tr><td>cnarb_us_1sig <td>none<td>CN(B)+Ar 

<tr><td><td><td>J. Han, M. C. Heaven, U. Schnupf, and M. H. Alexander, J. Chem. Phys. 128, 104308 (2008).

<tr><td>cnarxa_new, cnarx_refit, cnarxa_73, cnar_pi  <td>none<td>CN(A v=3,X v=7)-Ar [refit of Berning's original PES's]

<tr><td><td><td>A. Berning (1995)  Ph. D. thesis, Universitat Stuttgart.; M. H. Alexander, X. Yang, P. J. Dagdigian, A. Berning, and H.-J. Werner, J. Chem. Phys. 112, 781 (2000).

<tr><td>cnarx_us_1sig <td>none<td>CN(X)+Ar  

<tr><td><td><td>J. Han, M. C. Heaven, U. Schnupf, and M. H. Alexander, J. Chem. Phys. 128, 104308 (2008).

<tr><td>cnar_A_dbennett_v0  
<td>CNAAr_sum.data.tab, CNAAr_sum_coeff.tab, CNAAr_diff.data.tab, CNAAr_diff.coeff.tab
<td>CN(A<sup>2</sup>&Pi;)+Ar [original *ab initio* UCCSD(T) PES's - r averaged over v=0 probability distribution]

<tr><td><td><td>
S. J, McGurk, K. G. McKendrick, M. L. Costen, D.I.G. Bennett,J. Klos, M. H. Alexander, P. J. Dagdigian, J. Chem. Phys. 136, 164306 (2012).

<tr><td>cnar_A_dbennett_v4  
<td>CNAAr_sum.data.tab, CNAAr_sum_coeff.tab, CNAAr_diff.data.tab, CNAAr_diff.coeff.tab
<td>CN(A<sup>2</sup>&Pi;)+Ar [original *ab initio* UCCSD(T) PES's - r equals average value for v=4]

<tr><td><td><td>
S. J, McGurk, K. G. McKendrick, M. L. Costen, D.I.G. Bennett,J. Klos, M. H. Alexander, P. J. Dagdigian, J. Chem. Phys. 136, 164306 (2012).

<tr><td>cnar_X_dbennett  <td>CNXAr.data.tab, CNXAr.coeff.tab
<td>CN(X<sup>2</sup>&Sigma;)+Ar [original *ab initio* UCCSD(T) PES's]

<tr><td><td><td>
S. J, McGurk, K. G. McKendrick, M. L. Costen, D.I.G. Bennett,J. Klos, M. H. Alexander, P. J. Dagdigian, J. Chem. Phys. 136, 164306 (2012).


<tr><td>cnne  <td>none<td>CN(<sup>2</sup>&Sigma;-<sup>2</sup>&Pi;)+Ne &Sigma;-&Pi; atom scattering 

<tr><td><td><td>M. Yang and M. H. Alexander, J. Chem. Phys. 107, 7148 (1997).

<tr><td>cnnn_bennett  <td>none <td>CN(A<sup>2</sup>&Pi;)+N<sub>2</sub>  [UCCSD(T), v=3 geom for Rcn]

<tr><td><td><td>A. Khachatrian, P. J. Dagdigian, D. I. G. Bennett, F. Lique,  J. Klos and M. H. Alexander, J. Phys. Chem. A  113, 3922 (2009).

<tr><td>cph2 <td>none<td>C<sup>+</sup>(<sup>2</sup>P) + <i>p</i>,<i>o</i>-H<sub>2</sub>  [*ab initio* MRCI+Q PES's]

<tr><td><td><td>
F. Lique, G. Werfelli, P. Halvick, T. Stoecklin, A. Faure, L. Wiesenfeld,
and P. J. Dagdigian, J. Chem. Phys. 138, 204314 (2013).

<tr><td>fh2_aquilanti <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion]

<tr><td><td><td>V. Aquilanti, D. Cappelletti, S. Cavalli, F. Pirani, and A. Volpi, J. Phys. Chem. A 105, 2401 (2001).

<tr><td>fh2_asw <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Alexander-Stark-Werner PES's]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).

<tr><td>fh2_lwa and fh2_lwal <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Li-Werner-Alexander PES's]

<tr><td><td><td>G. Li, H.-J. Werner, F. Lique, and M. H. Alexander, J. Chem. Phys. 127, 174302 (2007); M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).

<tr><td>fh2_sw <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of Stark-Werner PES's r=1.4]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).

<tr><td>hcl<td>none<td>HCl (<i>X, a, b, </i><sup>3</sup>&Sigma;<sup>+</sup>) [*ab initio*]

<tr><td><td><td>M. H. Alexander, B. Pouilly, and T. Duhoo, J. Chem. Phys. 99, 1752 (1993).

<tr><td>hcl_a3pi and hcl_a3pinew  <td>none<td>HCl(X,A,a,triplet)

<tr><td><td><td>M. H. Alexander, B. Pouilly, and T. Duhoo,  J. Chem. Phys. 99, 1752 (1993).

<tr><td>hcl_t3sig   <td>none<td>H(<sup>2</sup>S)+Cl(<sup>2</sup>P)

<tr><td><td><td>M. H. Alexander, B. Pouilly, and T. Duhoo,  J. Chem. Phys. 99, 1752 (1993).

<tr><td>hcl_vib and hcl_vib_new  <td>none<td>H(<sup>2</sup>S)+Cl(<sup>2</sup>P)

<tr><td><td><td>M. H. Alexander, B. Pouilly, and T. Duhoo,  J. Chem. Phys. 99, 1752 (1993).

<tr><td>hecn_dgels  <td>hecn_fitmlv.dat <td>CN(X<sup>2</sup>&Sigma;)+He  *ab initio* RCCSD(T) PES

<tr><td><td><td>F. Lique, A. Spielfiedel, N. Feautrier, I. F. Schneider, J. Klos, and M. H. Alexander, J. Chem. Phys. 132, 024303, 2010.

<tr><td>heco_sapt  <td>none<td>CO(<sup>1</sup>&Sigma;)+He [SAPT PES]

<tr><td><td><td>
T. G. A. Heijmen, R. Moszynski, P. E. S. Wormer, and A. van der Avoird,
J. Chem. Phys. 107, 9921 (1997).

<tr><td>heno_klos <td>none<td>NO(X<sup>2</sup>&Pi;) + He  [original *ab initio* RCCSD(T) PES's]

<tr><td><td><td>
J. Klos, G. Chałasiński, M. T. Berry, R. Bukowski, and S. M. Cybulski, 
J. Chem. Phys. 112, 2195 (2000).

<tr><td>heno_2sg_klos<td>none<td>NO(A<sup>2</sup>&Sigma;)+He  [original *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, J. Chem. Phys. 129, 244303 (2008).

<tr><td>henp_klos_lique <td>henp_fitmlv.dat<td>PN(X<sup>1</sup>&Sigma;)+He  [original *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>R. Tobola, J. Klos, F. Lique, G. Chalasinski, and M. H. Alexander, Astronomy & Astrophysics, 468 (3), 1123-1127 (2007).

<tr><td>heohx_lmax10<td>none<td>OH(X<sup>2</sup>&Pi;)+He  [Cybulski's *ab initio* RCCSD(T) PES]

<tr><td><td><td>H.-S. Lee, A. B. McCoy, R. R. Toczylowski, and S. M. Cybulski, J. Chem. Phys. 113, 5736 (2000).

<tr><td>heohx_3d<td>pot_heohx_3d_data<td>OH(X<sup>2</sup>&Pi;)+He  [*ab initio* avtz-av5z RCCSD(T) PES, averaged over OH(X,v=0) wavefunction]

<tr><td><td><td>K. B. Gubbels, Q. Ma, M. H. Alexander, P. J. Dagdigian, D. Tanis, G. C. Gronenboom, A. van der Avoird, and S. Y. T. van de Meerakker, J. Chem. Phys. 136, 144308 (2012).

<tr><td>hepha  <td>none<td>PH(A<sup>3</sup>&Pi;)+He [original *ab initio* CEPA PES's ]

<tr><td><td><td>Ch. Kolczewski, K. Fink, V. Staemmler, and L. Neitsch, J. Chem. Phys. 106, 7637 (1997).

<tr><td>hepha_red  <td>none<td>PH(A<sup>3</sup>&Pi;)+He [PES modified to include just v00, v20, v22]

<tr><td><td><td>Ch. Kolczewski, K. Fink, V. Staemmler, and L. Neitsch, J. Chem. Phys. 106, 7637 (1997).

<tr><td>heshx_cybulski <td>none<td>SH(X<sup>2</sup>&Pi;)+He  [original *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>S. M. Cybulski, R. R. Toczylowski, H.-S. Lee, and A. B. McCoy, J. Chem. Phys. 113, 9549 (2000).

<tr><td>heno_2sg_klos <td>none<td>NO(A<sup>2</sup>&Sigma;) + He  [*ab initio* RCCSD(T) PES]

<tr><td><td><td>
J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, 
J. Chem. Phys. 129, 244303 (2008).

<tr><td>hfhf   <td>none<td>HF(X<sup>1</sup>&Sigma;)+HF(X<sup>1</sup>&Sigma;)

<tr><td><td><td>M. H. Alexander and A. E. DePristo, J. Chem. Phys. 65, 5009 (1976); P. F. Vohralik, R. O. Watts, and M. H. Alexander, J. Chem. Phys. 91, 7563 (1989); J. Chem. Phys. 93, 3983 (1990).

<tr><td>h2f_klos_ad <td>none<td>F(<sup>2</sup>P)+H<sub>2</sub> [Klos's CCSD(T) adiabats]

<tr><td><td><td>J. Klos, G. Chalasinski and M. M. Szczesniak, Int. J. Quant. Chem. 90, 1038 (2002).

<tr><td>h2ohe  <td>h2o_coefd.dat, h2o_coefi.dat, h2o_params.dat<td>H<sub>2</sub>O+He [*ab initio* SAPT PES]

<tr><td><td><td>K. Patkowski, T. Korona, R. Moszyniski, B. Jeziorski, and K. Szalewicz, J. Molec. Struct. (Theochem) 591, 231 (2002), fitted by P. Dagdigian.

<tr><td>h2a <td>none<td>H<sub>2</sub>(<sup>3</sup>&Sigma;<sub>u</sub>)  [*ab initio* MRCI]

<tr><td><td><td>P. J. Dagdigian, unpublished.

<tr><td>h2oh_c34 <td>none<td>H<sub>2</sub>O+H  [*ab initio* RCCSD(T) theory]

<tr><td><td><td>
P. J. Dagdigian and M. H. Alexander, to be published.

<tr><td>h2x <td>none<td>H<sub>2</sub>(X<sup>1</sup>&Sigma;<sub>g</sub>)  [*ab initio* MRCI]

<tr><td><td><td>P. J. Dagdigian, unpublished.

<tr><td>krno_klos<td>none<td>NO(X<sup>2</sup>&Pi;)+Kr  [Klos's *ab initio* RCCSD(T) PES]

<tr><td><td><td>B. Wen, H. Meyer, J. Klos, and M. H. Alexander, J. Phys. Chem. A 113, 7366 (2009); S. Marinakis, B. J. Howard, F. J. Aoiz, J. Klos, Chem. Phys. Lett. 512, 161 (2011).

<tr><td>krno_2sg_dkroll_klos <td>none<td>NO(A<sup>2</sup>&Sigma;)+Kr  [original *ab initio* RCCSD(T)+Douglas-Kroll PES for r=r<sub>e</sub>]

<tr><td><td><td>J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, J. Chem. Phys. 129, 244303 (2008).

<tr><td>kroh_2sg_re_klos<td>none<td>OH(A<sup>2</sup>&Sigma;)+Kr  [Klos's *ab initio* RCCSD(T) PES]

<tr><td><td><td>H. Chadwick, M. Brouard, Y.-P. Chang, C. J. Eyles, T. Perkins, S. A. Seamons, J. Klos, M. H. Alexander, and F. J. Aoiz, J. Chem. Phys. 137, 154305 (2012).

<tr><td>kroh_2sg_re_klos_rkhs2d_anyvllmax<td>none<td>OH(A<sup>2</sup>&Sigma;)+Kr  [Klos's *ab initio* RCCSD(T) PES]

<tr><td><td><td>H. Chadwick, M. Brouard, Y.-P. Chang, C. J. Eyles, T. Perkins, S. A. Seamons, J. Klos, M. H. Alexander, and F. J. Aoiz, J. Chem. Phys. 137, 154305 (2012).

<tr><td>li2ne  <td>LI2NE.MHA<td>Li<sub>2</sub>(A)+Ne 

<tr><td><td><td>M. H. Alexander and H.-J. Werner, J. Chem. Phys. 95, 6524 (1991).

<tr><td>he_nco_klos_tobola_11angles <td>none<td>NCO(X<sup>2</sup>&Pi;)+He  [original *ab initio* RCCSD(T) PES's]

<tr><td><td><td>J. Klos, R. Tobola and G. Chalasinski, J. Phys. Chem. A,  113, 14480 (2009).

<tr><td>nenhc   <td>none<td>NH(c<sup>1</sup>&Pi;)+Ne [*ab initio* potential]

<tr><td><td><td>G. Kerenskaya, U. Schnupf, and M. C. Heaven, J. Chem. Phys. 119, 8424 (2003).

<tr><td>nenh3pi   <td>none<td>NH(A<sup>3</sup>&Pi;)+Ne [*ab initio* potential]

<tr><td><td><td>G. Kerenskaya, U. Schnupf, M. C. Heaven, A. van der Avoird, and G. C. Groenenboom, Phys. Chem. Chem. Phys. 7, 846 (2005).

<tr><td>neno_ccsdt <td>none<td>NO(X<sup>2</sup>&Pi;)+Ne  [CCSD(T) avqz]

<tr><td><td><td>M. H. Alexander, P. Soldan, T. G. Wright, Y. Kim, H. Meyer, P. J. Dagdigian, and E. P. F. Lee, J. Chem. Phys. 114, 5588 (2001).

<tr><td>neno_2sg_rydb_ramon <td>none<td>NO(A<sup>2</sup>&Sigma;)+Ne  [R. Hernandez's *ab initio* PES]

<tr><td><td><td>P. Pajon-Suarez, G. Rojas-Lorenzo, J. Rubayo-Soneira, R. Hernandez-Lamoneda, Chem. Phys. Lett. 421,389-394 (2006).

<tr><td>neohx   <td>none<td>OH(<i>X</i><sup>2</sup>&Pi;)+Ne [*ab initio*]

<tr><td><td><td>M. Yang and M. H. Alexander, J. Chem. Phys. 103, 3400 (1995).

<tr><td>nhhe   <td>none<td>NH(<i>A</i><sup>3</sup>&Pi;)+He [*ab initio*]

<tr><td><td><td>R. Jonas and V. Staemmler, Z. Phys. D 14, 143 (1989); M. H. Alexander, P. J. Dagdigian, and D. Lemoine, J. Chem. Phys. 95, 5036 (1991).

<tr><td>noag_cem   <td>none<td>NO(X<sup>2</sup>&Pi;)+Ag(111) [CEM PES's]

<tr><td><td><td>A. E. DePristo and M. H. Alexander, J. Chem. Phys. 94, 8454 (1991).

<tr><td>nohe<td>none<td>NO(<i>X</i><sup>2</sup>&Pi;)+He [*ab initio*]

<tr><td><td><td>M. Yang and M. H. Alexander, J. Chem. Phys. 103, 6973 (1995).

<tr><td>nohe_klos<td>none<td>NO(X<sup>2</sup>&Pi;)+He  [original *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>J. Klos, G. Chalasinski, R. Bukowski, S.M. Cybulski, and M. T. Berry,  J. Chem. Phys. 112, 2195 (2000).

<tr><td>no_closed_he_klos<td>none<td>NO(X<sup>2</sup>&Pi;)+He  [original *ab initio* RCCSD(T) PES for r=r<sub>e</sub>, only V<sub>sum</sub> NO molecule treated as closed shell]

<tr><td><td><td>J. Klos, G. Chalasinski, R. Bukowski, S.M. Cybulski, and M. T. Berry,  J. Chem. Phys. 112, 2195 (2000).

<tr><td>nh3he_wheatley  <td>none<td>NH<sub>3</sub>(X)+He [scaled perturbation theory PES]

<tr><td><td><td>M. P. Hodges and R. J. Wheatley, J. Chem. Phys. 14 , 8836 (2001).

<tr><td>nh3h2_1993  <td>fitvij.h2<td>NH<sub>3</sub>(X)+H<sub>2</sub> [*ab initio* PES, Valiron expansion of the potential]

<tr><td><td><td>C. Rist, M. H. Alexander, and P. Valiron, J. Chem. Phys. 98, 4662 (1993).

<tr><td>nh3h2_2009  <td>fitvij_bf_62.h2<td>NH<sub>3</sub>(X)+H<sub>2</sub> [*ab initio* PES, Valiron expansion of the potential]

<tr><td><td><td>S. Maret, A. Faure, E. Scifoni, and L. Wiesenfeld, Mon. Not. R. Soc. 399, 425 (2009).

<tr><td>ohh_trip  <td>none<td>OH(X<sup>2</sup>&Pi;)+H(<sup>2</sup>S) [triplet potentials from  MRCI calculations]

<tr><td><td><td>M. H. Alexander, E. J. Rackham, and D. E. Manolopoulos, J. Chem. Phys. 121, 5221 (2004).

<tr><td>o3ph2   <td>none<td>O(<sup>3</sup>P)+H<sub>2</sub> [Dubernet-Hutson expansion of scaled Alexander PES's ]

<tr><td><td><td>M. H. Alexander, J. Chem. Phys. 108, 4467 (1998); M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994); M. H. Alexander, J. Chem. Phys. 99, 6014 (1993); M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).

<tr><td>ohh2   <td><i>To be specified in input files</i><td>General potential subroutine for the collision between a<sup>2</sup>&Pi; molecule and a <sup>1</sup>&Sigma; molecule (or a structureless atom)
<tr><td><td><td>OH+H<sub>2</sub>: To be published. [use the data files <i>pot_ohh2_mrci_core1.dat</i> for the MRCI PES, <i>pot_ohh2_mrci_core1_highsym.dat</i> for the MRCI High-symmetry PES, and <i>pot_ohh2_ccsdf12-avtzbf.dat</i> for the CCSD(T) PES]


<tr><td>stp1sg_qma   <td><i>To be specified in input files</i><td>General potential subroutine for the collision between a symmetric top molecule and a <sup>1</sup>&Sigma; molecule (or a structureless atom)
<tr><td><td><td>CH<sub>3</sub>+H<sub>2</sub>: To be published.<br>
NH<sub>3</sub>+H<sub>2</sub>: S. Maret, A. Faure, E. Scifoni, and L. Wiesenfeld, Mon. Not. R. Soc. 399, 425 (2009). (use the data file pot_nh3h2_cp.dat)

<tr><td>vfit   <td>N2PHXHE.BIN<td>N<sub>2</sub><sup>+</sup>+He

<tr><td><td><td>B. Follmeg, P. Rosmus, and H.-J. Werner, J. Chem. Phys. 93, 4687 (1990).

<tr><td>xeno_klos<td>none<td>NO(X<sup>2</sup>&Pi;)+Xe  [Klos's *ab initio* RCCSD(T) PES for r=r<sub>e</sub>]

<tr><td><td><td>
Theoretical studies of bound states and scattering of Xe atoms with NO(X) radical. J. Klos, F. J. Aoiz, M. H. Alexander, in preparation 2011.

<tr><td>xeno_klos<td>none<td>NO(X<sup>2</sup>&Pi;)+Xe  [*ab initio* RCCSD(T)/ECP PES's]

<tr><td><td><td>J. Klos, F. J. Aoiz, M. Menendez, M. Brouard, H. Chadwick, and C. J. Eyles, J. Chem. Phys. 137, 014312 (2012).

<tr><td>xeno_2sg_klos<td>none<td>NO(A<sup>2</sup>&Sigma;)+Xe  [original *ab initio* RCCSD(T)/ECP PES for r=r<sub>e</sub>]

<tr><td><td><td>J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, J. Chem. Phys. 129, 244303 (2008).


</table>