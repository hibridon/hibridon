cstart unix-ibm unix-aix
@process noopt
cend
cstart unix-hp
c;!$hp$optimize off
cend
      subroutine version(iunit)
* to output version number of hibridon code
* current revision date:  4-apr-2003
      write (iunit, 10)
10    format
     : (/' ---------------------------------------------------',
     :  '-----------------------',
     :  /,9x,
     :   '  HIBRIDON SCATTERING CODE V xdate',
     : //'     AUTHORS: M. ALEXANDER, D. MANOLOPOULOS,',
     :   ' H.-J. WERNER, B. FOLLMEG',
     :  /' CONTRIBUTORS: D. LEMOINE, P. VOHRALIK,',
     :  ' G. COREY, R. JOHNSON, T. ORLIKOWSKI',
     :  /'               A. BERNING, A. DEGLI-ESPOSTI,',
     :  ' C. RIST, P. DAGDIGIAN, B. POUILLY',/,
     :  '               G. VAN DER SANDEN, M. YANG, F. DE WEERD',
     :  ', S. GREGURICK',
     : /' ---------------------------------------------------',
     :  '-----------------------')
      return
      end
cstart unix-ibm unix-aix
@process noopt
cend
cstart unix-hp
c;!$hp$optimize off
cend
      subroutine acknow(iunit,ipos)
* to acknowledge authors of hibridon
* current revision date:  22-apr-1997
      logical ipos
      if (ipos) write (iunit, 10)
10    format
     : (/' --------------------------------------------------------',
     :  '----------------------------------------------------------',
     :  '-------------',
     :  /,' All publications resulting from use of the integrators',
     :  ' included in the Hibridon code must include',
     :   ' the following reference:',
     : //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);',
     :  ' M. H. Alexander and D. E. Manolopoulos, J. Chem. Phys.',
     :  ' 80, 2044 (1987).',
     :  //,' All publications resulting from use of the',
     :  ' Hibridon package must include',
     :  ' the following reference:',
     : //,' HIBRIDON is a package of programs for the time-independent',
     :  ' quantum treatment',
     :  ' of inelastic collisions and photodissociation',
     :  /,' written by',
     :  ' M. H. Alexander, D. E.  Manolopoulos, H.-J. Werner, and',
     : ' B. Follmeg,',
     :  ' with contributions by',
     :  ' P. F. Vohralik, D. Lemoine,',
     :  /,' G. Corey, B. Johnson, T. Orlikowski, A. Berning,',
     :  ' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,',
     :  ' G. van der Sanden, M. Yang, F. de Weerd, and S. Gregurick',
     : /,' --------------------------------------------------------',
     :  '-------------',
     :  '----------------------------------------------------------')
      if (.not.ipos) write (iunit, 20)
20    format
     : (/' --------------------------------------------------------',
     :  '---------------------',
     :  /,' All publications resulting from use of the integrators',
     :  ' included in the',/,' Hibridon code must include',
     :   ' the following reference:',
     : //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);',
     :  ' M. H. Alexander and D. E.',/,
     :    ' Manolopoulos, J. Chem. Phys.',
     :  ' 80, 2044 (1987).',
     :  //,' All publications involving the determination of',
     :     'photodissociation cross sections',/,' must also include',
     :     ' the following reference:',
     :  //,' M. H. Alexander, Comput. Phys. Commun, 75, 87 (1993).',
     :  //,' All publications investigating flux redistribution must',
     :  ' include the following references:',
     :  //, ' M. H. Alexander, J. Chem. Phys. 95, 8931 (1991); D. E.',
     :      ' Manolopoulos and',/,' M. H. Alexander, J. Chem. Phys.',
     :      '  97, 2527 (1992).',
     :  //,' All publications resulting from use of the',
     :  ' Hibridon package must include',
     :  /,' the following reference:',
     : //,' HIBRIDON is a package of programs for the time-independent',
     :  ' quantum treatment',
     :  /,' of inelastic collisions and photodissociation',
     :  ' written by',
     :  ' M. H. Alexander,',/,' D. E.  Manolopoulos, H.-J. Werner, and',
     : ' B. Follmeg,',
     :  ' with contributions by',
     :  /,' P. F. Vohralik, D. Lemoine,',
     :  ' G. Corey, B. Johnson, T. Orlikowski, A. Berning,',
     :  /,' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,',
     :  ' G. van der Sanden, M. Yang, F. de Weerd, and S. Gregurick',
     : /,' --------------------------------------------------------',
     :  '---------------------')
      return
      end
