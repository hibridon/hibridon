      subroutine version(iunit)
        implicit none
        integer :: iunit, i
        character(len=100) :: txt

        write(iunit,10) _BUILD_VERS_
10      format (/' ---------------------------------------------------',
     :  '-----------------------',/
     :   ,9x,
     :   '  HIBRIDON SCATTERING CODE V ',f4.2,//
     :   '      AUTHORS: M. ALEXANDER, D. MANOLOPOULOS,',
     :   ' H.-J. WERNER, B. FOLLMEG,'/
     :   '               P. DAGDIGIAN',/
     :   ' CONTRIBUTORS: D. LEMOINE, P. VOHRALIK,',
     :   ' G. COREY, R. JOHNSON, T. ORLIKOWSKI,',/
     :   '          A. BERNING, A. DEGLI-ESPOSTI,',
     :   ' C. RIST, B. POUILLY, J. KLOS, Q. MA,',/,
     :   '          G. VAN DER SANDEN, M. YANG, F. DE WEERD',
     :   ', S. GREGURICK, F. LIQUE',/
     :   ' ---------------------------------------------------',
     :   '-----------------------')

       write (iunit,15)
15     format(11x,'BUILD CONFIGURATION:')
        write(iunit,*)
        write(iunit,'(a18,a)') 'CMAKE BUILD TYPE: ', _BUILD_TYPE_
        write(iunit,'(a5,a,a1,a)')  'SYS: ', _BUILD_SYS1_ , ' ', _BUILD_SYS2_
        write(iunit,'(a11,a)') 'PROCESSOR: ', _BUILD_PROC_
        write(iunit,'(a10,a)') 'COMPILER: ', _BUILD_COMP_
        write(iunit,'(a17,a)') 'COMPILE OPTIONS: ', _BUILD_FLAGS_
        write (iunit,20)
20     format (' ---------------------------------------------------',
     :  '-----------------------')

      write (iunit,30)
30    format(11x,' CURRENT HARDWARE CONFIGURATION:')
        write(iunit,*)
        if (_BUILD_SYS1_ == "Darwin") then
            call system('system_profiler SPHardwareDataType 
     :      | grep -E "Model|Processor|Cores|Cache|Hyper|Memory" >/tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            do 
                read(666,'(a)', iostat=i) txt
                if(i.ne.0) exit
                write(iunit,'(a)') trim(txt)
            enddo 
            close(666)
      elseif (_BUILD_SYS1_ == "Linux") then
        
        call system('cat /sys/devices/virtual/dmi/product_name > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a12,a)') "Model Name: ", trim(txt)
        call system('cat /sys/devices/virtual/dmi/product_sku > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a18,a)') "Model Identifier: ", trim(txt)
        call system('lscpu | grep -m 1 "Model name" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a16,a)') "Processor Name: ", trim(txt)    
        call system('lscpu | grep -m 1 "CPU MHz" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a17,a)') "Processor speed: ", trim(txt)    
        call system('lscpu | grep -m 1 "Socket(s)" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a22,a)') "Number of Processors: ", trim(txt)   
        call system('lscpu | grep -m 1 "Core(s)" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a23,a)') "Total Number of Cores: ", trim(txt) 
        call system('lscpu | grep -m 1 "L2" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a21,a)') "L2 Cache (per Core): ", trim(txt) 
        call system('lscpu | grep -m 1 "L3" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a10,a)') "L3 Cache: ", trim(txt) 
        call system('cat /proc/cpuinfo | grep -m 1 -o "ht" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt 
            close(666)
        if(trim(txt)=='ht')write(iunit,'(a35)') "Hyper-Threading Technology: Enabled" 
        call system('cat /proc/meminfo | grep "MemTotal" > /tmp/hw.info')
            open(unit=666, file="/tmp/hw.info", status="old")
            read(666,'(a)', iostat=i) txt, txt
            close(666)
        write(iunit,'(a8,a)') "Memory: ", trim(txt) 
      endif
      write(iunit,31)
31    format(
     : ' ---------------------------------------------------',
     :  '-----------------------')

      end subroutine version

      subroutine acknow(iunit,ipos)
        implicit none
        integer :: iunit
        logical :: ipos
        if (ipos) write (iunit, 10)
10    format(/' --------------------------------------------------------',
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
     :  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,',
     :  ' J. Klos, and F. Lique',
     : /,' --------------------------------------------------------',
     :  '-------------',
     :  '----------------------------------------------------------')
      if (.not.ipos) write (iunit, 20)
20    format(/' --------------------------------------------------------',
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
     :  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,',
     :  ' J. Klos, and F. Lique',
     : /,' --------------------------------------------------------',
     :  '---------------------')
      end subroutine acknow