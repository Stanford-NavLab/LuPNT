c iritest.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c
c test program for the iri_web subroutine
c
c corrections
c-version-mm/dd/yy------------------------------------------------------
c 2000.01 05/07/01 initial version
c 2000.02 07/11/01 line 210: do i=1,100 instead of i=2,100 (K. Tokar)
c 2000.03 28/12/01 output oar(39) for IG12 (R. Conkright, NGDC, NOAA)
c 2000.04 28/10/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
c 2000.05 02/06/03 Ne(Te) only 300,400; foF1 and hmF1 output corr.
c 2000.06 01/19/05 (.not.jf(20)) instead of (..jf(2)) (G. Schiralli)
c 2005.01 05/06/06 included spread-F (jf(28)) and topside (jf(29)) options
c 2007.00 05/18/07 Release of IRI-2007
c
      INTEGER           pad1,pad2,pad3
      DIMENSION         outf(20,100),oar(50,100),jfi(6)
      LOGICAL		    jf(30)
      CHARACTER*2       timev(2)
      CHARACTER*3       uni(48),sopt
      CHARACTER*4       IMZ(8),MAP,xtex,coorv(2)
      CHARACTER*5       ITEXT(8)
      CHARACTER*6       pna(48)
      CHARACTER*8       bopt
      CHARACTER*9       topt,pname(6)
      CHARACTER*10      iopt
      CHARACTER*16      f1opt

      DATA  IMZ  /' km ','GEOD','GEOD','yyyy',' mm ',' dd ','YEAR',
     &      'L.T.'/, ITEXT/'  H  ',' LATI',
     &      ' LONG',' YEAR','MONTH',' DAY ','DAYOF',' HOUR'/

      data pna/'NmF2','hmF2','NmF1','hmF1','NmE','hmE','NmD','hmD',
     &  'h05','B0','NVmin','hVtop','Tpeak','hTpek','T300','T400','T600',
     &  'T1400','T3000','T120','Ti450','hTeTi','sza','sndec','dip',
     &  'dipla','modip','dela','Srise','Sset','seasn','nseas','Rz12',
     &  'cov','B1','M3000','TEC','TECtp','IG12','F1prb','F107d',
     &  'C1','daynr','vdrft','foF2r','F1noL','F1+L','sp_F'/
      data uni/'m-3','km','m-3','km','m-3','km','m-3','km','km','km',
     &   'm-3','km','K','km',7*'K','km',5*'deg',' ',2*'h',6*' ','m-2',
     &   '%',5*' ','m/s',4*' '/,timev/'LT','UT'/,coorv/'geog','geom'/

        data jfi/8,9,13,14,15,16/

c user input of IRI input parameters
c
1       print *,'jmag(=0/1,geog/geom),lati/deg,long/deg'
        read(5,*) jm,xlat,xlon
        print *,'year(yyyy),mmdd(or -ddd),iut(=0/1,LT/UT),hour'
        read(5,*) iy,imd,iut,hour
        print *,'height/km'
        print *,'(enter  0 for list of peak heights and densities)'
        print *,'(enter -1 for plasma frequencies, B0, M3000, ',
     &                        'valley, width and depth,)'
        print *,'(         F1 probability, equatorial vertical ',
     &         'ion drift, and)'
        print *,'(         foF2 storm/quiet ratio,)'
        print *,'(         or 3 parameter of your choice)'
        read(5,*) hx
        print *,'upper height [km] for TEC integration (0 for no TEC)'
        read(5,*) htec_max

        print *,'variable? (1/2/../8 for height/lat/long/year/month/',
     &                        'day/day of year/hour)'
        read(5,*) ivar
        print *,'begin, end, and stepsize for the selected variable'
        read(5,*) vbeg,vend,vstp

        print *,'Options: t(rue) or f(alse)'
        print *,'Standard: t,t,t,t,f,f,t,t,t,t,t,t,t,t,t,t,t,t,t,t,',
     &                        'f,t,f,t,t,t,t,f,f,f'
        print *,'Enter 0 to use standard or 1 to enter your own'
        read(5,*) jchoice
          do i=1,30
                jf(i)=.true.
                enddo
        if(jchoice.eq.0) then
c          jf(2)=.false.				  ! no temperatures
c          jf(3)=.false.				  ! no ion composition
          jf(5)=.false.               ! URSI foF2 model
          jf(6)=.false.               ! Newest ion composition model
          jf(21)=.false.			  ! no ion drift
          jf(23)=.false.              ! TTS Te model is standard
          jf(28)=.false.			  ! no spread-F
          jf(29)=.false.              ! New Topside options
          jf(30)=.false.              ! NeQuick topside
        else
          print *,'Compute Ne, T, Ni? (enter: t,t,t  if you want all)'
          read(5,*) jf(1),jf(2),jf(3)
        if(jf(1)) then
              print *,'LAY version: t=standard ver., f=LAY version. {t}'
              read(5,*) jf(11)
              print *,'Ne Topside: t=IRI-2001, f=new options {t}'
              read(5,*) jf(29)
              print *,'Ne Topside: t=IRI01_corrt, f=NeQuick {t}'
              read(5,*) jf(30)
              print *,'Ne Topside: t=F10.7<188, f=unlimited {t}'
              read(5,*) jf(7)
              print *,'foF2 model: t=CCIR, f=URSI-88 {standard:f}'
              read(5,*) jf(5)
              print *,'foF2: t=with storm model, f=without {t}'
              read(5,*) jf(26)
              print *,'F2 peak density or foF2: t=model, ',
     &              'f=user input {t}'
              read(5,*) jf(8)
              print *,'F2 peak height or M3000F2: t=model, ',
     &              'f=user input {t}'
              read(5,*) jf(9)
              print *,'Bottomside thickness B0: t=Table-option, ',
     &            'f=Gulyaeva {t}.'
              read(5,*) jf(4)
              print *,'F1 peak density or foF1: t=model, ',
     &            'f=user input {t}'
              read(5,*) jf(13)
            if(.not.jf(11)) then
              print *,'F1 peak height: t=model, f=user input {t}'
              read(5,*) jf(14)
            endif
              print *,'F1: t=with probability model, f=without   {t}'
              read(5,*) jf(19)
              print *,'F1: t=standard probability, f=with L ',
     &            'condition {t}'
              read(5,*) jf(20)
              print *,'E peak density or foE: t=model, f=user input {t}'
              read(5,*) jf(15)
              print *,'E peak height: t=model, f=user input {t}'
              read(5,*) jf(16)
              print *,'D: t=old model, f=new options {t}'
              read(5,*) jf(24)
          endif
        if(jf(2)) then
              print *,'Te(Ne) model: t=not used, f=correlation is',
     &          ' used. {t}'
              read(5,*) jf(10)
              print *,'Te: t=AEROS/ISIS model, f=InterKosmos model {f}'
              read(5,*) jf(23)
          endif
        if(jf(3)) then
              print *,'Ion comp. model: t=DS78/DY85, f=DS95/TTS05 {f}'
              read(5,*) jf(6)
              print *,'Ni: t=ion composition in %, f=ion densities',
     &             'in cm-3 {t}'
              read(5,*) jf(22)
              endif
            print *,'Equat. Vert. Ion Drift: t=computed, ',
     &            'f=not computed {t}'
            read(5,*) jf(21)
            print *,'Spread-F probability: t=computed, ',
     &            'f=not computed {t}'
            read(5,*) jf(28)
            print *,'Sunspot index: t=from file, f=user input.  {t}'
            read(5,*) jf(17)
            print *,'Ionospheric index: t=from file, f=user input. {t}'
            read(5,*) jf(27)
            print *,'F10.7D Index: t=from file, f=user input {t}'
            read(5,*) jf(25)
            print *,'UT/LT computation: t=no date change, f=ut_lt ',
     &             'subroutine {t}'
            read(5,*) jf(18)
            print *,'Message output unit: t=(UNIT=6), f=(UNIT=12). {t}'
            read(5,*) jf(12)
       endif

c option to enter three additional parameters
c
      if(hx.lt.0) then
        print *,'Three additional output parameters (number:1-48)'
        print *,(pna(j),j=1,10)
        print *,(pna(j),j=11,20)
        print *,(pna(j),j=21,30)
        print *,(pna(j),j=31,40)
        print *,(pna(j),j=41,48)
        print *,'or 0,0,0 for default (F1 probability[40], equ. vert. ',
     &   'ion drift[44], foF2_st/_quiet[45])'
        read(5,*) pad1,pad2,pad3
        if(pad1.eq.0) pad1=40     ! F1 probability
        if(pad2.eq.0) then
             pad2=44              ! equatorial vertical ion drift
             jf(21)=.true.
             endif
        if(pad3.eq.0) pad3=45     ! fof2_storm/foF2_quiet
      endif

c option to enter measured values for NmF2, hmF2, NmF1, hmF1, NmE, hmE,
c N(300), N(400), N(600) if available;
c
          print *,' '
          print *,' '
          print *,' '
          numstp=int((vend-vbeg)/vstp)+1
              if(ivar.eq.1) numstp=1
        if(jf(1)) then
         if(.not.jf(8).or..not.jf(9).or..not.jf(13).or..not.jf(14).or.
     &      .not.jf(15).or..not.jf(16)) then
          var=vbeg
          i=1
2234              if(.not.jf(8)) then
                jf(26)=.false.    ! storm model off, if user input
                        print *,'foF2/Mhz or NmF2/m-3 for ',itext(ivar),
     &                      '=',var
                        read(5,*) oar(1,i)
                    pname(1)='foF2/MHz'
                    if(oar(1,i).gt.30.) pname(1)='NmF2/m-3'
                        endif
              if(.not.jf(9)) then
                        print *,'hmF2/km or M3000F2 for ',itext(ivar),
     &                      '=',var
                        read(5,*) oar(2,i)
                    pname(2)='M(3000)F2'
                    if(oar(2,i).gt.50.) pname(2)='hmF2/km'
                        endif
              if(.not.jf(13)) then
                        print *,'foF1/MHz or NmF1/m-3 for ',itext(ivar),
     &                  '=',var
                        read(5,*) oar(3,i)
                    pname(3)='foF1/MHz'
                    if(oar(3,i).gt.30.) pname(3)='NmF1/m-3'
                        endif
              if(.not.jf(14)) then
                        print *,'hmF1/km for ',itext(ivar),'=',var
                        read(5,*) oar(4,i)
                    pname(4)='hmF1/km'
                        endif
              if(.not.jf(15)) then
                        print *,'foE/MHz or NmE/m-3 for ',itext(ivar),
     &                     '=',var
                        read(5,*) oar(5,i)
                    pname(5)='foE/MHz'
                    if(oar(5,i).gt.30.) pname(5)='NmE/m-3'
                        endif
              if(.not.jf(16)) then
                        print *,'hmE/km for ',itext(ivar),'=',var
                        read(5,*) oar(6,i)
                    pname(6)='hmE/km'
                        endif
           i=i+1
           var=var+vstp
           if(ivar.gt.1.and.var.le.vend) goto 2234
        endif
        endif

c option to enter Ne for Te-Ne relationship
c
        if(jf(2).and..not.jf(10)) then
          var=vbeg
          do 1235 i=1,numstp
                        print *,'Ne(300km),Ne(400km)/m-3',
     &                     ' for ',itext(ivar),'=',var,' [-1 if not]'
                        read(5,*) oar(15,i),oar(16,i)
1235            var=var+vstp
          endif

c option to enter F107D
c
            if(.not.jf(25)) then
                        print *,'User input for F107D:'
                        read(5,*) f107d
                        do i=1,100
                                    oar(41,i)=f107d
                                    enddo
                        endif

c option to enter Rz12 and/or IG12
c
            if(.not.jf(17)) then
                        print *,'User input for Rz12'
                        read(5,*) oar(33,1)
                        do i=2,100
                                    oar(33,i)=oar(33,1)
                                    enddo
                        endif

            if(.not.jf(27)) then
                        print *,'User input for IG12'
                        read(5,*) oar(39,1)
                        do i=2,100
                                    oar(39,i)=oar(39,1)
                                    enddo
                        endif

c end of user input
c

        num1=(vend-vbeg)/vstp+1
        numstp=iabs(num1)
        if(numstp.GT.100) numstp=100

        map='URSI'
        if(jf(5)) map='CCIR'

        bopt='Gulyaeva'
        if(jf(4)) bopt='B0-Table'

        iopt='DS95&TTS04'
        if(jf(6)) iopt='DS78&DY85 '

        sopt='off'
        if(jf(26)) sopt='on '

        topt='TTSA-2000'
        if(jf(23)) topt='IRI-95'

        if(jf(19)) then
              f1opt='Scotto-97 no L'
              if(.not.jf(20)) f1opt='Scotto-97 with L'
        else
              f1opt='IRI-95'
        endif

        hxx=hx
        jmag=jm
        mmdd=imd

c calling IRI subroutine
c

        call iri_web(jmag,jf,xlat,xlon,iy,mmdd,iut,hour,
     &          hxx,htec_max,ivar,vbeg,vend,vstp,outf,oar)

c preparation of results page
c

        write(7,3991) iy,mmdd,hour,timev(iut+1),coorv(jmag+1),xlat,xlon
        if(jf(1)) then
       		if(jf(29)) then
        		if(jf(30)) then
                	write(7,3314)
        		else
                	write(7,3315)
                endif
            else
        		if(jf(30)) then
                	write(7,3316)
        		else
                	write(7,3317)
                endif
            endif
                if(jf(8)) write(7,301) map
                if(jf(9)) write(7,303)
                if(.not.jf(24)) write(7,3081)
                write(7,309) bopt
                write(7,3291) sopt
                write(7,3295) f1opt
                numi=numstp
                if(ivar.eq.1) numi=1
                do j=1,6
                  ij=jfi(j)
                  if(.not.jf(ij)) then
                     write(7,302) pname(j)
                     write(7,402) (oar(j,i),i=1,numi)
                     endif
                  enddo
                endif

        if(jf(2)) write(7,3292) topt
        if(jf(3)) write(7,329) iopt

        if(ivar.eq.1) then
                if(oar(3,1).lt.1.) oar(4,1)=0.
                yp2=0
                if(oar(3,1).gt.0.0) yp2=oar(3,1)/1.e6
                write(7,213) oar(1,1)/1.E6,yp2,oar(5,1)/1.E6
                write(7,214) oar(2,1),oar(4,1),oar(6,1)
        else
                write(7,307)
        endif

        write(7,211) oar(23,1),oar(25,1),oar(27,1)

        if(.not.jf(17)) then
                write(7,223) oar(33,1)
        else
                write(7,212) oar(33,1)
        endif

        if(.not.jf(27)) then
                write(7,2231) oar(39,1)
        else
                write(7,2121) oar(39,1)
        endif


        if(htec_max.gt.50.0) write(7,3914) htec_max
        if(pad1.eq.44.or.pad2.eq.44.or.pad3.eq.44) write(7,3915)
        if(pad1.eq.40.or.pad2.eq.40.or.pad3.eq.40) write(7,4915)
        if(pad1.eq.45.or.pad2.eq.45.or.pad3.eq.45) write(7,4916)
        if(pad1.eq.48.or.pad2.eq.48.or.pad3.eq.48) write(7,4196)

3991    format(///'DATE(yyyy/-ddd or mmdd/hh.h):',I4,'/',I4,'/',F4.1,
     &      A2,2X,A4,' Lat/Long=',F5.1,'/',F6.1/)
3914    format(/'TEC [1.E16 m-2] is obtained by numerical integration',
     &     ' in 1km steps'/'  from 50 to ',f6.1,' km.  t is the',
     &     ' percentage of TEC above the F peak.')
3915    format(/'vdrft: equatorial vertical F-region drift.')
4915    format(/'F1_pb: F1-layer occurrence probability')
4916    format(/'foF2r: foF2_storm/foF2_quiet')
4196    format(/'sp_F: spread-F occurrence probability')
3916    format(/'M3000F2: Propagation factor related to hmF2'/
     &     'B0: bottomside thickness parameter.')
301     format(A4,' maps are used for the F2 peak density (NmF2)')
302     format(A9,' provided by user:')
402     format(7(1PE10.3))
303     format('CCIR maps are used for the F2 peak height (hmF2)')
304     format(A25,'=',F5.1,') provided by user')
307     format(1x/'Solar and magnetic parameter for the 1st profile',
     &          ' point:')
309     format(A8,' option is used for the bottomside thickness ',
     &          'parameter B0')
3081    format('Special D-region output with all options')
329     format(A10,' option is used for ion composition')
3291    format('The foF2 STORM model is turned ',A3)
3292    format(A9,' option is used for the electron temperature')
3293    format(A8,' option is used for the D-region Ne')
3295    format(A16,' option is used for the F1 occurrence probability')
211     format('Solar Zenith Angle/degree',28X,F6.1/
     &          'Dip (Magnetic Inclination)/degree',20X,
     &          F6.2/'Modip (Modified Dip)/degree',26X,F6.2)
212     format('Solar Sunspot Number (12-months running mean) Rz12',
     &          4X,F5.1)
223     format('Solar Sunspot Number (12-months running mean) Rz12',
     &          4X,F5.1,'{user provided input}')
2121    format('Ionospheric-Effective Solar Index IG12',16X,F5.1)
2231    format('Ionospheric-Effective Solar Index IG12',16X,
     &          F5.1,'{user provided input}')
213     format(/'Peak Densities/cm-3: NmF2=',F9.1,'   NmF1=',F9.1,
     &          '   NmE=',F9.1)
214     format('Peak Heights/km:     hmF2=',F9.2,'   hmF1=',F9.2,
     &          '   hmE=',F9.2/)
3314    format('IRI-2001 is used for topside Ne profile')
3315    format('TTS model is used for topside Ne profile')
3316    format('Corrected IRI01 is used for topside Ne profile')
3317    format('NeQuick is used for topside Ne profile')
c
c table head .......................................................
c

        agnr=7          !output unit number
        xtex=imz(ivar)
        if(jmag.gt.0.and.(ivar.eq.2.or.ivar.eq.3)) xtex='GEOM'
        if(iut.gt.0.and.ivar.eq.8) xtex='U.T.'

        PIKTAB=0
        IF(IVAR.NE.1) THEN
                IF(HX.LT.1.0) PIKTAB=1
                IF(HX.LT.0.0) PIKTAB=2
                ENDIF
        if(.not.jf(24)) piktab=3

        IF(PIKTAB.EQ.3) WRITE(7,8199) ITEXT(IVAR),xtex
        IF(PIKTAB.EQ.2) WRITE(7,8194) ITEXT(IVAR),pna(pad1),
     &    pna(pad2),pna(pad3),xtex,uni(pad1),uni(pad2),uni(pad3)
        IF(PIKTAB.EQ.1) WRITE(7,8192) ITEXT(IVAR),xtex
        IF(PIKTAB.EQ.0) WRITE(7,8193) ITEXT(IVAR),xtex

8192  FORMAT(/'-'/2X,A5,6X,'PEAK ALTITUDES IN KM',8X,'PEAK DEN',
     &  'SITIES IN cm-3  TEC top/%'/3X,A4,'    hmF2  hmF1   hmE   ',
     &  'hmD      NmF2   NmF1    NmE    NmD  1E16m-2')
8194  FORMAT(/'-'/2X,A5,15X,'E-VALLEY',3X,'PLAS FREQ/',
     & ' MHz  ',A6,3X,A6,3X,A6/3X,A4,' M3000 B0/km  W/km Depth ',
     &  ' foF2 foF1  foE   ',A3,6X,A3,6X,A3)
8193  FORMAT(/'-'/1X,A5,' ELECTRON DENSITY   TEMPERATURES ',
     &  8X,'ION PERCENTAGES/%',5x,'1E16m-2'/2X,A4,' Ne/cm-3 Ne/NmF2',
     &  ' Tn/K  Ti/K  Te/K  O+  N+  H+ He+ O2+ NO+ Clust TEC t/%')
8199  FORMAT(/'-'/1X,A5,12X,' D-REGION ELECTRON DENSITY IN CM-3'/
     &  2X,A4,7x,'Mechtley',3x,'Friedrich',7X,'Danilov et al. 1995'/
     &  13x,'-Bilitza',4X,'-Torkar',5x,'standard',2x,'major SW',2x,
     &  'strong WA')

        xcor=vbeg

        do 1234 li=1,numstp

c
c special output: peak densities and altitudes
c

      IF(PIKTAB.eq.1) THEN
        if(oar(3,li).lt.1.) oar(4,li)=0.
        iyp1=int(oar(1,li)/1.e6+.5)
        iyp2=0
        if(oar(3,li).gt.0.0) iyp2=int(oar(3,li)/1.e6+.5)
        iyp3=int(oar(5,li)/1.e6+.5)
        iyp4=int(oar(7,li)/1.e6+.5)
            tec=oar(37,li)
        if(tec.gt.0.0) then
            tec=tec/1.e16
            itopp=int(oar(38,li)+.5)
        else
            tec=-1.0
            itopp=-1
        endif
        WRITE(7,3910) XCOR,oar(2,li),oar(4,li),oar(6,li),oar(8,li),
     &    iyp1,iyp2,iyp3,iyp4,tec,itopp
3910    FORMAT(F7.1,2X,4F6.1,1X,I9,3I7,1x,f6.2,i4)
        GOTO 1234
      ENDIF

      IF(PIKTAB.eq.2) THEN
        yp1=sqrt(oar(1,li)/1.24e10)
        yp2=0.0
        if(oar(3,li).gt.0.0) yp2=sqrt(oar(3,li)/1.24e10)
        yp3=sqrt(oar(5,li)/1.24e10)
c        if(pad1.eq.45.and.oar(pad1,li).le.0.0) oar(pad1,li)=-1.
c        if(pad2.eq.45.and.oar(pad2,li).le.0.0) oar(pad2,li)=-1.
c        if(pad3.eq.45.and.oar(pad3,li).le.0.0) oar(pad3,li)=-1.
        WRITE(7,3919) XCOR,oar(36,li),oar(10,li),oar(12,li)-oar(6,li),
     &    oar(11,li)/oar(5,li),yp1,yp2,yp3,oar(pad1,li),oar(pad2,li),
     &    oar(pad3,li)
3919    FORMAT(F7.1,2X,f4.2,1x,f5.1,1x,f5.1,1x,f5.3,F6.2,2F5.2,
     &    1PE9.2,1PE9.2,1PE9.2)
        GOTO 1234
      ENDIF

c special output for D-region options
c
      IF(PIKTAB.eq.3) THEN
        hxx=hx
        if(ivar.ne.1.and.hx.lt.20.) hxx=85.
        ip1=int(outf(1,li)/1.e6+.5)
        ip2=int(outf(13,li)/1.e6+.5)
        if(ivar.eq.1) hxx=xcor
        ihx=int((hxx-60)/5.)+1
        ixne=-1
        ixn1=-1
        ixn2=-1
        if(ihx.gt.0.and.ihx.lt.7) then
              hihx=60+(ihx-1)*5
              ixne=int((outf(14,ihx)+(outf(14,ihx+1)-outf(14,ihx))/
     &                   5*(hxx-hihx))/1.e6+.5)
              ih1=ihx+7
              ixn1=int((outf(14,ih1)+(outf(14,ih1+1)-outf(14,ih1))/
     &                   5*(hxx-hihx))/1.e6+.5)
              ih2=ihx+14
              ixn2=int((outf(14,ih2)+(outf(14,ih2+1)-outf(14,ih2))/
     &                   5*(hxx-hihx))/1.e6+.5)
              endif
        WRITE(7,3819) XCOR,ip1,ip2,ixne,ixn1,ixn2
3819    FORMAT(F7.1,2X,2i10,5X,3I10)
        GOTO 1234
      ENDIF
c
c normal output
c
        if(ivar.eq.1) then
                oar(1,li)=oar(1,1)
                oar(37,li)=oar(37,1)
                oar(38,li)=oar(38,1)
                endif
        jne=int(outf(1,li)/1.e6+.5)
        xner=outf(1,li)/oar(1,li)
        jtn=int(outf(2,li)+.5)
        jti=int(outf(3,li)+.5)
        jte=int(outf(4,li)+.5)
        jio=INT(OUTF(5,li)+.5)
        jih=INT(OUTF(6,li)+.5)
        jihe=INT(OUTF(7,li)+.5)
        jino=INT(OUTF(8,li)+.5)
        jio2=INT(OUTF(9,li)+.5)
        jicl=INT(OUTF(10,li)+.5)
        jin=INT(OUTF(11,li)+.5)
        if(outf(1,li).lt.0) jne=-1
        if(outf(1,li).lt.0) xner=-1.
        if(outf(2,li).lt.0) jtn=-1
        if(outf(3,li).lt.0) jti=-1
        if(outf(4,li).lt.0) jte=-1
        if(outf(5,li).lt.0) jio=-1
        if(outf(6,li).lt.0) jih=-1
        if(outf(7,li).lt.0) jihe=-1
        if(outf(8,li).lt.0) jino=-1
        if(outf(9,li).lt.0) jio2=-1
        if(outf(10,li).lt.0) jicl=-1
        if(outf(11,li).lt.0) jin=-1
        if(tec.gt.0.0) then
            tec=oar(37,li)/1.e16
            itopp=int(oar(38,li)+.5)
        else
            tec=-1.0
            itopp=-1
        endif
        WRITE(7,7117) XCOR,jne,xner,jtn,jti,jte,jio,jin,
     &        jih,jihe,jino,jio2,jicl,tec,itopp
7117    FORMAT(F6.1,I8,1x,F6.3,3I6,7I4,f6.1,i4)

1234    xcor=xcor+vstp

        print *,'Enter 0 to exit or 1 to generate another profile?'
        read(5,*) icontinue
        if (icontinue.gt.0) goto 1

            stop
            end
