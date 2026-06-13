c iritec.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
C
C contains IRITEC, IONCORR subroutines to computed the
C total ionospheric electron content (TEC) and the ionospheric
C correction caused by TEC, respectively.
C
c-----------------------------------------------------------------------
C Corrections
C
C  3/25/96 jmag in IRIT13 as input
C  8/31/97 hu=hr(i+1) i=6 out of bounds condition corrected
C  9/16/98 JF(17) added to input parameters; OUTF(11,50->100)
C  ?/ ?/99 Ne(h) restricted to values smaller than NmF2 for topside
C 11/15/99 JF(20) instead of JF(17)
C 10/16/00 if hr(i) gt hend then hr(i)=hend
C 12/14/00 jf(30),outf(20,100),oarr(50)
C
C Version-mm/dd/yy-Description (person reporting correction)
C 2000.01 05/07/01 current version
c 2000.02 10/28/02 replace TAB/6 blanks, enforce 72/line   [D. Simpson]
c 2000.03 11/08/02 common block1 in iri_tec with F1reg
c 2007.00 05/18/07 Release of IRI-2007
c 2007.02 10/31/08 outf(.,100) -> outf(.,500)
c
C 2012.00 10/05/11 IRI-2012: bottomside B0 B1 model (SHAMDB0D, SHAB1D),
C 2012.00 10/05/11    bottomside Ni model (iriflip.for), auroral foE
C 2012.00 10/05/11    storm model (storme_ap), Te with PF10.7 (elteik),
C 2012.00 10/05/11    oval kp model (auroral_boundary), IGRF-11(igrf.for),
C 2012.00 10/05/11    NRLMSIS00 (cira.for), CGM coordinates, F10.7 daily
C 2012.00 10/05/11    81-day 365-day indices (apf107.dat), ap->kp (ckp),
C 2012.00 10/05/11    array size change jf(50) outf(20,1000), oarr(100).
C
C 2020.00 03/15/23 Inclusion of plasmasphere
C 2020.00 03/15/23 Revised numerical integration and stepsizes
C 2020.01 03/23/23 Revised IRIT13 to IRITEC
C 2020.02 05/10/23 Added JFF(50)
C 2020.02 12/04/23 Deleted subroutine iri_tec, no longer needed
C 2020.03 03/03/24 if(iisect.. moved to after tecb=0         [R. Skantz]
C 2020.03 03/03/24 hmF2,xnmF2,xnorm defined for last segment [R. Skantz]
C
c-----------------------------------------------------------------------
C
C
        subroutine IRITEC_2020(ALATI,ALONG,jmag,jf,iy,md,hour,hbeg,hend,
     &                          hstep,oarr,tecbo,tecto)
c-----------------------------------------------------------------------
c Program for numerical integration of IRI profiles from h=hbeg
C to h=hend.
C
C  INPUT:  ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C          jmag         =0 geographic   =1 geomagnetic coordinates
C          jf(1:50)     =.true./.false. flags; explained in IRISUB.FOR
C          iy,md        date as yyyy and mmdd (or -ddd)
C          hour         decimal hours LT (or UT+25)
c          hbeg,hend    upper and lower integration limits in km
c          hstep        stepsize in km
C
C  OUTPUT: tecbo,tecto  Total Electron Content in m-2 below hmF2 (tecb)
C                       and above hmF2 (tect)
c-----------------------------------------------------------------------

        dimension       outf(20,1000),oarr(100)
        logical         jf(50),jff(50)

c
c turning off computations that are not needed for integration
c

          do 2938 i=1,50
2938         jff(i) = jf(i)
          jff(2)=.false.       ! f=no temperatures (t)
          jff(3)=.false.       ! f=no ion composition (t)
          jff(21)=.false.      ! t=ion drift computed f=not comp.(f)
          jff(28)=.false.	  ! t=spread-F computed f=not comp. (f)
          jff(33)=.false. 	  ! t=auroral boundary   f=off (f)
          jff(34)=.false. 	  ! t=messages on f= off (t)
          jff(35)=.false. 	  ! t=auroral E-storm model on f=off (f)
          jff(47)=.false. 	  ! t=CGM on  f=CGM off (f)

        iisect = int(((hend-hbeg)/hstep)/1000.0)
	hsect = 1000.0 * hstep
	hastep = hstep/2.0
        tect= 0.
        tecb= 0.
	if(iisect.lt.1) goto 2345
c
c calculate IRI densities from hbeg+hstep/2 to hend-hstep/2
c
	do j=1,iisect
          abeg = hbeg + (j-1)* hsect + hastep
          aend = abeg + hsect - hstep
          call IRI_SUB_2020(JFF,JMAG,ALATI,ALONG,IY,MD,HOUR,
     &        abeg,aend,hstep,OUTF,OARR)
	  if(j.lt.2) then
            hmF2 = oarr(2)
            xnmF2 = oarr(1)
            xnorm = xnmF2/1000.
            endif
c
c Numerical integration for each 1000 point segment.
c (xnorm is divided by 1000 to account for heights in km)
c
          hx = abeg + hastep
          do jj=1,1000
            yyy = outf(1,jj) * hstep / xnorm
            if (hx.le.hmF2) then
              tecb = tecb + yyy
            else
              tect = tect + yyy
            endif
            hx = hx + hstep
            enddo
          enddo
c
c Numerical integration for the last segment
c
2345    hlastbeg = hbeg + iisect * hsect
        ilast = int((hend-hlastbeg)/hstep)
        abeg = hlastbeg + hastep
        aend = hlastbeg + ilast * hstep - hastep
        call IRI_SUB_2020(JF,JMAG,ALATI,ALONG,IY,MD,HOUR,
     &        abeg,aend,hstep,OUTF,OARR)
        hmF2 = oarr(2)
        xnmF2 = oarr(1)
        xnorm = xnmF2/1000.
        hx = abeg + hastep
        do jj=1,ilast
          yyy = outf(1,jj) * hstep / xnorm
          if (hx.le.hmF2) then
            tecb = tecb + yyy
          else
            tect = tect + yyy
          endif
          hx = hx + hstep
          enddo

        tecto = tect * xnmF2
        tecbo = tecb * xnmF2

        return
        end
c
c
        real function ioncorr_2020(tec,f)
c-----------------------------------------------------------------------
c computes ionospheric correction IONCORR (in m) for given vertical
c ionospheric electron content TEC (in m-2) and frequency f (in Hz)
c-----------------------------------------------------------------------
        ioncorr_2020 = 40.3 * tec / (f*f)
        return
        end
c
c
