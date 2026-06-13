c
c	Program runs the subroutine iri_sm to obtain IRI13 densities
c	along an L-shell.
c
c     dlg June 3, 2009  fixed issue with trying to calculate bridge for locations
c                       below the F2 peak along the selected L-shell
c     dlg June 11, 2009 added switchon feature to field aligned bridge function
c                       so that the equatorial density would be reached without
c                       having to use hugh power-law function factors when the
c                       topside fitted power-law was above equatorial density
c                       at the equator.
c
	subroutine iri_ps_bridge(rr,al,alatr,amlt,itime,eq_iri_ps_trough,
     &		transh,rf2,alpha,dno,co,switchh,switchw,istat)
c
	real re,tot_delh,delh,rr,rstart,amltrad
	real rsample1,rsample2
	parameter (re=6371.0,tot_delh=600.0/re,delh=5.0)
	parameter (r=260.0/re+1.0)
c	parameter (r=260.0/re+1.0,delr=delh/re)
	parameter (amltrad=3.1415927/12.0)
c
	real outf(20,100),oarr(50),alatr,along
	real dens,hs,dens_old,rf3
	real delrr,refden
	real amlt,rs,al
	real ro,transh,delhh,alpha,term1
	real eqh,ano,eq_iri_ps_trough,co,fract
	real dent
	real diffr,delr,cosr1,alatr1,cosr2,alatr2
	real ansample1,ansample2,rloc,rpos,rf2
	real diffold,rlocold,an1old,an2old,dl,ahemisphere
	real switchh,switchw
	real*8 dno,dntransh,dtransh,dalpha
	integer*4 itime(2)
	integer istat,icount,iflag

	common /irioutput/ rz12,f107,neiri,nhoiri,nheiri,noiri

d	type *,'entering iri_ps_bridge',rr,al,amlt,itime,eq_iri_ps_trough
c istat must be either 0 or -1 as it is used in an equation later
	istat=0
	dl=al
c !Trevor Garner found error assuming north only, now pass latitude
	ahemisphere=sign(1.0,alatr)
c  get height and densiy of the f2 peak
c	cosrl=amin1(sqrt(r/al),1.0)
c	alatr=acos(cosrl)  !Trevor Garner found error assuming north only, now pass lat
	along=amod((amlt+12.0),24.0)*amltrad
	cosrl=amin1(sqrt(r/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r,itime,outf,oarr)
      r2=oarr(2)/re+1.0
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r2,itime,outf,oarr)
      r2=oarr(2)/re+1.0
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,r2,itime,outf,oarr)

c approximate the F2 peak along the L-shell=al
	rf2=oarr(2)/re+1.0

c  If L-shell is at or below "r", the starting radial distance
c  for searching for the maximum negative slope above the f2 peak,
c  then the L-shell provided is exclusively an ionospheric issue
c  and we need to pass back parameters that will minimize the hassle
c  associated with the rest of the calculation for density, which
c  will necessarily exclude the bridge density anyway.
d     type *,'f2 peak at:',rf2,al

	if(rr.le.rf2) then
	  istat=-1
d	  type *,'No bridge required, istat=-1 ',rs,al
	  return
	endif

c In an effort to reduce the cals to iri2007 the following is used
c to approximate the point of maximum negative slope in the topside
c ionosphere. This has been obtained from a linear fit to this location
c (derived from the search algorithm above) as a function of returned
c rz12 value from IRI2007. That analysis obtained this relationship:
c     ro = (1.05454+-0.000102) + (8.62678e-5+-1.20975e-6)*rz12
      ro = 1.05454 + 8.62678e-5*rz12
d     type *,'fieldaligned_bridge:',rz12,ro,rf2
      if (ro .le. rf2) ro=rf2+0.01

	transh=(ro-1.0)*re
      diffh=1.0
      diffr=diffh/re
      ah1=transh-diffh
      ah2=transh+diffh
	r1=ah1/re+1.0
	r2=ah2/re+1.0
c get the density at the maximum slope height
	cosrl=amin1(sqrt(ro/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	  call iri_sm(alatrl,along,ro,itime,outf,oarr)
	antransh=outf(1,1)

c setup for use of densities and heights of the locations
c on either side of the point of maximum negative slope.
c Since only calculating ro from a fitted function, need to separately
c determine the ionospheric densities above and below to support initial
c calculation of the power law function.
	cosrl=amin1(sqrt(r1/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	call iri_sm(alatrl,along,r1,itime,outf,oarr)
	an1=outf(1,1)
d     type *,'an1: ',alatrl,along,r1,an1,al
	cosrl=amin1(sqrt(r2/al),1.0)
	alatrl=acos(cosrl)*ahemisphere
	call iri_sm(alatrl,along,r2,itime,outf,oarr)
	an2=outf(1,1)
d     type *,'an2: ',alatrl,along,r2,an2,al

	if(al.le.r2) then
	  istat=-1
d	  type *,'No bridge required, istat=-1 ',al,r2
	  return
	endif

	eqh=(al-1.0)*re
d      type *,'bridge=',ah1,ah2,eqh,transh,antransh
d      type *,'      =',an1,an2,eq_iri_ps_trough
      alpha=-alog10(an1/an2)/alog10(ah1/ah2)
      ano=an1*ah1**alpha
d      type *,'intial alpha,ano:',alpha,ano

      an3=ano*eqh**(-alpha)
d      type *,'setup:',an3,eq_iri_ps_trough

c set up use of switch term that will not function by default
      switchh=eqh*2.0
      switchw=eqh/10.0

      if (eq_iri_ps_trough .ge. an3) then
        if(an2.le.eq_iri_ps_trough) then
d     type *,'inverse IRI-eq:'
          alpha=alog10(antransh/eq_iri_ps_trough)/alog10(transh/eqh)
          ano=antransh*transh**alpha
          dno=ano
        else
d      type *,'greater than or equal too'
          co=eq_iri_ps_trough - an3
          alpha=-alog10((an1-co)/(an2-co))/alog10(ah1/ah2)
          ano=(an1-co)*ah1**alpha
          dno=ano
        endif
      else
d      type *,'less than'
c  keep initial alpha and ano values
c  provide switch values that bring the bridge function to the equatorial
c  density at the equator
        switchh=transh+(eqh-transh)/2.0
        switchw=(eqh-transh)/2.0
        dno=ano
        co=0.0
      endif


d	type *,'final=',dno,alpha,co
d	type *,'    =',ah1,ah2,eqh,an1,an2,eq_iri_ps_trough

d	type *,dens_old,dens,delh,hs
d	type *,'leaving iri_ps_bridge',alpha,ano,transh,switchh,switchw,co
	return
	end
