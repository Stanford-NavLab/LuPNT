c
c	Program determines the fit parameters for a power-law function
c	that is matched to the ionosphere at the point of maximum slope
c     and to the slope of the plasmaspheric interior density. This used
c     to be matched to the slope at the point of the maximum density
c     but that sometimes results in the power law function not falling
c     to the interior plasmaspheric densities. That won't work, so now
c     what is being done is to find the point of maximum slope, calculate
c     a power law function that matches the slope at that point, then
c     find where that function has the same slope as the interior plasmaspheric
c     density profile, then use that point and the point of maximum slope
c     to recalculate the power function. That will result in the power law
c     having a slope close to that of the topside ionosphere at the maximum slope
c     but still drop to interior plasmaspheric densities.
c
c     First, then, find the location
c	in the topside ionosphere where the slope is a maximum.  It
c	starts looking at the F2 peak and stops when the slope starts
c	to decrease.  It steps upward with an interval of delh.
c
c	This version of the code only looks along the magnetic equator
c     for producing a smooth equatorial density profile.
c	d.l. gallagher    March 13, 1999.
c
c     d.l.g. November 24, 2007 changed to calculate a power law that insures
c     it falls to the interior plasmaspheric density.
c
	subroutine iri_ps_eq_bridge(al,amlt,itime,transh,alpha,ano,
     &                            am1,b1,x234,psL)
c
	real re,delh,delr,amltrad,tot_delh,r
	parameter (re=6371.0,tot_delh=600.0/re)
	parameter (r=200.0/re+1.0)
	parameter (amltrad=3.1415927/12.0)
c
	real outf(20,100),oarr(50),alatr,along
	real pi,hs,amlt,rs,diff,diff_old,dens_old,dens
	real ano,alpha,transh,x234
      real rz12,f107,neiri,nhoiri,nheiri,noiri
	integer*4 itime(2)
	data pi/3.1415927/

	common /irioutput/ rz12,f107,neiri,nhoiri,nheiri,noiri

d	type *,'entering iri_ps_eq_bridge',amlt,itime,am1,b1,x234

	alatr=0.0
	along=(amlt+12.0)*amltrad
	along=along - (1.0-sign(1.0,(12.0-amlt)))*pi

c  get height and densiy of the f2 peak
d	  type *,'iri_ps_eq_bridge calling iri_sm @r',r,along
	call iri_sm(alatr,along,r,itime,outf,oarr)
	rf2=oarr(2)/re+1.0
d     type *,'f2:',rf2

c In an effort to reduce the cals to iri2007 the following is used
c to approximate the point of maximum negative slope in the topside
c ionosphere. This has been obtained from a linear fit to this location
c (derived from the search algorithm above) as a function of returned
c rz12 value from IRI2007. That analysis obtained this relationship:
c     ro = (1.05454+-0.000102) + (8.62678e-5+-1.20975e-6)*rz12
      ro = 1.05454 + 8.62678e-5*rz12
d     type *,'eq_bridge:',rz12,ro,rf2
c      if (ro .le. rf2) ro=rf2+0.01
      ro=amax1((rf2+0.01),ro)
	transh=(ro-1.0)*re
c	ah1=(rposold-1.0)*re
c	ah2=(rposold+diffr-1.0)*re
      ah1=transh-1.0
      ah2=transh+1.0

c get the density at the maximum slope height
	call iri_sm(alatr,along,ro,itime,outf,oarr)
	dens=outf(1,1)

c since only calculating ro from a fitted function, need to separately
c determine the ionospheric densities above and below to support initial
c calculation of the power law function.
	call iri_sm(alatr,along,(ah1/re+1.0),itime,outf,oarr)
	an1old=outf(1,1)
	call iri_sm(alatr,along,(ah2/re+1.0),itime,outf,oarr)
	an2old=outf(1,1)

d     type *,'power law foder:',an1old,an2old,ah1,ah2

c calculate the initial equatorial power law transition function
	alphao=-log(an1old/an2old)/log(ah1/ah2)
	ano=dens/transh**(-alphao)
d     type *,'bridge init:',alphao,ano,dens,transh
c calculate where this initial power law intersects the plasmasphere profile
      psh=2000.0
      do ii=1,5
d       type *,'psh interation:',ii,psh,ano,alphao,am1,b1
        psh=10.0**((am1*(psh/re+1.0)+b1+x234+6.0-alog10(ano))/(-alphao))
                          ! inner plasmasphere density calculation
      enddo
      psL=psh/re+1.0
c Check to see whether psh is too large, meaning it was continously increasing
c during the 5 iterations. If it is too large, then these two curves did not
c intersect and we have to force a more steep topside ionosphere.
      if (psh .ge.0.5*re) then
c Instead calculate the location along the equator where this power law function
c matches the slope of the interior plasmaspheric density
        psL=1.0 - alphao/am1/alog(10.0)
        psh=(psL-1.0)*re
      endif
d     type *,'psh final:',psh
      psden=10.0**(am1*psL+b1+x234+6.0) !inner plasmasphere density calculation

c new power law value, alpha, needs to match the ionosphere at the point
c of maximum slope and the interior plasmaspheric density where the initial
c power law slope matches the plasmasphere interior density slope
      alpha=-alog10(dens/psden)/alog10(transh/psh)
      ano=dens/transh**(-alpha)
d     type *,'reworked bridge:',psL,psh,psden,am1,b1
d	type *,'leaving iri_ps_eq_bridge',alpha,ano,dens,ro,transh

	return
	end
