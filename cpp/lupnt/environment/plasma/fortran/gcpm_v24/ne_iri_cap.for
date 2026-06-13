
c
c	This subroutine provides polar cap densities based on IRI
c     and on the GCPM polar cap model.
c
c     D. Gallagher 08/06/2007
c
c INPUT:
c r	= radial distance (RE)
c alatr	= geomagnetic latitude in radians
c alatd	= geomagnetic latitude in degrees
c amlt	= geomagnetic local time in hours
c
c OUTPUT:
c The density bridge has the form:
c	density = exp(-2.8618*log(r)+refn)
c		where r is height in km
c refn	= this is obtained from the IRI density at 350km altitude
c		and at the latitude, L-shell, magnetic local time,
c		and local time of the point provided by the
c		empirical model user.  It is units of m-3, like the IRI
c		total electron density.
c powern = -2.8618
c refalt = reference altitude where polar cap profile and IRI are made
c          to agree
c
c The mathematical form of the density model comes from Persoon et al. and
c Chandler et al, 1991.  It approximates Alouette/ISIS, DE 1 RIMS, and DE 1
c PWI observations, while allowing for the density to rise and fall with
c the IRI model through the refn parameter.
c
	function ne_iri_cap(r,alatr,amlt,itime)
c
	real*4 alatr,along,rz12,outf(20,100),oarr(50)
	real*4 nb1,ahcrit, r
	real*4 refn,powern
	real*4 ne_iri_cap,aheight,edensity_iri
	integer*4 itime(2)
	logical jf(12)
c
c  ahcrit = center height for transition between iri and polar cap profile
c  overlap = +- range over which transition occurs
	data chrrad/2.6179939e-1/
	data ahcrit/350.0/,overlap/50.0/
	data re/6371.0/
c
	aheight=(r-1.0)*re
c setup input parameters to IRI
	along=(amlt-12.0)*chrrad
c
d     type *,'into cap:',r,alatr,amlt,aheight,along

	if(aheight.lt. (ahcrit-overlap) ) then
c  if a low enough altitude, then will only need to get density from IRI
	  call iri_sm(alatr,along,r,itime,outf,oarr)
	  ne_iri_cap=outf(1,1)
	else
c
	  call iri_sm(alatr,along,(ahcrit+re)/re,itime,outf,oarr)
	  nb1=outf(1,1)
c
c  compute the refn parameter
c    refn=log(nb1)+2.8618*log(350)
	  refn=log(nb1)+16.764
	  powern=-2.8618
c
c aheight = height in km
c bridge density = exp(powern*log(aheight)+refn)
	  ne_iri_cap=exp(powern*log(aheight)+refn)+0.001
d     type *,'cap setup:',nb1,refn,powern,ne_iri_cap
c
	  if(aheight.le. (ahcrit+overlap) ) then
	    call iri_sm(alatr,along,r,itime,outf,oarr)
	    edensity_iri=outf(1,1)
c
	    refhh=ahcrit
	    spred=-0.16
	    widthh=overlap
	    refh2=refhh-spred
	    refh3=refhh+spred
	    switch2=switchon(aheight,refh2,widthh)
	    switch3=switchon(aheight,refh3,widthh)
	    ne_iri_cap=edensity_iri*(1.0-switch3)+ne_iri_cap*switch2
d     type *,'cap height trans:',edensity_iri,ne_iri_cap,switch2,switch3
	  end if
	end if
c
d     type *,'leaving cap:',ne_iri_cap
	return
	end
