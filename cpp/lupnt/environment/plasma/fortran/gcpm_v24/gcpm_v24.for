ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c		Global Core Plasma Model ***Version 2.3***
c		Original Version 1.0, January 1, 2000
c
c     modified by dlg 8/3/2007 to include the polar cap model
c     modified by dlg 1/6/2009 to include seasonal and solar
c           cycle variations in inner plasmasphere from C&A 1992
c     corrected by dlg 6/3/2009 to fix ionosphere-plasmasphere bridge
c           code that was trying to make bridge below the F2 peak.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine gcpm_v24(itime,r,amlt,alatr,akp,outn)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Input parameters:
c
c	itime	integer*4	dimensions=2
c		(1) = yeardoy, e.g. 2001093
c		(2) = miliseconds of day
c	r		real*4		dimension=1
c		geocentric radial distance in Re
c	amlt	real*4		dimension=1
c		solar magnetic local time in hours
c	alatr	real*4		dimension=1
c		solar magnetic latitude in radians
c	akp		real*4		dimension=1
c		planetary Kp index
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Output parameters:
c	outn	real*4		dimensions=4
c			(1) = total electron density in 1/cm^3
c			(2) = total hydrogen density in 1/cm^3
c			(3) = total helium density in 1/cm^3
c			(4) = total oxygen density in 1/cm^3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real*4 r,amlt,alatr,outn(4)
	real*4 clat,al,akp,ne_iri_ps_trough
	real*4 pn(72,10),ps(72,10),ne_iri_cap
	integer*4 itime(2)
	data re/6371.0/,degrad/0.01745329/

c poleward auroral edge
       DATA PN/73.6,73.8,74.0,74.2,74.6,74.8,75.0,75.1,75.3,75.4,75.5,75
     1.6,75.5,75.2,75.1,74.9,74.8,74.7,74.6,74.5,74.5,74.5,74.5,74.5,74.
     25,74.9,75.3,76.1,76.8,77.3,78.0,78.5,78.9,79.2,79.5,79.7,79.9,80.0
     3,80.1,80.2,80.2,80.2,80.2,80.1,80.1,80.0,79.8,79.7,79.4,79.3,79.1,
     478.9,78.6,78.3,78.1,77.8,77.3,76.9,76.4,76.0,75.5,75.1,74.9,74.6,7
     54.3,74.1,73.9,73.8,73.8,73.7,73.6,73.6,
     &72.9,73.0,73.2,73.4,73.7,73.9,74.1,74.4,74.6,74.8,75.0,75
     1.1,75.2,75.3,75.4,75.4,75.4,75.3,75.2,75.1,75.0,74.9,74.8,74.7,74.
     26,74.8,75.0,75.5,76.0,76.5,77.0,77.5,78.0,78.3,78.7,78.9,79.0,79.1
     3,79.2,79.2,79.2,79.1,79.0,78.9,78.8,78.7,78.5,78.2,78.0,77.8,77.6,
     477.4,77.2,77.0,76.8,76.4,76.1,76.0,75.8,75.6,75.3,75.1,74.9,74.7,7
     54.4,74.2,74.0,73.8,73.4,73.3,73.1,73.0,
     &72.3,72.2,72.3,72.6,72.8,73.0,73.2,73.6,73.9,74.1,74.3,74
     1.8,75.0,75.2,75.3,75.4,75.5,75.4,75.3,75.1,75.0,74.9,74.8,74.8,74.
     29,75.1,75.6,76.0,76.4,76.8,77.2,77.4,77.8,77.9,78.0,78.0,78.0,77.9
     3,77.8,77.7,77.6,77.5,77.4,77.2,77.0,76.9,76.8,76.7,76.5,76.4,76.3,
     476.2,76.1,76.1,76.0,75.9,75.8,75.7,75.6,75.4,75.2,75.0,74.8,74.6,7
     54.2,74.0,73.7,73.3,73.0,72.9,72.7,72.4,
     &72.8,72.9,73.0,73.2,73.4,73.6,73.8,74.0,74.2,74.4,74.6,74
     1.8,74.9,74.9,75.0,75.0,75.0,75.0,75.0,75.0,75.0,75.0,75.0,75.1,75.
     21,75.2,75.4,75.6,75.8,76.0,76.1,76.2,76.4,76.6,76.7,76.8,76.9,77.0
     3,77.0,77.0,76.9,76.8,76.6,76.4,76.2,76.1,76.0,75.9,75.8,75.7,75.6,
     475.5,75.5,75.5,75.5,75.5,75.4,75.4,75.3,75.3,75.2,75.1,75.0,74.9,7
     54.6,74.3,74.0,73.9,73.7,73.4,73.2,73.0,
     &73.6,73.8,74.0,74.1,74.2,74.4,74.6,74.7,74.8,74.9,75.0,75
     1.1,75.2,75.2,75.2,75.1,75.0,74.9,74.8,74.8,74.7,74.7,74.7,74.8,74.
     28,74.9,75.0,75.1,75.2,75.3,75.5,75.6,75.7,75.8,75.9,75.9,76.0,76.0
     3,76.0,76.0,76.0,75.9,75.8,75.7,75.6,75.5,75.3,75.2,75.1,75.0,75.0,
     474.9,74.9,74.9,74.9,74.9,75.0,75.0,75.0,75.0,75.0,75.0,74.9,74.8,7
     54.7,74.6,74.3,74.1,74.0,73.9,73.8,73.7,
     &74.2,74.3,74.4,74.6,74.7,74.8,74.9,75.0,75.1,75.2,75.3,75
     1.4,75.5,75.4,75.2,75.1,75.0,74.8,74.6,74.2,74.0,73.9,73.8,73.7,73.
     26,73.7,73.8,73.9,74.0,74.2,74.5,74.9,75.0,75.1,75.2,75.2,75.2,75.2
     3,75.1,75.1,75.0,75.0,75.0,74.9,74.8,74.8,74.7,74.7,74.7,74.7,74.7,
     474.7,74.7,74.7,74.7,74.7,74.7,74.6,74.6,74.6,74.5,74.5,74.5,74.4,7
     54.4,74.3,74.3,74.2,74.2,74.2,74.2,74.2,
     &75.0,75.1,75.2,75.3,75.4,75.5,75.6,75.6,75.7,75.7,75.8,75
     1.8,75.8,75.7,75.5,75.3,75.0,74.8,74.3,74.0,73.5,73.0,72.7,72.5,72.
     25,72.6,72.8,73.0,73.2,73.4,73.6,73.8,74.0,74.1,74.2,74.2,74.2,74.2
     3,74.2,74.2,74.2,74.2,74.2,74.2,74.2,74.1,74.1,74.1,74.1,74.0,74.0,
     474.0,74.0,74.0,74.0,74.0,74.0,74.0,74.0,74.0,74.0,74.0,74.1,74.2,7
     54.3,74.5,74.6,74.7,74.8,74.9,75.0,75.0,
     &75.5,75.5,75.5,75.5,75.5,75.5,75.5,75.4,75.4,75.3,75.2,75
     1.1,75.0,75.0,74.9,74.8,74.7,74.4,74.0,73.5,73.0,72.5,72.0,71.7,71.
     25,71.6,71.8,71.9,72.0,72.2,72.5,72.7,72.8,73.0,73.1,73.2,73.3,73.4
     3,73.6,73.6,73.7,73.7,73.7,73.7,73.7,73.8,73.8,73.8,73.8,73.8,73.8,
     473.8,73.7,73.7,73.6,73.6,73.5,73.5,73.5,73.6,73.8,73.9,74.0,74.1,7
     54.3,74.6,74.8,74.9,75.0,75.1,75.3,75.4,
     &76.0,75.9,75.8,75.7,75.6,75.5,75.4,75.2,75.0,74.9,74.8,74
     1.5,74.3,74.0,73.6,73.1,72.9,72.4,72.0,71.6,71.1,70.8,70.4,70.0,70.
     20,69.9,69.9,70.0,70.2,70.4,70.8,71.0,71.3,71.8,72.0,72.3,72.7,73.0
     3,73.1,73.2,73.4,73.7,73.8,73.8,73.7,73.5,73.2,73.0,72.8,72.6,72.4,
     472.2,72.1,72.1,72.1,72.2,72.3,72.5,72.8,72.9,73.1,73.6,74.0,74.3,7
     54.7,75.0,75.2,75.4,75.8,75.9,76.0,76.0,
     &76.0,75.9,75.8,75.7,75.6,75.5,75.4,75.2,75.0,74.9,74.8,74
     1.5,74.3,74.0,73.6,73.1,72.9,72.4,72.0,71.6,71.1,70.8,70.4,70.0,70.
     20,69.9,69.9,70.0,70.2,70.4,70.8,71.0,71.3,71.8,72.0,72.3,72.7,73.0
     3,73.1,73.2,73.4,73.7,73.8,73.8,73.7,73.5,73.2,73.0,72.8,72.6,72.4,
     472.2,72.1,72.1,72.1,72.2,72.3,72.5,72.8,72.9,73.1,73.6,74.0,74.3,7
     54.7,75.0,75.2,75.4,75.8,75.9,76.0,76.0/
c
c equatorial auroral edge
      data ps/
     &65.5,65.6,65.8,66.0,66.2,66.6,66.8,66.9,67.0,67.2,67.4,67
     1.7,68.0,68.6,69.0,69.6,70.0,70.2,70.6,70.9,71.1,71.3,71.8,72.0,72.
     25,72.8,73.1,73.6,74.0,74.3,74.8,75.0,75.4,75.7,76.0,76.1,76.2,76.4
     3,76.5,76.6,76.7,76.8,76.8,76.8,76.7,76.6,76.4,76.2,76.0,75.8,75.4,
     475.0,74.7,74.2,73.7,73.2,72.8,72.3,71.9,71.4,71.0,70.6,70.0,69.6,6
     59.0,68.4,67.8,67.2,66.8,66.2,65.9,65.7,
     &64.5,64.6,64.8,65.0,65.2,65.5,65.9,66.0,66.2,66.6,67.0,67
     1.2,67.6,68.0,68.5,69.0,69.4,69.9,70.1,70.5,70.8,70.9,71.1,71.3,71.
     27,72.0,72.4,72.9,73.2,73.7,74.0,74.3,74.6,74.8,75.0,75.1,75.2,75.2
     3,75.2,75.1,75.1,75.0,75.0,75.0,74.9,74.8,74.8,74.7,74.6,74.4,74.2,
     473.9,73.6,73.2,72.7,72.1,71.6,70.9,70.3,69.6,69.0,68.6,68.2,67.8,6
     57.3,67.0,66.4,66.0,65.4,65.0,64.8,64.7,
     &63.5,63.6,63.8,64.0,64.2,64.5,64.7,64.8,64.8,64.8,64.9,65
     1.0,65.0,65.0,65.1,65.3,65.8,66.2,66.8,67.4,68.0,68.6,69.1,69.7,70.
     21,70.8,71.2,71.8,72.2,72.8,73.1,73.6,73.8,74.0,74.1,74.2,74.2,74.2
     3,74.2,74.2,74.2,74.1,74.0,74.0,73.9,73.8,73.6,73.3,73.0,72.8,72.3,
     472.0,71.5,71.0,70.4,69.9,69.1,68.5,68.0,67.6,67.1,66.8,66.3,66.0,6
     55.8,65.3,65.0,64.4,64.0,63.9,63.7,63.6,
     &62.5,62.5,62.5,62.6,62.7,62.8,62.9,62.9,63.0,63.1,63.3,63
     1.6,63.8,63.9,64.2,64.8,65.4,66.1,66.9,67.6,68.2,68.9,69.4,70.0,70.
     23,70.8,71.1,71.5,71.9,72.2,72.4,72.7,72.9,73.0,73.1,73.1,73.0,72.9
     3,72.8,72.7,72.6,72.4,72.1,72.0,71.9,71.8,71.7,71.5,71.4,71.2,71.0,
     470.8,70.3,70.0,69.5,69.0,68.3,67.6,66.9,66.0,65.2,64.8,64.3,64.0,6
     53.8,63.7,63.6,63.5,63.3,63.1,63.0,62.8,
     &61.4,61.4,61.3,61.2,61.2,61.1,61.2,61.3,61.4,61.6,61.8,62
     1.0,62.3,63.0,63.7,64.3,65.1,66.0,66.8,67.3,68.1,68.8,69.1,69.6,69.
     28,70.0,70.2,70.5,70.8,71.0,71.3,71.7,71.9,72.0,72.0,72.0,72.0,72.0
     3,72.0,72.0,72.0,72.0,72.0,71.8,71.6,71.4,71.1,70.8,70.3,70.0,69.5,
     469.0,68.4,67.9,67.3,66.7,66.0,65.8,65.2,64.9,64.5,64.1,64.0,63.8,6
     53.5,63.2,63.0,62.8,62.5,62.0,61.9,61.6,
     &60.5,60.5,60.4,60.3,60.2,60.1,60.1,60.2,60.3,60.4,60.8,61
     1.0,61.5,62.0,62.9,63.6,64.4,65.1,65.9,66.4,67.0,67.5,68.0,68.2,68.
     25,68.8,69.0,69.3,69.7,69.9,70.0,70.3,70.6,70.8,70.9,70.9,71.0,71.0
     3,71.0,70.9,70.7,70.5,70.2,70.0,69.6,69.2,68.9,68.6,68.4,68.1,67.9,
     467.7,67.4,67.1,66.8,66.4,66.0,65.5,65.0,64.4,63.6,63.0,62.5,62.0,6
     51.7,61.3,61.0,60.9,60.9,60.9,60.7,60.6,
     &59.6,59.7,59.6,59.6,59.5,59.4,59.4,59.5,59.6,59.8,60.0,60
     1.2,60.8,61.3,62.0,62.8,63.5,64.2,64.9,65.4,66.0,66.4,66.8,67.1,67.
     23,67.5,67.8,68.1,68.4,68.8,69.1,69.4,69.7,69.9,70.0,70.1,70.0,70.0
     3,70.0,69.8,69.5,69.1,68.8,68.3,67.9,67.3,66.8,66.3,66.0,65.8,65.5,
     465.3,65.1,65.0,64.9,64.7,64.2,64.0,63.8,63.5,63.0,62.5,62.1,61.8,6
     51.4,61.0,60.7,60.2,59.9,59.7,59.6,59.6,
     &58.6,58.8,58.9,58.9,58.9,58.9,58.9,58.9,59.0,59.2,59.7,59
     1.9,60.2,60.7,61.2,61.8,62.2,62.8,63.2,63.9,64.5,65.0,65.5,66.0,66.
     24,66.8,67.0,67.3,67.8,68.0,68.2,68.3,68.6,68.8,68.9,68.9,68.8,68.7
     3,68.5,68.3,68.0,67.7,67.2,66.8,66.2,65.7,65.1,64.8,64.2,64.0,63.8,
     463.7,63.6,63.5,63.5,63.4,63.2,63.1,63.0,62.9,62.5,62.2,62.0,61.7,6
     51.2,60.9,60.5,60.0,59.8,59.3,59.0,58.8,
     &57.9,58.0,58.1,58.3,58.6,58.8,59.0,59.3,59.6,59.8,59.9,60
     1.0,60.0,60.4,61.0,61.6,62.1,62.8,63.3,63.8,64.1,64.6,64.9,65.0,65.
     21,65.3,65.6,65.8,66.0,66.2,66.4,66.7,66.8,66.9,67.0,67.0,67.0,67.0
     3,66.8,66.5,66.1,65.7,65.2,64.8,64.1,63.7,63.1,62.8,62.4,62.4,62.5,
     462.7,62.8,62.9,63.0,63.0,63.1,63.0,62.8,62.5,62.0,61.6,61.0,60.5,6
     50.0,59.4,59.0,58.6,58.3,58.2,58.1,58.0,
     &57.9,58.0,58.1,58.3,58.6,58.8,59.0,59.3,59.6,59.8,59.9,60
     1.0,60.0,60.4,61.0,61.6,62.1,62.8,63.3,63.8,64.1,64.6,64.9,65.0,65.
     21,65.3,65.6,65.8,66.0,66.2,66.4,66.7,66.8,66.9,67.0,67.0,67.0,67.0
     3,66.8,66.5,66.1,65.7,65.2,64.8,64.1,63.7,63.1,62.8,62.4,62.4,62.5,
     462.7,62.8,62.9,63.0,63.0,63.1,63.0,62.8,62.5,62.0,61.6,61.0,60.5,6
     50.0,59.4,59.0,58.6,58.3,58.2,58.1,58.0/
c
	data oldmlt/-1.0/,oldkp/-1.0/

	common /irioutput/ rz12,f107
c
c  altrans = the half width in L-shell of over which the transition
c            takes place between the trough and polar cap models.
c            The L-shell at which this transition is centered is
c            obtained from PN, which hold empirical locations for
c            the polarward edge of the auroral zone for several
c            values of Kp and magnetic local time.
	altrans=2.0
c
c  Will execute this section if we need new a new value for the
c  invarient latitude of the polarward edge of the auroral zone.
c  This location is determined as a function of MLT and Kp from
c  the array PN.
	if(oldmlt.ne.amlt .or. oldkp.ne.akp) then
	  oldmlt=amlt
	  oldkp=akp
	  bmlt=amlt*3.0+1.0
	  imlt=bmlt
	  diffmlt=bmlt-float(imlt)
	  if(imlt.gt.72) imlt=1
	  jmlt=imlt+1
	  if(jmlt.gt.72) jmlt=1
c
	  ikp=akp+1.0
	  diffkp=akp-aint(akp)
	  if(ikp.gt.10) ikp=10
	  jkp=ikp+1
	  if(jkp.gt.10) jkp=10
c
	  pn1=(pn(jmlt,ikp)-pn(imlt,ikp))*diffmlt+pn(imlt,ikp)
	  pn2=(pn(jmlt,jkp)-pn(imlt,jkp))*diffmlt+pn(imlt,jkp)
	  ps1=(ps(jmlt,ikp)-ps(imlt,ikp))*diffmlt+ps(imlt,ikp)
	  ps2=(ps(jmlt,jkp)-ps(imlt,jkp))*diffmlt+ps(imlt,jkp)
c
c  Here we must determine the L-shell locations over which the transition
c  will be made from the trough/plasmasphere model to the polar cap model.
	  alatcrits=(ps2-ps1)*diffkp+ps1
	  alatcritn=(pn2-pn1)*diffkp+pn1

	  alcrit=1.0/cos(alatcritn*degrad)**2
	  aurora_mlat=alatcritn
	  tranlow=alcrit-altrans
	  tranhigh=alcrit+altrans
d	  type *,'auroral zone:',imlt,jmlt,ikp,jkp,pn1,pn2,ps1,ps2
	endif

c  We need to obtain the L-shell of the location given, while limiting
c  the maximum L-shell that will be used.  Higher latitudes and L-shells
c  technically extending to infinity will be described by the maximum
c  acceptable L-shell value, without causing problems.
	clat=cos(alatr)**2
	if(clat.lt.1.0e-5) clat=1.0e-5
	al=r/clat
	aheight=(r-1.0)*re
d	  type *,'setup:',tranlow,tranhigh,altrans,alcrit

c  The model contains elements for the plasmasphere and plasmapause, the
c  trough, and the polar cap.  The IRI is used for the ionosphere.
c  The polar cap - trough boundary is chosen to be at L=alcrit.  The function
c  subroutine switchon is used to transition between these regions.  This
c  function is also used to transition between the IRI and the other exterior
c  models.
c
	if(al.lt.tranlow) then
c  execute this section if we are equatorward of the polar cap
	  edensity=ne_iri_ps_trough(r,al,alatr,amlt,akp,itime)
d      type *,'low n:',al,edensity
	else
	  if(al.le.tranhigh) then
c  execute this section if we are in the transition region between
c  the trough and polar cap regions
	    ps_edensity=ne_iri_ps_trough(r,al,alatr,amlt,akp,itime)
	    cap_edensity=ne_iri_cap(r,alatr,amlt,itime)
	    switch=switchon(al,alcrit,altrans)
	    edensity=ps_edensity*(1.0-switch)+
     &		cap_edensity*switch
d      type *,'mid n:',al,ps_edensity,cap_edensity,edensity,switch
	  else
c  execute this section if we are polarward of the trough
	    edensity=ne_iri_cap(r,alatr,amlt,itime)
d      type *,'polar n:',al,edensity
	  end if
	end if
c
c  convert back to #particles per cm^3
	den=edensity/1.0e6

c compute He+ to H+ density ratio in the plasmasphere
	aHeH=10.0**(-1.541-0.176*r+8.557e-3*f107
     &		-1.458e-5*f107*f107)
c Helium concentration drops dramatically with transition to high latitudes
c and open field lines
d	type *,'intermediate aheh',aheh
	aHeH=aHeH*(1.0-switchon(al,alcrit,altrans))
c compute relative O+ density
d	type *,'aheh',aheh
	aheight=(r-1.0)*re
	alphaO=0.995/(1.0+(aheight-350.0)**2/281250.0)**3+0.005
d	type *,'alphaO',alphaO
c compute relative He+ concentration in the plasmasphere
	if(aHeH.ne.0.0) then
	  alphaHeP=(1.0-alphaO)/(1.0+1.0/aHeH)
	  alphaHe=amax1(0.0,alphaHeP*(1.0-exp(-(aheight-400.0)/600.0)))
	else
	  alphaHe=0.0
	end if
c compute densities of H+, He+, and O+
d	type *,'alphaHe',alphaHe
	outn(1)=den
	outn(3)=alphaHe*den
	outn(4)=alphaO*den
	outn(2)=outn(1)-outn(3)-outn(4)

	return
	end
