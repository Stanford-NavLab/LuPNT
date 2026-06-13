c
c This function determines the trough density in the equatorial
c plane in particles/cm^3.
c
c	written	dlg	07/14/99
c     modified dlg 1/6/2009 to include solar cycle and seasonal variations
c               in inner plasmasphere density from C&A 1992
c
c INPUT:
c	al    = L-shell (real*4)
c	amlt = solar magnetic local time in hours (real*4)
c	akp  = Global Kp index as a real number, e.g 1+=1.3 (real*4)
c
c OUTPUT:
c	ne_eq_trough = trough density in the equatorial plane
c					in cm-3 units (real*4)
c
	function ne_eq_trough(al,amlt,akp)
c
	real ne_eq_trough,al,amlt,akp
	real phitp,antp,damping_time,damping
	real dendamp,dengrow,rz12,f107
	real switch0,switch1,switch2,switch3
	real down_time,del,center,diff,aminden,width,denmin
	real sdel,shift,switchon
	real ne_inner_ps,a6,a7,x234,doy
	real dummy,am1,b1,trough_density
      real zl,oldl,a8,trough
	real outf(20,100),oarr(50),alatr,along
      integer*4 icount,itime(2),iyear,itime1_o,itime2_o
      data a6old/0.0/,a7old/0.0/,a8old/0.0/

	common /irioutput/ rz12,f107

d     type *,'ne_eq_trough caled with',al,amlt,akp

c
c  compute MLT where trough density peaks assuming constant Kp
	phitp=0.145*akp**2-2.63*akp+21.86
c
c  assuming constant density increase of 0.56+-0.08 cm-3h-1 before phitp
c  and constant density decrease of -0.83+-0.15 cm-3h-1 after that and before
c  1 hour MLT, where the decrease is replaced by zero growth.  Growth begins
c  again at 3.5 hours MLT.  A minimum density of 0.18 cm-3 is used.
c  Transitions between growth rates are made over an hour or more, as shown below.
c
c  compute unsmoothed peak density at phitp
	antp=(phitp-3.5)*0.56
c  how long with it take to loose what was gained?
c    minimum loss rate required to get back to 0.18 cm-3 by 2.0 hours MLT
	damping_time=amin1(26.0-phitp,antp/0.83)
c  given the loss rate above as a minimum.
	damping=-antp/damping_time
	down_time=phitp+damping_time
	del=3.5-(down_time-24.0)
	center=3.5-del/2.0
	if(center.lt.0.0) center=24.0+center
	diff=amlt-center
	if(diff.lt.-12.0) diff=24.0+diff
	if(diff.gt.12.0) diff=diff-24.0
	aminden=0.18
c	width=25.0
	width=2.0*del
	denmin=aminden+diff**2/(del*width)
c
c  must deal with computation differently based on current mlt
	dengrow=0.56*(amlt-3.5)+aminden
	sdel=0.4
	shift=0.5
	switch1=switchon(amlt,3.5+shift,sdel)
	switch2=switchon(amlt,phitp,0.5)
	if(amlt.lt.8.0) then
	  dendamp=antp+damping*(amlt+24.0-phitp)
	  switch0=switchon(amlt,down_time-24.0-shift,sdel)
	  geosync_trough= denmin*switch0*(1.0-switch1)
     &		+ dendamp*(1.0-switch0)
     &	           + dengrow*switch1*(1.0-switch2)
d     type *,'lt.8.0:',denmin,switch0,switch1,dendamp,dengrow
	else
	  dendamp=antp+damping*(amlt-phitp)
	  switch3=switchon(amlt,down_time-shift,sdel)
	  geosync_trough= denmin*switch3 + dengrow*switch1*(1.0-switch2)
     &	           + dendamp*switch2*(1.0-switch3)
d     type *,'ge.8.0:',denmin,switch3,dengrow,switch1,switch2,dendamp
	end if
c scale density from trough density at geosynchronous orbit
c to L-shell of interest "l" using a power law of -4.5 in L
c	    =ne_trough[at L=6.6]*al**(-4.5)/6.6**(-4.5)
	trough_density=geosync_trough*al**(-4.5)/2.0514092e-4
d     type *,'trough:',geosync_trough,trough_density,al
	ne_eq_trough=trough_density
c
	return

c
c
c  This program calcuates the density profile for the interior
c  plasmasphere and for the emipirical modeling project.
c
c  The Carpenter and Anderson, 1992 plasmasphere is used only for a saturated plasmasphere.
c  For a more typical plasmaspheric density profile, we will use a profile obtained from
c  Gallagher et al. 1988
c
c  D. Gallagher 4/15/91
c               10/3/94  modified back to Carpenter and Anderson
c                        to represent steady state model
c		    7/10/97  changed to more typical profile based on DE1/RIMS
c					a6=-0.7024 & a7=4.9412
c					these come from RIMS; launch through 82012
c					from the empirical plasmaspheric database
c			1/2/99 find for 3<L-shell<9 and Kp daily sum<168 which is
c					the mean for the DE1/RIMS density samples
c					a6=-0.55 and a7=4.44
c          1/6/2009 include C&A 1992 seasonal and solar cycle variations
c                   in inner plasmasphere.
c          6/11/2009 initialized the value of r prior to calling iri_sm
c
	entry ne_inner_ps(al,amlt,itime,am1,b1,x234)
	data itime1_o/-99/,itime2_o/-99/,pi/3.1415927/
c
c	data a6/-0.49/,a7/4.7/
c	data a6/-0.3145/,a7/3.9043/
c  rims 1981 data set with mean mlt of 20.13 hours gives
c	density=10**(-0.6061*al+4.8596)
c  rims 1981data set with mean mlt of 7.61 hours gives
c	density=10**(-0.758*al+4.975)
c  rims 1981 data set that includes both day and night sides gives
c	density=10**(-0.7*al+4.94)
c	data a6/-0.6061/,a7/4.8596/
c	data a6/-0.55/,a7/4.44/
c  from the JGR [2000] GCPM paper
c      data a6/-0.79/,a7/5.3/
c  from JGR [2000] GCPM paper, but corrected to accept solar cycle and seasonal
c     terms such that the resulting function still fits the DE 1 RIMS data used
c     to obtain the paper results. The mean date on which the paper is based
c     is January 27, 1983. For the inner plasmaspheric density to match the paper
c     values it is necessary to adjust the constant term a7 by 0.092
      data a6/-0.79/,a7/5.208/,re/6371.0/
c
	dummy=amlt		!reserved for future use
c
c calculate seasonal and solar cycle terms
c   Now we're supposed to use the 13 month average sunspot number,
c   but we don't have that available, so we're using the 12 month average value, rz12
c   to get rz12 we'll need to set r to a reasonable value that doesn't cause trouble
      r=1000.0/re + 1.0
c if the date changes, then we need to call IRI to get the current rz12 and
c calculate x234
	if(itime(1).ne.itime1_o .or. itime(2).ne.itime2_o) then

c  call IRI to get the value of rz12 loaded
	  alatr=0.0
	  along=0.0
	  call iri_sm(alatr,along,r,itime,outf,oarr)

        iyear=itime(1)/1000
        doy=float(itime(1) - iyear*1000)
        doy_factor=pi*(doy+9.0)/365.0
        x234=( 0.15*(cos(2.0*doy_factor) - 0.5*cos(4.0*doy_factor))
     &    + (0.00127*rz12-0.0635) ) * exp(-(al-2.0)/1.5)
	  itime1_o=itime(1)
	  itime2_o=itime(2)
c      type *,'season_cycle factor: ',doy,rz12,al,x234
	end if

c calculate the inner plasmaspheric equatorial density
c   same calculation done in check_crossing
c    calculation used also in iri_ps_eq_bridge
	ne_inner_ps=10**(a6*al+a7 + x234)

	am1=a6
	b1=a7
c
	return

c This routine checks to make sure that a8 is not beyond the point where
c the inner plasmasphere and trough density models would cross were there
c no plasmapause. If a8 is beyond that point, then the crossing point is
c substituted for a8.
      entry check_crossing(a8,am1,b1,x234,amlt,akp)

d      type *,'initial crossing=',a8,am1,b1,amlt,akp
c Determine where the inner plasmasphere plus plasmapause profile
c intersects with the trough density profile.
      stepl=0.5
      zl=a8
      a=10**(am1*zl+b1 + x234) !same calculation done in ne_inner_ps
      b=pp_profile(zl,amlt,akp,a8)
      c=geosync_trough*zl**(-4.5)/2.0514092e-4
      diff=a*b - c
d     type *,'crossing:',zl,stepl,diff,a,b,c,geosync_trough,a8,akp,amlt
      icount=0
      do while (abs(stepl).gt.0.05)
        if ((diff.lt.0.0).and.(stepl.gt.0.0) .or.
     &     (diff.gt.0.0).and.(stepl.lt.0.0)) stepl=-stepl/2.0
        zl=zl+stepl
c same calculation done in ne_inner_ps
        diff=10**(am1*zl+b1 + x234)*pp_profile(zl,amlt,akp,a8)
     &     -geosync_trough*zl**(-4.5)/2.0514092e-4
d     type *,'crossing:',zl,stepl,diff
        icount=icount+1
        if (icount.gt.100) then
          temp=pp_profile(zl,amlt,akp,a8)
          type *,'check_crossing is loop-bound:',
     &                am1,b1,zl,x234,amlt,akp,a8,
     &                temp,geosync_trough,stepl
          type *,'STOPPING***********'
          stop
        endif
      enddo

      check_crossing=zl
      return
	end
