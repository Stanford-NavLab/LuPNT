c
c	This subroutine is responsible for deteriming the total electron
c	density as a function of position within the ionosphere and
c	plasmasphere.
c
c     dlg June 3, 2009  fixed a problem with obtaining density below the F2 peak
c
c     dlg June 11, 2009 added switchon feature to field aligned bridge function
c                       so that the equatorial density would be reached without
c                       having to use hugh power-law function factors when the
c                       topside fitted power-law was above equatorial density
c                       at the equator.
c
	function ne_iri_ps_trough(r,al,alatr,amlt,akp,itime)
c
	real pi,amltrad,re,amlt_o,akp_o,al_o
	parameter (pi=3.1415927,amltrad=pi/12.0)
	real r,al,amlt,akp,ne_iri_ps_trough_eq
	real ne_iri_ps_trough,aheight,alatr,transh,alpha
	real outf(20,100),oarr(50),eq_iri_ps_trough
	real ano,co,along,eq_bridge,swtch1,switchon
	real switchh,switchw
	real*8 dno,rf2
	integer*4 itime(2),itime1_o,itime2_o,istat
	data re/6371.0/
	data amlt_o/-99.0/akp_o/-99.0/,al_o/-99.0/
	data itime1_o/-99/,itime2_o/-99/

d     type *,'entering ne_iri_ps_trough with',r,al,alatr
d     type *,'           and ',amlt,akp,itime
c
c  We don't do inside the earth; we're just not that kind of model, humph!
	if(r.le.1.0) then
	  ne_iri_ps_trough=0.0
	  return
	end if
c
c
c  altitude of point of interest
	aheight=(r-1.0)*re
	ahemisphere=sign(1.0,alatr)
c  density at the equator for given L-shell
	eq_iri_ps_trough=ne_iri_ps_trough_eq(al,amlt,akp,itime)
d	type *,'back from ne_iri_ps_trough_eq=',eq_iri_ps_trough
d     type *,'bridge? ',amlt,amlt_o,akp,akp_o,al,al_o
d     type *,'      ? ',itime(1),itime1_o,itime(2),itime2_o
	if(amlt.ne.amlt_o .or. akp.ne.akp_o .or.
     &	itime(1).ne.itime1_o .or. itime(2).ne.itime2_o .or.
     &	abs(al-al_o).gt.1.0e-5 .or. istat.lt.0 .or.
     &  ahemisphere.ne.ahemisphere_o) then
d     type *,'Evaluate field aligned power law function'
d     type *,'showit:',(amlt-amlt_o),(akp-akp_o),(al-al_o)
d     type *,(itime(1)-itime1_o),(itime(2)-itime2_o)
c
c  determine the height power law fit parameters that connect the
c  topside ionosphere to the equatorial density along given L-shell
	  call iri_ps_bridge(r,al,alatr,amlt,itime,eq_iri_ps_trough,
     &		transh,rf2,alpha,dno,co,switchh,switchw,istat)
	  if (istat.lt.0) then
	    along=(amlt-12.0)*amltrad
	    call iri_sm(alatr,along,r,itime,outf,oarr)
	    ne_iri_ps_trough=outf(1,1)
d     type *,'back from iri_ps_bridge with stat<0'
d     type *,'     ',r,alatr,along,ne_iri_ps_trough
	    return
	  endif
	  al_o=al
	  amlt_o=amlt
	  akp_o=akp
	  ahemisphere_o=ahemisphere
	  itime1_o=itime(1)
	  itime2_o=itime(2)
d	type *,'ne_iri_ps_trough transh=',transh,alpha,dno,co,
d    &                  switchh,switchw,istat
d      type *,'   >',r,al,alatr,amlt,itime,eq_iri_ps_trough
	end if

c  compute density as given by the bridge function
c        swtchb=switchon(aheight,switchh,switchw)
         if (aheight .le. (switchh-switchw)) then
           swtchb=0.0
         else
           if (aheight .ge. (switchh+switchw)) then
             swtchb=1.0
           else
             swtchb=(aheight-(switchh-switchw))/(2.0*switchw)
           endif
         endif
	  eq_bridge=(dno*aheight**(-alpha) + co)*(1.0-swtchb)
     &                 + swtchb*eq_iri_ps_trough
d      type *,'bridge ne=',eq_bridge,dno,aheight,alpha,co

c  to call IRI or not to call IRI, that is the question...
	  if (aheight.lt.1000.0) then
	    along=(amlt-12.0)*amltrad
	    call iri_sm(alatr,along,r,itime,outf,oarr)
	    if (r .le. rf2) then
	      ne_iri_ps_trough=outf(1,1)
	    else
	      swtch1=switchon(aheight,transh,5.0)
	      ne_iri_ps_trough=outf(1,1)*(1.0-swtch1) + eq_bridge*swtch1
	      ne_iri_ps_trough=ne_iri_ps_trough*float(istat+1)
     &						+ outf(1,1)*float(iabs(istat))
d	type *,'swtch',outf(1,1),eq_bridge,swtch1,ne_iri_ps_trough
          endif
	  else

	    ne_iri_ps_trough=eq_bridge
	  end if
d     type *,'leaving ne_iri_ps_trough=',ne_iri_ps_trough
d     type *,'                        =',r,al,alatr,amlt,akp
c
d      type *,'ne_iri_ps_trough:',eq_bridge,ne_iri_ps_trough
	return
	end
