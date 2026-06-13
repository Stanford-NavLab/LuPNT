c
c	This subroutine is responsible for deteriming the total electron
c	density as a function of position within the ionosphere and
c	plasmasphere at the magnetic equator.
c
c     modified by dlg 1/6/2009 to send itime to ne_inner_ps routine
c
	function ne_iri_ps_trough_eq(al,amlt,akp,itime)
c
	real pi,amltrad,fact,amlt_o,akp_o
	parameter (pi=3.1415927,amltrad=pi/12.0)
	parameter (fact=0.013/pi)
	real al,amlt,akp,ne_iri_ps_trough_eq,ne_inner_ps
	real outf(20,100),oarr(50),a8,ne_eq_trough,re
	real aheight,pp_factor,pp_profile,ps_inner,am1
	real b1,transh,alpha,ano,x1,x2,y1,y2,am2,b2
	real xintercept,rintercept,ps_bridge,scale,off
	real swtch1,swtch2,switchon,swtch3,alatr,along
	real swtch4,swtch5,ne_eq_trough_den,zl,x234
	integer*4 itime(2),itime1_o,itime2_o
	data re/6371.0/
	data amlt_o/-99.0/akp_o/-99.0/
	data itime1_o/-99/,itime2_o/-99/

d     type *,'entering ne_iri_ps_trough_eq:',al,amlt,akp,itime
c
c  We don't do inside the earth; we're just not that kind of model, humph!
	if(al.le.1.0) then
	  ne_iri_ps_trough_eq=0.0
d       type *,'leaving ne_iri_ps_trough_eq early',ne_iri_ps_trough_eq
	  return
	end if
c
c
c  altitude of point of interest
	aheight=(al-1.0)*re
	pp_factor=pp_profile(al,amlt,akp,a8)
d     type *,'from pp_profile:',pp_factor,al,amlt,akp,a8

	ps_inner=ne_inner_ps(al,amlt,itime,am1,b1,x234)*1.0e6
d     type *,'from ne_inner_ps:',ps_inner,al,amlt,am1,b1,x234

	if(amlt.ne.amlt_o .or. akp.ne.akp_o .or.
     &	itime(1).ne.itime1_o .or. itime(2).ne.itime2_o) then
c
c  determine the height power law fit parameters that connect the
c  topside ionosphere to the plasmasphere
	  call iri_ps_eq_bridge(al,amlt,itime,transh,
     &              alpha,ano,am1,b1,x234,rintercept)
d     type *,'from iri_ps_eq_bridge:',
d    &       r,amlt,itime,transh,alpha,ano,rintercept

	  amlt_o=amlt
	  akp_o=akp
	  itime1_o=itime(1)
	  itime2_o=itime(2)
	end if
d     type *,'values 1:',a8,am1,b1

	  ps_bridge=ano*aheight**(-alpha)
d       type *,'bridge-ps n:',ps_inner,ps_bridge
        off=0.0
        transwidth=0.02
	  swtch2=switchon(al,rintercept+off,transwidth)
	  swtch3=switchon(al,rintercept-off,transwidth)

	    alatr=0.0
	    along=(amlt-12.0)*amltrad
c	    along=along - (1.0-sign(1.0,(12.0-amlt)))*pi
	    call iri_sm(alatr,along,al,itime,outf,oarr)
	    swtch1=switchon(aheight,transh,5.0)
d     type *,'entering ne_eq_trough:',swtch1,aheight,transh
          ne_eq_trough_den=ne_eq_trough(al,amlt,akp)
d     type *,'entering check_crossing:',ne_eq_trough_den
          zl=check_crossing(a8,am1,b1,x234,amlt,akp)
          diff=a8-zl
          offset=(0.0166513
     &            - 0.0450188*diff)*(1.0-switchon(diff,0.3698744,0.05))
c          offset=amax1(0.0,0.02*((zl-a8-0.2)/0.15))
d          type *,'offset:',akp,a8,zl,(a8-zl),offset
	    swtch4=switchon(al,zl+offset,0.3)
	    swtch5=switchon(al,zl-offset,0.3)
	    ne_iri_ps_trough_eq=outf(1,1)*(1.0-swtch1) +
     &  ((ps_bridge*(1.0-swtch2)*swtch1 +
     &	ps_inner*swtch3)*pp_factor)*(1.0-swtch4) +
     &	ne_eq_trough_den*1.0e6*swtch5
d	type *,'** ',swtch1,swtch2,swtch3,swtch4,swtch5
d	type *,al,rintercept,off,aheight,transh,pp_factor
d	type *,outf(1,1),ps_bridge,ps_inner,ne_eq_trough_den
d     type *,ne_iri_ps_trough_eq
c

	return
	end
