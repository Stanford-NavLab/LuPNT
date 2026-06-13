c
c  This subroutine models the plasmapause location and slope using a
c  modified Lorentzian function as the basis.  It returns a factor
c  that is used elsewhere to effect inclusion of the plasmapause
c  transition from the inner plasmasphere to the trough.

c  The factor returned varies from a value of 1 well inside a8
c  to a value of 0 well outside a8.
c
c  The rotation of the bulge with variation with Kpmax has been included
c
c  A plasmapause profile is obtained from Carpenter & Anderson [1991] combined
c  with Higel & Wu [1984] and Moldwin et al, [1994].
c  The profile is included in a8.
c
c  An average Kpmax dependence for the plasmapause slope is also included in
c  a9.
c
	function pp_profile(al,amlt,akp,a8)
c
	real pp_profile,al,amlt,a8,a9,factor
	real centroid,akp_old,akp,amlt_old
	data akp_old/-99.0/,amlt_old/-99.0/
d     type *,'into pp_profile:',al,amlt,akp,a8
c
c  Allow for mlt rotation of the buldge with Kpmax
	if((akp.ne.akp_old) .or. (amlt.ne.amlt_old))
     &            call bulge(amlt,akp,a8,a9,centroid)
d	type *,'Recalled subroutine bulge'
	akp_old=akp
	amlt_old=amlt
c
c  compute 10**factor in such a way to avoid floating overflow
	factor=amin1(27.75,2.0*(a9-1.0)*alog10(al/a8))
	pp_profile=(1.0+10.0**factor)**(-a9/(a9-1.0))

d     type *,'leaving pp_profile:',pp_profile,factor,a8,a9,centroid
c
	return
	end
