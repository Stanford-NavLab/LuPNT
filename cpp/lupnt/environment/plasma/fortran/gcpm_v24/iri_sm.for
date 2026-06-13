c
c	This subroutine is used to call the IRI model from GCPM
c
c  modified to use iri2007 on 30Oct2007 by dlg
c
	subroutine iri_sm(alatr,along,r,itime,outf,oarr)
c
	real alatr,along,outf(20,100),oarr(50)
	real dhour,r,aheight,blongd,ne,nh,nhe,no
	real delh,re,rz12,f107,raddeg
	real pos_sm(3),pos_geo(3),blatr,blatd,blong,rtemp
	integer*4 ddd,jmag,itime(2),yyyy,i
	logical jf(30)

	common /irioutput/ rz12,f107,ne,nh,nhe,no

	data re/6371.0/,raddeg/57.295780/
	data jf/.true.,.true.,.true.,.true.,.false.,.false.,
     &		  .true.,.true.,.true.,.true.,.true.,.false.,
     &		  .true.,.true.,.true.,.true.,.true.,.true.,
     &		  .true.,.true.,.false.,.true.,.false.,.true.,
     &		  .true.,.true.,.true.,.false.,.false.,.false./
	data delh/1.0/
c
	yyyy=itime(1)/1000
	ddd=itime(1)-yyyy*1000
	dhour=float(itime(2))/3600000.0+25.0
c
	aheight=(r-1.0)*re
d	type *,'iri_sm:',aheight,alatr,along,r,itime
	if(aheight.gt.3000.0) then
	  outf(1,1)=0.0
	  oarr(2)=0.0
	  return
	end if

c
	call pol_to_cart(alatr,along,r,pos_sm)
	call sm_to_geo(itime,pos_sm,pos_geo)
	call cart_to_pol(pos_geo,blatr,blong,rtemp)
	blatd=blatr*raddeg
	blongd=blong*raddeg
c
	jmag=0
d	type *,'iri called with:'
d	type *,blatd,blongd,yyyy,-ddd,dhour,aheight
	call iri_sub(jf,jmag,blatd,blongd,yyyy,-ddd,dhour,
     &		aheight,aheight,delh,outf,oarr)
d	type *,'density returned=',outf(1,1),oarr(2)
	  outf(1,1)=amax1(0.0,outf(1,1))

	  rz12=oarr(33)
c	  f107=63.75+rz12*(0.728+rz12*0.00089)	!this is from iris12.for
c this is from Hathaway, March 22, 2000
c	  f107=49.4 + 0.97*rz12 + 17.6*exp(-0.035*rz12)
        f107=oarr(41)
	  ne=outf(1,1)
	  nh=outf(6,1)
	  nhe=outf(7,1)
	  no=outf(5,1)
c
	return
	end
